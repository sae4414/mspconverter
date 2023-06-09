import configparser
from collections import OrderedDict
from os import listdir
from os.path import isfile, join
import re
from pathlib import Path
import sys

class Config:
    """Define constants"""
    AAMass = OrderedDict([('A', 71.037114), ('C', 103.009185), ('D', 115.026943), ('E', 129.042593),
                      ('F', 147.068414), ('G', 57.021464), ('H', 137.058912), ('I', 113.084064),
                      ('K', 128.094963), ('L', 113.084064), ('M', 131.040485), ('N', 114.042927),
                      ('P', 97.052764), ('Q', 128.058578), ('R', 156.101111), ('S', 87.032028),
                      ('T', 101.047679), ('V', 99.068414), ('W', 186.079313), ('Y', 163.0633),
                      ('p', 79.97), ('o', 15.99)
                      ,('h', 0.98), ('c', 57.02), ('a', 42.01),
                      ('r', -17.03), ('y', 43.01), ('d', -18.01), ('t', 26.02)])

    ModMass = {"Oxidation": 15.994915, "CAM": 57.02146, "Carbamidomethyl": 57.02146, "ICAT_light": 227.12,
               "ICAT_heavy": 236.12, "AB_old_ICATd0": 442.20, "AB_old_ICATd8": 450.20, "Acetyl": 42.0106,
               "Deamidation": 0.9840, "Pyro-cmC": -17.026549, "Pyro-glu": -17.026549, "Pyro_glu": -18.010565,
               "Amide": -0.984016, "Phospho": 79.9663, "Methyl": 14.0157, "Carbamyl": 43.00581, "Delta:H[2]C[2]":26.015650}

    H2O = 18.015

    @classmethod
    def load_mods_from_file(cls, filepath):
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(filepath)

        if config.has_section("ModMass"):
            cls.ModMass = dict(config.items("ModMass"))

def GetAAMass(AA):
    return Config.AAMass[AA] + 57.021464 if AA == 'C' else Config.AAMass[AA]

def GetPepMass(pep):
    return sum(Config.AAMass[aa] for aa in pep) + Config.H2O

def msp_to_mgf(msp_file, out_directory):
    with open(msp_file) as f:
        lines = f.readlines()

    count = 0
    spec_counter = 0
    
    found_name = False
    found_mw = False

    prev = 0
    i = 0

    while i < len(lines):
        line = lines[i]
        i += 1
        mods_list = []

        if line.startswith('Name:'):
            name_groups = re.search(r"Name:\s(?P<pep>[a-zA-Z]+)/(?P<charge>\d+)?", line) # only match pep and charge; modifications are in the comment line

            if not name_groups:
                found_name = found_mw = found_num_peaks = False
                continue

            if "iTRAQ" in line or "TMT" in line:
                found_name = found_mw = found_num_peaks = False
                continue

            pep = name_groups['pep']
            l_charge = int(name_groups['charge'])

            mass = GetPepMass(pep)

            found_name = True

        if found_name and line.startswith('Comment:'):
            mass = round(float(re.findall(r"Parent=(\d+.\d+)", line)[0]), 4)
            found_mw = True

            mods_result = re.search(r"Mods=(?P<mod_num>\d+)(?P<mods>\S*)?", line)

            num_mods = int(mods_result['mod_num'])

            if num_mods > 0:
                modstring = mods_result['mods']

                mods_list = []
                if modstring.startswith("("): # assume nist format
                    mods_list = modstring.replace(")", "").split("(")[1:]
                elif modstring.startswith("/"): # assume prosit format 
                    mods_list = modstring.split("/")[1:]
                else:
                    print(f"Cannot parse modifications string: {modstring}; not processing modifications")

                mods_list = [mod.split(',') for mod in mods_list]
                mods_list = [[int(mod[0]), mod[1], mod[2]] for mod in mods_list]
                mods_list = sorted(mods_list, reverse=True, key=lambda x: x[0])
                for mod in mods_list:
                    mod_index = mod[0] + 1
                    assert mod_index <= len(pep)

                    temp_mod = mod[2]

                    if temp_mod in Config.ModMass:
                        mod_mass = Config.ModMass[mod[2]]
                    else:
                        print(f"Found unrecognized PTM: {temp_mod}; ignoring. If this is unexpected, you may need to add the mass of this PTM to your config.ini file")
                        mod_mass = 0

                    if mod_index == 1:
                        pep = str(mod_mass) + pep
                    else:
                        pep = pep[:mod_index] + str(mod_mass) + (pep[mod_index:] if mod_index < len(pep) else "")
                    mass += mod_mass if not found_mw else 0.0

        if found_name and found_mw and line.startswith('Num peaks'):
            title_text = ".".join(msp_file.split('/')[-1].split('.')[:-1])
            mgf_file = title_text + ".mgf"
            with open(join(out_directory, mgf_file), 'a') as f:
                f.write("BEGIN IONS\n")
                f.write("TITLE={}.{}.{}.{}\n".format(title_text, i, i, l_charge))
                f.write("PEPMASS={}\n".format(mass))
                f.write("CHARGE={}+\n".format(l_charge))
                f.write("SEQ={}\n".format(pep))
                spec_counter += 1
                while i < len(lines) and lines[i] != '\n' and not lines[i].startswith("Name"): # prosit output doesn't have newlines between spectra
                    mz_line = lines[i]
                    mz_splits = mz_line.split('\t')
                    moz, intensity = float(mz_splits[0]), float(mz_splits[1])
                    f.write('{} {}\n'.format(moz, intensity))
                    i += 1
                f.write("END IONS\n")
                f.write("\n")

            
            count = count + 1
            pep = 0
            new = int((i / len(lines)) * 100)
            if new > prev + 10:
                print(str(new) + '%')
                prev = new
            found_name = found_mw = False
            del pep, l_charge, num_mods, mods_list, mass

def verify_in_dir(dir_path, ext):
    in_path = Path(dir_path)
    assert in_path.exists() and in_path.is_dir()
    
    files = [join(dir_path, f) for f in listdir(dir_path) if
                 isfile(join(dir_path, f)) and not f.startswith('.') and f.split('.')[-1] == ext]
    assert len(files) > 0
    return files

def main():
    if len(sys.argv) < 3:
        print("The input MSP file and output directory must be provided as command-line arguments")
        sys.exit(1)

    msp_dir = sys.argv[1]

    msp_files = verify_in_dir(msp_dir, "msp")

    out_dir = sys.argv[2]
    out_path = Path(out_dir)
    assert out_path.exists() and out_path.is_dir()

    if len(sys.argv) == 4:
        config_file_path = sys.argv[3]
        assert(Path(config_file_path).exists())
        Config.load_mods_from_file(config_file_path)

    for msp_file in msp_files:
        print(f"Processing {msp_file}\n")        
        msp_to_mgf(msp_file, out_dir) 

if __name__ == "__main__":
    main()