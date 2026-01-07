import os
import sys
import argparse

def prRed(skk):  print("\033[91m {}\033[00m" .format(skk), end="")
def prCyan(skk): print("\033[96m {}\033[00m" .format(skk), end="")


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-n", "--no_material_expansion", action="store_true")
    parser.add_argument("input", nargs='*', help="<FLUKA file> [expanded FLUKA file]")
    args = parser.parse_args()

    try:
        inputfile = args.input[0]
    except:
        parser.print_help()
        return 1
    assert inputfile.endswith(".inp"), "Input file should have .inp extension"

    if len(args.input) > 1:
        outputfile = args.input[1]
    else:
        outputfile = inputfile.replace(".inp","_exp.inp")

    lines = []
    file_expandable = False
    with open(inputfile) as ifile:
        for line in ifile:
            # file to be included: open and read it
            if line.startswith("#include"):
                data = line.strip().split()
                include_file = data[1]

                # do not expand materials.inp
                no_material_expansion = include_file.endswith("materials.inp") and args.no_material_expansion
                if os.path.isfile(include_file) and no_material_expansion:
                    lines.append(line)
                elif os.path.isfile(include_file):
                    prCyan("Including: ")
                    print(include_file)
                    with open(include_file, 'r') as fin:
                        lines.append(fin.read())
                        file_expandable = True
                else:
                    lines.append(line)
            # normal line
            else:
                lines.append(line)

    if not file_expandable:
        prRed(f"{inputfile} not expandable\n")
        return 0

    with open(outputfile, 'w') as ofile:
        for line in lines:
            if not line.strip(): continue
            if '\n' in line:
                ofile.write(f"{line}")
            else:
                ofile.write(f"{line}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
