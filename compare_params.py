import parser

from combine_params import MMParmContainer, MMParm
from read_graphs import read_geom, read_charges
from molecule_graph import ConnectionGraph
import sys

ch1, coord1 = read_geom("struct.xyz")
cg = ConnectionGraph()
cg.add_nodes_from_geometry(ch1, coord1)

atom_indexs = []
for i in range(len(cg.nodes)):
    atom_indexs.append(cg.nodes[i].Atom + str(i))

num1 = int(sys.argv[1]) - 1
num2 = int(sys.argv[2]) - 1
num3 = int(sys.argv[3]) - 1
num4 = int(sys.argv[4]) - 1

symbl_to_symbl = {}

for i in range(num2 - num1):
    symbl_to_symbl[atom_indexs[i + num1]] = atom_indexs[i + num3]

params = []


def main():
    container = MMParmContainer()
    container2 = MMParmContainer()
    with open("ff.par", "r") as f:
        for line in f:
            if "NONBONDED  NBXMOD 5" in line:
                break
            if "!" in line:
                continue
            len_line = len(line.split())
            prev_line = line
            is_replaced = False
            if len_line > 3:
                for key, value in symbl_to_symbl.items():
                    if key + " " in prev_line or key + "\t" in prev_line:
                        is_replaced = True
                        line = line.replace(key, value)
                mm_par = MMParm(line)
                if is_replaced:
                    container2.append(mm_par)
                else:
                    container.append(mm_par)

    for param in container2.to_list():
        param2 = container.get_param(param.print_line())
        if param2 != param:
            print("------------------------------")
            print("Don't same\n")
            print("1: " + param.print_line())
            print("2: " + param2.print_line())


if __name__ == "__main__":
    main()


