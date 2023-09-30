from read_graphs import read_geom
from molecule_graph import ConnectionGraph
from combine_params import MMParm, MMParmContainer


def read_charges(file: str):
    result = []
    with open(file, 'r') as f:
        for line in f:
            result.append(float(line))
    return result


def main():
    charges, coords = read_geom("struct.xyz")
    cg = ConnectionGraph()
    cg.add_nodes_from_geometry(charges, coords)
    cg.set_bonds()
    atom_indexs = {}
    n_atoms = len(cg.nodes)
    for i in range(n_atoms):
        atom_indexs[i] = cg.nodes[i].Atom + str(i) + "Z"
    mmcharges = read_charges("charges.txt")
    result = []
    result.append("%nprocs=8\n")
    result.append("%mem=2gb\n")
    result.append("#p amber=softfirst nosymm opt geom=connectivity\n")
    result.append("\n")
    result.append("mm opt\n")
    result.append("\n")
    result.append("0 1\n")
    for i in range(n_atoms):
        mm_line = cg.nodes[i].Atom + "-" + atom_indexs[i]
        charge = "-{0:>2.6f}".format(mmcharges[i])
        charge = charge.replace(" ", "")
        mm_line += charge + "    0    "
        mm_line += "{0:>6.3f} {1:>6.3f} {2:>6.3f}".format(cg.nodes[i].x[0], cg.nodes[i].x[1], cg.nodes[i].x[2])
        mm_line += "  L\n"
        result.append(mm_line)
    result.append("\n")
    for i in range(n_atoms):
        con = cg.connections[i]
        clean_con = []
        line = "" + str(i + 1) + " "
        for atom_i in con:
            if atom_i > i:
                clean_con.append(atom_i)
        for atom_i in clean_con:
            line += str(atom_i + 1) + " 1.0 "

        line += "\n"
        result.append(line)

    container = MMParmContainer()
    with open("ff.par", "r") as f:
        for i, line in enumerate(f):
            if "NONBONDED  NBXMOD 5" in line:
                continue
            if "CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0 " in line:
                break
            if "!" in line:
                continue
            len_line = len(line.split())
            if len_line > 3:
                mm_par = MMParm(line)
                container.append(mm_par)
        for line in f:
            if len(line.split()) < 3:
                break
            mm_par = MMParm.read_as_vdw(line)
            container.append(mm_par)
    result.append("\n")
    for b in container.bonds:
        result.append(b.print_gauss_line())
    for a in container.angles:
        result.append(a.print_gauss_line())
    dhs = container.get_combined_dh()
    result.extend(dhs)
    for v in container.vdws:
        result.append(v.print_gauss_line())
    with open('opt.inp', 'w') as f:
        f.writelines(result)
        f.write("\n")
        f.write("\n")


if __name__ == '__main__':
    main()
