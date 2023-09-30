from read_graphs import read_geom
from molecule_graph import ConnectionGraph
from combine_params import MMParm, MMParmContainer
from files_generatos import generate_par, generate_psf
from pathlib import Path
import numpy as np
import sys

atom_allias = {
    "X": "Pt",
    "T": "Cl"}


def read_pdb(file: str):
    charges = []
    coords = []
    with open(file, 'r') as f:
        for line in f:
            if "ATOM" not in line:
                continue
            coords.append([float(v) for v in line[30:56].split()])
            charge = line.split()[2][0]
            if charge in atom_allias.keys():
                charge = atom_allias[charge]
            charges.append(charge)
    return charges, coords


def check_atom_on_border(coord: list, percent: float, max_x, min_x, max_y, min_y, max_z, min_z) -> bool:
    bord_x1 = (max_x - min_x) * ((1 - percent) / 2) + min_x
    bord_x2 = (max_x - min_x) * (1 - (1 - percent) / 2) + min_x
    if min_x - 1 < coord[0] < bord_x1 or bord_x2 < coord[0] < max_x + 1:
        return True
    bord_y1 = (max_y - min_y) * ((1 - percent) / 2) + min_y
    bord_y2 = (max_y - min_y) * (1 - (1 - percent) / 2) + min_y
    if min_y - 1 < coord[1] < bord_y1 or bord_y2 < coord[1] < max_y + 1:
        return True
    bord_z1 = (max_z - min_z) * ((1 - percent) / 2) + min_z
    bord_z2 = (max_z - min_z) * (1 - (1 - percent) / 2) + min_z
    if min_z - 1 < coord[2] < bord_z1 or bord_z2 < coord[2] < max_z + 1:
        return True
    return False


def check_atom_in_border(coord: list, l_size: float, max_x, min_x, max_y, min_y, max_z, min_z) -> bool:
    lx = (max_x - min_x)

    bord_x1 = lx / 2 + min_x - l_size / 2
    bord_x2 = lx / 2 + min_x + l_size / 2
    if not(bord_x2 > coord[0] > bord_x1):
        return False
    ly = (max_y - min_y)
    bord_y1 = ly / 2 + min_y - l_size / 2
    bord_y2 = ly / 2 + min_y + l_size / 2
    if not(bord_y2 > coord[1] > bord_y1):
        return False
    lz = max_z - min_z
    bord_z1 = lz / 2 + min_z - l_size / 2
    bord_z2 = lz / 2 + min_z + l_size / 2
    if not(bord_z2 > coord[2] > bord_z1):
        return False
    return True

try:
    int("1")
except:
    pass


def read_charges_from_psf(file: str):
    result = []
    with open(file, 'r') as f:
        for line in f:
            if "!NATOM" in line:
                break
        for line in f:
            if len(line.split()) < 3:
                break
            result.append(float(line.split()[6]))
    return result


def read_indexs_from_psf(file: str):
    result = []
    with open(file, 'r') as f:
        for line in f:
            if "!NATOM" in line:
                break
        for line in f:
            if len(line.split()) < 3:
                break
            result.append(line.split()[5])
    return result


def read_pdb_indexes(file: str):
    indexs = []
    category = []
    with open(file, 'r') as f:
        for line in f:
            if "ATOM" not in line:
                continue
            indexs.append(line.split()[2])
            category.append(line.split()[3])
    return indexs, category


def main():
    if len(sys.argv) != 4:
        print("You need to set last qm number, all solvent box size, inner solvent box size")
        print("Example 139 34 0.8")
        exit()

    last_qm_atom = int(sys.argv[1])
    percent_inner = float(sys.argv[3])


    if len(sys.argv) > 2:
        boorder_cutoff = float(sys.argv[2])
        charges, coords = read_pdb("mol.pdb")
        max_x = max([v[0] for v in coords])
        min_x = min([v[0] for v in coords])
        max_y = max([v[1] for v in coords])
        min_y = min([v[1] for v in coords])
        max_z = max([v[2] for v in coords])
        min_z = min([v[2] for v in coords])
        leave_atoms = []
        all_atoms = []
        n_atoms = len(charges)
        for i in range(n_atoms):
            if i < last_qm_atom:
                continue
            if charges[i] == "O":
                if check_atom_in_border(coords[i], boorder_cutoff, max_x, min_x, max_y, min_y, max_z, min_z):
                    leave_atoms.append(i)
        for i in range(last_qm_atom):
            all_atoms.append(i)
        for i in leave_atoms:
            all_atoms.append(i)
            all_atoms.append(i + 1)
            all_atoms.append(i + 2)
            # all_atoms.append(i + 3)
            # all_atoms.append(i + 4)
            # all_atoms.append(i + 5)
        coords1 = []
        charges1 = []
        with open("cut.xyz", "w") as f:
            for i in all_atoms:
                line = charges[i] + " "
                line += str(coords[i][0]) + " "
                line += str(coords[i][1]) + " "
                line += str(coords[i][2]) + " \n"
                f.write(line)
        for i in all_atoms:
            charges1.append(charges[i])
            coords1.append(coords[i])
        charges, coords = charges1, coords1
    else:
        charges, coords = read_pdb("mol.pdb")
    max_x = max([v[0] for v in coords])
    min_x = min([v[0] for v in coords])
    max_y = max([v[1] for v in coords])
    min_y = min([v[1] for v in coords])
    max_z = max([v[2] for v in coords])
    min_z = min([v[2] for v in coords])
    n_atoms = len(charges)
    choosen_ozygens = []
    for i in range(n_atoms):
        if i < last_qm_atom:
            continue
        if charges[i] == "O":
            if check_atom_on_border(coords[i], percent_inner, max_x, min_x, max_y, min_y, max_z, min_z):
                choosen_ozygens.append(i)
    choosen_atoms = []
    for i in choosen_ozygens:
        choosen_atoms.append(i)
        choosen_atoms.append(i + 1)
        choosen_atoms.append(i + 2)
        # choosen_atoms.append(i + 3)
        # choosen_atoms.append(i + 4)
        # choosen_atoms.append(i + 5)
    qm_part = [i for i in range(last_qm_atom)]
    with open("qm.xyz", "w") as f:
        for i in range(last_qm_atom):
            line = charges[i] + " "
            line += str(coords[i][0]) + " "
            line += str(coords[i][1]) + " "
            line += str(coords[i][2]) + " \n"
            f.write(line)
    with open("inner_water.xyz", 'w') as f:
        for i in range(n_atoms):
            if i < last_qm_atom:
                continue
            if i not in choosen_atoms:
                line = charges[i] + " "
                line += str(coords[i][0]) + " "
                line += str(coords[i][1]) + " "
                line += str(coords[i][2]) + " \n"
                f.write(line)
    with open("outer_water.xyz", 'w') as f:
        for i in choosen_atoms:
            line = charges[i] + " "
            line += str(coords[i][0]) + " "
            line += str(coords[i][1]) + " "
            line += str(coords[i][2]) + " \n"
            f.write(line)
    cg = ConnectionGraph()
    if not Path("mol.psf").is_file():
        cg = ConnectionGraph()
        cg.add_nodes_from_geometry(charges, np.array(coords).reshape(-1))
        charges = []
        with open("charges.txt", 'r') as f:
            for line in f:
                charges.append(float(line))
        atom_indexs, category = read_pdb_indexes("mol.pdb")
        print("generate full bonds")
        cg.set_bonds()
        print("end full bonds genratetion")
        qm_cg = ConnectionGraph()
        for i in range(last_qm_atom):
            qm_cg.add_nodes(cg.nodes[i].Atom, cg.nodes[i].x, i)
        qm_cg.set_bonds()
        generate_psf("mol.psf", cg, charges=charges, is_three_ind=True, force_indexes=atom_indexs)
        generate_par("ff.par", qm_cg, to_zero=False, is_three_ind=True, force_indexs=atom_indexs)
        print("You need to correct ff.par and restart script")
        exit()
    cg.read_from_psf("mol.psf")
    n_atoms = len(charges)
    atom_indexs = read_indexs_from_psf("mol.psf")
    mmcharges = read_charges_from_psf("mol.psf")

    result = []
    result.append("%nprocs=8\n")
    result.append("%mem=2gb\n")
    result.append("#p amber=softfirst nosymm opt geom=connectivity\n")
    result.append("\n")
    result.append("mm opt\n")
    result.append("\n")
    result.append("0 1\n")
    for i in range(n_atoms):
        if cg.nodes[i].Atom in atom_allias.keys():
            mm_line = atom_allias[cg.nodes[i].Atom] + "-" + atom_indexs[i]
        else:
            mm_line = cg.nodes[i].Atom + "-" + atom_indexs[i]
        charge = "-{0:>2.6f}".format(mmcharges[i])
        charge = charge.replace(" ", "")
        if i in qm_part:
            mm_line += charge + "   -1    "
        elif i in choosen_atoms:
            mm_line += charge + "   -1    "
        else:
            mm_line += charge + "    0    "
        mm_line += "{0:>6.3f} {1:>6.3f} {2:>6.3f}".format(coords[i][0], coords[i][1], coords[i][2])
        if i in qm_part:
            mm_line += "  H\n"
        else:
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
        result.append(b.print_gauss_line(with_out_Z=True))
    for a in container.angles:
        result.append(a.print_gauss_line(with_out_Z=True))
    dhs = container.get_combined_dh(with_out_Z=True)
    result.extend(dhs)
    for v in container.vdws:
        result.append(v.print_gauss_line(with_out_Z=True))
    with open('opt.inp', 'w') as f:
        f.writelines(result)
        f.write("\n")


if __name__ == '__main__':
    main()
