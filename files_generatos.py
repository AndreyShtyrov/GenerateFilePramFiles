from molecule_graph import ConnectionGraph
from combine_params import MMParmContainer, MMParm
from constants import CODE_TO_CHARGE, CHARGE_TO_MASS, VDW_A_EC, VDW_A_ER, alphabet, atom_alias
import numpy as np


def generate_un_three_symbol_par(atom: str, i: int, rep_param: int) -> str:
    if atom in atom_alias.keys():
        res = atom_alias[atom]
    else:
        res = atom
    i = i - (i // rep_param) * rep_param
    v = i // 36
    res += alphabet[v]
    v = i % 36
    res += alphabet[v]
    return res


def generate_psf(file: str, cg: ConnectionGraph, charges: list =None, is_three_ind=False,
                 imprs=None, rep_param=None, force_indexes=None):
    with open(file, 'w') as f:
        f.write("PSF\n")
        f.write("\n")
        f.write('       3 !NTITLE\n')
        f.write(' REMARKS original generated structure x-plor psf file\n')
        f.write(' REMARKS topology mol.top\n')
        f.write(' REMARKS segment A { first NONE; last NONE; auto angles dihedrals }\n')
        f.write('\n')
        f.write("{0:>8}".format(len(cg.nodes)) + " !NATOM\n")
        n_atoms = len(cg.nodes)
        for i in range(len(cg.nodes)):
            line = "{0:>8}".format(i + 1) + " A    1    UNL"
            atom = cg.nodes[i].Atom
            if cg.nodes[i].Atom in atom_alias.keys():
                atom = atom_alias[cg.nodes[i].Atom]
            if force_indexes is not None:
                line += "{0:>5}".format(cg.nodes[i].Atom) + "{0:>5}".format(force_indexes[i])
            else:
                if is_three_ind:
                    atom = generate_un_three_symbol_par(cg.nodes[i].Atom, i, n_atoms)
                else:
                    atom = cg.nodes[i].Atom + str(i)

                if is_three_ind:
                    if rep_param is not None:
                        atom_ind = generate_un_three_symbol_par(cg.nodes[i].Atom, i, n_atoms)
                    else:
                        atom_ind = generate_un_three_symbol_par(cg.nodes[i].Atom, i, rep_param)
                    line += "{0:>5}".format(atom) + "{0:>5}".format(atom_ind)
                else:
                    line += "{0:>5}".format(atom) + "{0:>5}".format(cg.nodes[i].Atom + str(i))
            if charges is None:
                line += "     0.000000"
            else:
                line += "  {0:>9.6f}".format(charges[i])
            line += "{0:>14.4f}".format(CHARGE_TO_MASS[CODE_TO_CHARGE[cg.nodes[i].Atom]])
            line += "{0:>12}".format(0) + "\n"
            f.write(line)
        f.write("\n")
        i_vectors = cg.get_all_pairs()

        acc = 0
        n = len(i_vectors)
        f.write("{0:>8}".format(n) + " !NBOND: bonds\n")
        while True:
            if acc >= n:
                break
            if acc + 4 > n:
                c_bonds = i_vectors[acc:]
            else:
                c_bonds = i_vectors[acc: acc + 4]
            acc += 4
            line = ""
            for b in c_bonds:
                line += "{0:>8}".format(b[0]) + "{0:>8}".format(b[1])
            line += "\n"
            f.write(line)
        f.write("\n")
        acc = 0
        i_vectors = cg.get_all_triplets()
        n = len(i_vectors)
        f.write("{0:>8}".format(n) + " !NTHETA: angles\n")
        while True:
            if acc >= n:
                break
            if acc + 3 > n:
                c_bonds = i_vectors[acc:]
            else:
                c_bonds = i_vectors[acc: acc + 3]
            acc += 3
            line = ""
            for b in c_bonds:
                line += "{0:>8}".format(b[0]) + "{0:>8}".format(b[1]) + "{0:>8}".format(b[2])
            line += "\n"
            f.write(line)
        acc = 0
        i_vectors = cg.get_all_quartet()
        f.write("\n")
        f.write("{0:>8}".format(len(i_vectors)) + " !NPHI: dihedrals\n")
        n = len(i_vectors)
        while True:
            if acc >= n:
                break
            if acc + 2 > n:
                c_bonds = i_vectors[acc:]
            else:
                c_bonds = i_vectors[acc: acc + 2]
            acc += 2
            line = ""
            for b in c_bonds:
                line += "{0:>8}{1:>8}{2:>8}{3:>8}".format(b[0], b[1], b[2], b[3])
            line += "\n"
            f.write(line)
        f.write("\n")
        if imprs is not None:
            expand_imprs = []
            if rep_param is not None:
                for impr in imprs:
                    i = 0
                    v_impr = np.array(impr)
                    while True:
                        if max(v_impr + i * rep_param) < n_atoms:
                            expand_imprs.append(v_impr + i * rep_param)
                        else:
                            break
                        i += 1
            if rep_param is not None and len(expand_imprs) > 0:
                imprs = expand_imprs
        if imprs is None:
            f.write("{0:>8}".format(0) + " !NIMPHI: impropers\n")
        else:
            f.write("{0:>8}".format(len(imprs)) + " !NIMPHI: impropers\n")
            n = len(imprs)
            acc = 0
            while True:
                if acc >= n:
                    break
                if acc + 2 > n:
                    c_bonds = imprs[acc:]
                else:
                    c_bonds = imprs[acc: acc + 2]
                acc += 2
                line = ""
                for b in c_bonds:
                    line += "{0:>8}{1:>8}{2:>8}{3:>8}".format(b[0], b[1], b[2], b[3])
                line += "\n"
                f.write(line)
        f.write("\n")
        f.write("{0:>8}".format(0) + " !NDON: donors\n")
        f.write("\n")
        f.write("{0:>8}".format(0) + " !NACC: acceptors\n")
        f.write("\n")
        f.write("{0:>8}".format(0) + " !NNB\n")
        f.write("\n")
        acc = 0
        n = len(cg.nodes)
        while True:
            if acc >= n:
                break
            if acc + 8 > n:
                c_n = n - acc
            else:
                c_n = 8
            acc += 8
            line = ""
            for _ in range(c_n):
                line += "{0:>8}".format(0)
            f.write(line + "\n")
        f.write("\n")
        f.write("       1       0 !NGRP\n")
        f.write("       0       0       0\n")
        f.write("\n")


def generate_par(file: str, cg: ConnectionGraph, to_zero: bool, is_three_ind=False, imprs=None,
                 rep_param=None, force_indexs=None):
    atom_indexs = []
    st_indexs = []
    n_atoms = len(cg.nodes)
    if is_three_ind:
        for i in range(len(cg.nodes)):
            if rep_param is None:
                atom_indexs.append(generate_un_three_symbol_par(cg.nodes[i].Atom, i, n_atoms))
            else:
                atom_indexs.append(generate_un_three_symbol_par(cg.nodes[i].Atom, i, rep_param))
        for i in range(len(cg.nodes)):
            st_indexs.append(cg.nodes[i].Atom + str(i))
    else:
        for i in range(len(cg.nodes)):
            atom_indexs.append(cg.nodes[i].Atom + str(i))
    if force_indexs is not None:
        atom_indexs = force_indexs
    with open(file, 'w') as f:
        f.write("BONDS\n")
        bonds = cg.get_all_pairs()
        for b in bonds:
            line = "{0:>5}".format(atom_indexs[b[0] - 1]) + " "
            line += "{0:>5}".format(atom_indexs[b[1] - 1]) + " "
            if to_zero:
                line += "{0:>7.3f}".format(0.0) + " " + "{0:>6.3f}".format(0.0)
            else:
                line += "{0:>7.3f}".format(110.0) + " " + "{0:>6.3f}".format(1.4)
            line += "\n"
            f.write(line)
        f.write("\n")
        f.write("THETAS\n")
        bonds = cg.get_all_triplets()
        for b in bonds:
            line = "{0:>5}".format(atom_indexs[b[0] - 1]) + " "
            line += "{0:>5}".format(atom_indexs[b[1] - 1]) + " "
            line += "{0:>5}".format(atom_indexs[b[2] - 1]) + " "
            if to_zero:
                line += "{0:>7.3f}".format(0.0) + " " + "{0:>6.3f}".format(0.0)
            else:
                line += "{0:>7.3f}".format(50.0) + " " + "{0:>6.3f}".format(110.0)
            line += "\n"
            f.write(line)
        f.write("\n")
        f.write("PHI\n")
        bonds = cg.get_all_quartet()
        for b in bonds:
            for i in range(4):
                line = "{0:>5}".format(atom_indexs[b[3] - 1]) + " "
                line += "{0:>5}".format(atom_indexs[b[2] - 1]) + " "
                line += "{0:>5}".format(atom_indexs[b[1] - 1]) + " "
                line += "{0:>5}".format(atom_indexs[b[0] - 1]) + " "
                if i % 2 == 0:
                    if to_zero:
                        line += "{0:>7.4f}".format(1.0) + " " + str(i + 1) + " " + "{0:>6.4f}".format(0.0)
                    else:
                        line += "{0:>7.4f}".format(1.0) + " " + str(i + 1) + " " + "{0:>6.4f}".format(0.0)
                else:
                    line += "{0:>7.4f}".format(0.0) + " " + str(i + 1) + " " + "{0:>6.4f}".format(180.0)
                line += "\n"
                f.write(line)
        f.write("\n")
        if imprs is not None:
            f.write("IMPHI\n")
            for impr in imprs:
                line = "{0:>5}".format(atom_indexs[impr[0] - 1]) + " "
                line += "{0:>5}".format(atom_indexs[impr[1] - 1]) + " "
                line += "{0:>5}".format(atom_indexs[impr[2] - 1]) + " "
                line += "{0:>5}".format(atom_indexs[impr[3] - 1]) + " "
                line += " 7.50000 2 180.00000\n"
                f.write(line)
            f.write("\n")

        f.write("NONBONDED  NBXMOD 5  GROUP SWITCH CDIEL -\n")
        f.write("     CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4\n")
        for i in range(len(cg.nodes)):
            line = atom_indexs[i]
            line += " 0.00 "
            line += "{0:>9.6f}".format(VDW_A_EC[cg.nodes[i].Atom]) + " " + "{0:>9.6f}".format(VDW_A_ER[cg.nodes[i].Atom])
            line += " 0.00 "
            line += "{0:>9.6f}".format(VDW_A_EC[cg.nodes[i].Atom] / 2) + " " + "{0:>9.6f}".format(VDW_A_ER[cg.nodes[i].Atom])
            line += "\n"
            f.write(line)
    if is_three_ind:
        with open('falias.txt', 'w') as f:
            for i in range(len(st_indexs)):
                f.write(st_indexs[i] + "   " + atom_indexs[i] + "\n")


def generate_par_from_MMContainer(file: str, mm_container: MMParmContainer, imprs=None, atom_indexes=None):
    with open(file, 'w') as f:
        with open(file, 'w') as f:
            f.write("BONDS\n")
            for bond in mm_container.bonds:
                f.write(bond.print_line())
            f.write("\n")
            f.write("THETAS\n")
            for angle in mm_container.angles:
                f.write(angle.print_line())
            f.write("\n")
            f.write("PHI\n")
            for dh in mm_container.dih:
                f.write(dh.print_line())
            if imprs is not None and atom_indexes is not None:
                f.write("IMPHI\n")
                for impr in imprs:
                    line = "{0:>5}".format(atom_indexes[impr[0] - 1]) + " "
                    line += "{0:>5}".format(atom_indexes[impr[1] - 1]) + " "
                    line += "{0:>5}".format(atom_indexes[impr[2] - 1]) + " "
                    line += "{0:>5}".format(atom_indexes[impr[3] - 1]) + " "
                    line += " 7.50000 2 180.00000\n"
                    f.write(line)
                f.write("\n")
            f.write("NONBONDED  NBXMOD 5  GROUP SWITCH CDIEL -\n")
            f.write("     CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4\n")
            for vdw in mm_container.vdws:
                f.write(vdw.print_line())


def generate_top(file: str, cg: ConnectionGraph, charges: list= None, is_three_ind=False, rep_param=None):
    atom_indexs = []
    n_atoms = len(cg.nodes)
    if is_three_ind:
        for i in range(len(cg.nodes)):
            if rep_param is None:
                atom_indexs.append(generate_un_three_symbol_par(cg.nodes[i].Atom, i, n_atoms))
            else:
                atom_indexs.append(generate_un_three_symbol_par(cg.nodes[i].Atom, i, rep_param))
    else:
        for i in range(len(cg.nodes)):
            atom_indexs.append(cg.nodes[i].Atom + str(i))
    with open(file, 'w') as f:
        for i in range(len(cg.nodes)):
            line = "MASS" + " " + str(i + 1) + " " + atom_indexs[i]
            line += "{0:>10.4f}".format(CHARGE_TO_MASS[CODE_TO_CHARGE[cg.nodes[i].Atom]]) + " "
            if cg.nodes[i].Atom in atom_alias.keys():
                line += atom_alias[cg.nodes[i].Atom]
            else:
                line += cg.nodes[i].Atom
            f.write(line + "\n")
        f.write("AUTO ANGLES DIHE\n")
        if charges is None:
            f.write("RESI   UNL 0.000\n")
        else:
            sum_charge = sum(charges)
            f.write("RESI   UNL {0:>4.3f}\n".format(sum_charge))
        for i in range(len(cg.nodes)):
            line = "ATOM " + "{0:>5} {0:>5}".format(atom_indexs[i], atom_indexs[i])
            if charges is None:
                line += " 0.0000\n"
            else:
                line += " {0:>6.4f}\n".format(charges[i])
            f.write(line)
        bonds = cg.get_all_pairs()
        for bond in bonds:
            line = "BOND " + atom_indexs[bond[1] - 1] + " " + atom_indexs[bond[0] - 1]
            f.write(line + "\n")
        f.write("PATCH FIRST NONE LAST NONE\n")
        f.write("END\n")


def generate_pdb(file: str, cg: ConnectionGraph, is_three_ind=False, rep_param=None, force_indexes=None):
    atom_indexs = []
    n_atoms = len(cg.nodes)
    if force_indexes is not None:
        atom_indexs = force_indexes
    else:
        if is_three_ind:
            for i in range(len(cg.nodes)):
                if rep_param is not None:
                    atom_indexs.append(generate_un_three_symbol_par(cg.nodes[i].Atom, i, rep_param))
                else:
                    atom_indexs.append(generate_un_three_symbol_par(cg.nodes[i].Atom, i, n_atoms))
        else:
            for i in range(len(cg.nodes)):
                atom_indexs.append(cg.nodes[i].Atom)
                if cg.nodes[i].Atom in atom_alias.keys():
                    atom_indexs.append(atom_alias[cg.nodes[i].Atom])
                else:
                    atom_indexs.append(cg.nodes[i].Atom)
    x, y, z = [], [], []
    with open(file, 'w') as f:
        for i in range(len(cg.nodes)):
            line = "ATOM{0:>7}".format(i + 1) + " {0:>3}   MOL A   1     ".format(atom_indexs[i])
            line += "{0:>7.3f} {1:>7.3f} {2:>7.3f}".format(cg.nodes[i].x[0], cg.nodes[i].x[1], cg.nodes[i].x[2])
            line += "  0.00  0.00      A"
            x.append(cg.nodes[i].x[0])
            y.append(cg.nodes[i].x[1])
            z.append(cg.nodes[i].x[2])
            f.write(line + "\n")
        f.write("TER\n")
        # bonds = cg.get_all_pairs()
        # for bond in bonds:
        #     line = "CONECT {0:>5} {1:>5}".format(bond[0], bond[1]) + "\n"
        #     f.write(line)
        # f.write("END\n")
    center_x = (max(x) - min(x)) / 2
    center_y = (max(y) - min(y)) / 2
    center_z = (max(z) - min(z)) / 2
    print('cellOrigin ', center_x, center_y, center_z)
