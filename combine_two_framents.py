from molecule_graph import ConnectionGraph
from constants import CHARGE_TO_CODE
from files_generatos import generate_top, generate_psf, \
    generate_par, generate_pdb, \
    generate_un_three_symbol_par, generate_par_from_MMContainer
from rotateXYZ import RX, RY, RZ
from combine_params import MMParmContainer, MMParm
import numpy as np
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--to_zero', type=bool, default=False)
parser.add_argument('--only_three', type=bool, default=False)
parser.add_argument('--rep_param', type=int, default=-1)


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

def read_impropers(file: str):
    result = []
    with open(file, 'r') as f:
        for line in f:
            col = line.split()
            result.append([int(col[0]), int(col[1]), int(col[2]), int(col[3])])
    return result


def read_charges(file: str):
    result = []
    with open(file, 'r') as f:
        for line in f:
            result.append(float(line))
    return result


def read_geom(file: str):
    charges = []
    coords = []
    with open(file, 'r') as f:
        for line in f:
            if len(line.split()) < 3:
                continue
            try:
                charge = int(line.split()[0])
                charge = CHARGE_TO_CODE[charge]
            except:
                charge = line.split()[0]
            charges.append(charge)
            coords.extend([float(v) for v in line.replace('\n', '').split()[1:]])
    return charges, np.array(coords)


if __name__ == '__main__':
    with open('config.cfg', 'r') as f:
        geoms = []
        fragment_names = [None, None]
        par_files = []
        connections = []
        charges = []
        psf_par = [None, None]
        for line in f:
            if "XYZ" in line:
                _charge = []
                coord = []
                for line in f:
                    if len(line.split()) < 3:
                        break
                    try:
                        charge = int(line.split()[0])
                        charge = CHARGE_TO_CODE[charge]
                    except:
                        charge = line.split()[0]
                    _charge.append(charge)
                    coord.extend([float(v) for v in line.replace('\n', '').split()[1:]])
                geoms.append(np.array(coord))
                charges.append(_charge)
            if "CONNECTION" in line:
                connections.append([int(line.split()[-2]), int(line.split()[-1])])
            if "PARAMETORS" in line:
                par_files.append(line.split()[-1].replace("\n", ""))
            if "FRAG_NAME1" in line:
                fragment_names[0] = line.split()[-1].replace("\n", "")
            if "FRAG_NAME2" in line:
                fragment_names[1] = line.split()[-1].replace("\n", "")
            if "PSF_FILE1" in line:
                psf_par[0] = line.split()[-1].replace("\n", "")
            if "PSF_FILE2" in line:
                psf_par[1] = line.split()[-1].replace("\n", "")
    if len(charges) != 2 and len(geoms) != 2 and len(fragment_names) != 2 and len(connections) != 2:
        print("You need provide full infromations for 2 fragments")

    fg1_abr = fragment_names[0]
    fg2_abr = fragment_names[1]
    charges1 = charges[0]
    charges2 = charges[1]
    coords1 = geoms[0]
    coords2 = geoms[1]
    a1, a2 = connections[0][0] - 1, connections[0][1] - 1
    a3, a4 = connections[1][0] - 1, connections[1][1] - 1
    improp = None
    charge_file = Path.cwd() / "charges1.txt"
    if charge_file.is_file():
        mmcharges1 = read_charges(str(charge_file))
    else:
        mmcharges1 = None
    if psf_par[0] is not None:
        mmcharges1 = read_charges_from_psf(psf_par[0])
    charge_file = Path.cwd() / "charges2.txt"
    if charge_file.is_file():
        mmcharges2 = read_charges(str(charge_file))
    else:
        mmcharges2 = None
    if psf_par[1] is not None:
        mmcharges2 = read_charges_from_psf(psf_par[1])
    if (Path.cwd()/"imps.txt").is_file():
        improp = read_impropers(str(Path.cwd()/"imps.txt"))

    cg1 = ConnectionGraph()
    cg1.add_nodes_from_geometry(charges1, coords1)
    cg1.set_bonds()
    cg2 = ConnectionGraph()
    cg2.add_nodes_from_geometry(charges2, coords2)
    cg2.set_bonds()

    cg1.remove_bond(a1, a2)
    cg1_i = []
    for i in range(len(cg1.nodes.values())):
        cg1_i.append(cg1.nodes[i].Atom + str(i))
    cg2_i = []
    for i in range(len(cg2.nodes.values())):
        cg2_i.append(cg2.nodes[i].Atom + str(i))
    cg2.remove_bond(a3, a4)

    cgs1 = cg1.split()
    cgs2 = cg2.split()
    if cgs1[0].has_old_index(a1):
        frg1 = cgs1[0]
    else:
        frg1 = cgs1[1]
    if cgs2[0].has_old_index(a4):
        frg2 = cgs2[0]
    else:
        frg2 = cgs2[1]

    IndM = np.zeros((3, 3))
    IndM[0, 0] = 1
    IndM[1, 1] = 1
    IndM[2, 2] = 1
    vector1 = cg1.nodes[a2].x - cg1.nodes[a1].x
    vector2 = cg2.nodes[a4].x - cg2.nodes[a3].x
    ox1 = cg1.nodes[a1].x
    ox2 = cg2.nodes[a3].x

    alpha1 = np.arccos((np.array([vector1[0], vector1[1]])/np.linalg.norm(np.array([vector1[0], vector1[1]])))[0])
    v1 = RZ(-alpha1).dot(vector1)
    if abs(v1[1]) > 0.0001:
        v1 = RZ(alpha1).dot(vector1)
        alpha1 = -alpha1
    beta1 = np.arccos((v1 / np.linalg.norm(v1))[2])
    v1 = RY(beta1).dot(v1)
    alpha2 = np.arccos((np.array([vector2[0], vector2[1]])/np.linalg.norm(np.array([vector2[0], vector2[1]])))[0])
    v2 = RZ(-alpha2).dot(vector2)
    if abs(v2[1]) > 0.0001:
        v2 = RZ(alpha2).dot(vector2)
        alpha2 = -alpha2
    beta2 = np.arccos((v2 / np.linalg.norm(v2))[2])
    v2 = RY(beta2).dot(v2)
    if True:
        IndM = np.zeros((3, 3))
        IndM[0, 0] = 1
        IndM[1, 1] = 1
        IndM[2, 2] = -1

    at_n1 = [frg1.nodes[i].Atom for i in range(len(frg1.nodes.values()))]
    at_n2 = [frg2.nodes[i].Atom for i in range(len(frg2.nodes.values()))]
    xyz1 = [frg1.nodes[i].x for i in range(len(frg1.nodes.values()))]
    xyz2 = [frg2.nodes[i].x for i in range(len(frg2.nodes.values()))]
    mid1_system = []
    mid2_system = []
    mid3_system = []
    for i in range(len(xyz2)):
        xyz = xyz2[i] - ox2
        new_line = at_n2[i] + " " + str(xyz[0]) + " " + str(xyz[1])
        new_line += " " + str(xyz[2]) + "\n"
        mid3_system.append(new_line)
        mxyz = RY(beta2).dot(RZ(-alpha2).dot(xyz))
        new_line = at_n2[i] + " " + str(mxyz[0]) + " " + str(mxyz[1])
        new_line += " " + str(mxyz[2]) + "\n"
        mid1_system.append(new_line)
        myyz = IndM.dot(mxyz)
        mxyz = RZ(alpha1).dot(RY(-beta1).dot(mxyz))

        new_line = at_n2[i] + " " + str(mxyz[0]) + " " + str(mxyz[1])
        new_line += " " + str(mxyz[2]) + "\n"
        mid2_system.append(new_line)
        xyz = mxyz + ox1
        xyz2[i] = xyz
    xyz1.extend(xyz2)
    at_n1.extend(at_n2)
    with open("result.xyz", 'w') as f:
        for i in range(len(at_n1)):
            new_line = at_n1[i] + " " + str(xyz1[i][0]) + " " + str(xyz1[i][1])
            new_line += " " + str(xyz1[i][2]) + "\n"
            f.write(new_line)

    xyz1 = np.array(xyz1).reshape(-1)
    new_full_graph = ConnectionGraph()
    new_full_graph.add_nodes_from_geometry(at_n1, xyz1)
    shift = 0
    end_index_to_frg1 = []
    rev_index_to_frg1 = {}
    end_index_to_frg2 = []
    rev_index_to_frg2 = {}
    for i in range(len(frg1.nodes.values())):
        for con in frg1.connections[i]:
            new_full_graph.add_connection(i, con, 1)
        shift += 1
        end_index_to_frg1.append(frg1.nodes[i].index)
        rev_index_to_frg1[frg1.nodes[i].index] = i
        end_index_to_frg2.append(None)
    for j in range(len(frg2.nodes.values())):
        if frg2.nodes[j].index == a4:
            new_full_graph.add_connection(rev_index_to_frg1[a1], j + shift, 1)
    for i in range(len(frg2.nodes.values())):
        for con in frg2.connections[i]:
            new_full_graph.add_connection(i + shift, con + shift, 1)
        end_index_to_frg1.append(None)
        end_index_to_frg2.append(frg2.nodes[i].index)
        rev_index_to_frg2[frg2.nodes[i].index] = i + shift

    n_atoms = len(frg2.nodes.values())
    new_charges = None
    if mmcharges1 != None and mmcharges2 != None:
        new_charges = []
        for i in range(len(frg1.nodes.values())):
            new_charges.append(mmcharges1[frg1.nodes[i].index])
        for i in range(len(frg2.nodes.values())):
            new_charges.append(mmcharges2[frg2.nodes[i].index])

    
    frg1_par_container = MMParmContainer()
    with open(par_files[0], 'r') as f:
        for line in f:
            if "NONBONDED  NBXMOD 5" in line:
                break
            if "!" in line:
                continue
            if len(line.split()) < 3:
                continue
            frg1_par_container.append(MMParm(line))
        for i, line in enumerate(f):
            if i == 0:
                continue
            if len(line.split()) < 3:
                continue
            frg1_par_container.append(MMParm.read_as_vdw(line))
    frg2_par_container = MMParmContainer()
    with open(par_files[1], 'r') as f:
        for line in f:
            if "NONBONDED  NBXMOD 5" in line:
                break
            if "!" in line:
                continue
            if len(line.split()) < 3:
                continue
            frg2_par_container.append(MMParm(line))
        for i, line in enumerate(f):
            if i == 0:
                continue
            if len(line.split()) < 3:
                continue
            frg2_par_container.append(MMParm.read_as_vdw(line))

    cg1_a_i = []
    fg1_a_i = {}
    fg1_a_i_temp = {}
    fg1_a_i_temp_rev = {}
    cg2_a_i = []
    fg2_a_i = {}
    fg2_a_i_temp = {}
    fg2_a_i_temp_rev = {}
    end_indexs = []
    temp_indexs = []
    if psf_par[0] is None:
        for i in range(len(cg1.nodes.values())):
            cg1_a_i.append(generate_un_three_symbol_par(cg1.nodes[i].Atom, i, len(cg1.nodes.values())))
    else:
        cg1_a_i = read_indexs_from_psf(psf_par[0])
    if psf_par[1] is None:
        for i in range(len(cg2.nodes.values())):
            cg2_a_i.append(generate_un_three_symbol_par(cg2.nodes[i].Atom, i, len(cg2.nodes.values())))
    else:
        cg2_a_i = read_indexs_from_psf(psf_par[1])
    for i in range(len(frg1.nodes.values())):
        par_line = cg1_a_i[frg1.nodes[i].index]
        fg1_a_i[par_line] = par_line
        if fg1_abr is not None:
            _str = fg1_a_i[par_line][1:]
            fg1_a_i[par_line] = fg1_abr + _str
        fg1_a_i_temp[par_line] = "a" + par_line
        fg1_a_i_temp_rev["a" + par_line] = fg1_a_i[par_line]
        end_indexs.append(fg1_a_i[par_line])
        temp_indexs.append(fg1_a_i_temp[par_line])
    for i in range(len(frg2.nodes.values())):
        par_line = cg2_a_i[frg2.nodes[i].index]
        fg2_a_i[par_line] = par_line
        if fg2_abr is not None:
            _str = fg2_a_i[par_line][1:]
            fg2_a_i[par_line] = fg2_abr + _str
        fg2_a_i_temp[par_line] = "b" + par_line
        fg2_a_i_temp_rev["b" + par_line] = fg2_a_i[par_line]
        end_indexs.append(fg2_a_i[par_line])
        temp_indexs.append(fg2_a_i_temp[par_line])

    generate_par("temp_full.par", new_full_graph, to_zero=True, force_indexs=temp_indexs)
    temp_par_container = MMParmContainer()

    with open("temp_full.par", 'r') as f:
        for line in f:
            if "NONBONDED  NBXMOD 5" in line:
                break
            if "!" in line:
                continue
            if len(line.split()) < 3:
                continue
            temp_par_container.append(MMParm(line))
        for i, line in enumerate(f):
            if i == 0:
                continue
            if len(line.split()) < 3:
                continue
            temp_par_container.append(MMParm.read_as_vdw(line))

    if not (Path.cwd() / "stitching.par").is_file():
        stitching_container = MMParmContainer()
        for bond in temp_par_container.bonds:
            if (bond.atoms[0][0] == "a" or bond.atoms[1][0] == "a") and (
                    bond.atoms[0][0] == "b" or bond.atoms[1][0] == "b"):
                stitching_container.append(bond)
        for angle in temp_par_container.angles:
            if (angle.atoms[0][0] == "a" or angle.atoms[1][0] == "a" or angle.atoms[2][0] == "a") and (
                    angle.atoms[1][0] == "b" or angle.atoms[1][0] == "b" or angle.atoms[2][0] == "b"):
                stitching_container.append(angle)
        for dh in temp_par_container.dih:
            if (dh.atoms[0][0] == "a" or dh.atoms[1][0] == "a"
            or dh.atoms[2][0] == "a" or dh.atoms[3][0] == "a") and (
                    dh.atoms[0][0] == "b" or dh.atoms[1][0] == "b"
            or dh.atoms[2][0] == "b" or dh.atoms[3][0] == "b"):
                stitching_container.append(dh)
        generate_par_from_MMContainer("stitching.par", stitching_container)
        print("You need to correct stitching par and rerun program")
    else:
        stitching_container = MMParmContainer()
        with open("stitching.par", 'r') as f:
            for line in f:
                if "NONBONDED  NBXMOD 5" in line:
                    break
                if "!" in line:
                    continue
                if len(line.split()) < 3:
                    continue
                stitching_container.append(MMParm(line))

    end_par_container = MMParmContainer.generate_from_psf(new_full_graph, to_zero=True, force_indexs=end_indexs)
    frg1_par_container.rename_all_atoms(fg1_a_i)
    frg2_par_container.rename_all_atoms(fg2_a_i)
    stitching_container.rename_all_atoms(fg1_a_i_temp_rev)
    stitching_container.rename_all_atoms(fg2_a_i_temp_rev)
    result_par_container = MMParmContainer()
    for bond in end_par_container.bonds:
        par_line = frg1_par_container.get_param(bond.print_line())
        par_line = frg2_par_container.get_param(par_line.print_line())
        par_line = stitching_container.get_param(par_line.print_line())
        result_par_container.append(par_line)
    for angle in end_par_container.angles:
        par_line = frg1_par_container.get_param(angle.print_line())
        par_line = frg2_par_container.get_param(par_line.print_line())
        par_line = stitching_container.get_param(par_line.print_line())
        result_par_container.append(par_line)
    for dh in end_par_container.dih:
        par_line = frg1_par_container.get_param(dh.print_line())
        par_line = frg2_par_container.get_param(par_line.print_line())
        par_line = stitching_container.get_param(par_line.print_line())
        result_par_container.append(par_line)
    for vdw in end_par_container.vdws:
        par_line = frg1_par_container.get_param(vdw.print_line())
        par_line = frg2_par_container.get_param(par_line.print_line())
        par_line = stitching_container.get_param(par_line.print_line())
        result_par_container.append(par_line)

    generate_psf("result.psf", new_full_graph, new_charges, is_three_ind=True, force_indexes=end_indexs)
    generate_pdb("result.pdb", new_full_graph, new_charges, force_indexes=end_indexs)
    generate_par_from_MMContainer("result.par", result_par_container, imprs=improp, atom_indexes=end_indexs)

