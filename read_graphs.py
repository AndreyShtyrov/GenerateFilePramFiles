from molecule_graph import ConnectionGraph
from constants import CHARGE_TO_CODE
from files_generatos import generate_top, generate_psf, generate_par, generate_pdb
import numpy as np
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--to_zero', type=bool, default=False)
parser.add_argument('--check_connections', type=bool, default=False)
parser.add_argument('--only_three', type=bool, default=False)
parser.add_argument('--only_pdb', type=bool, default=False)
parser.add_argument('--rep_param', type=int, default=-1)
parser.add_argument('--force_index', type=str, default="")


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


def read_geom_sl(file: str):
    charges = []
    coords = []

    with open(file, 'r') as f:
        for line in f:
            if len(line.split()) < 3:
                continue
            charge = line.split()[0].split("-")[0]
            charges.append(charge)
            coords.extend([float(v) for v in line.replace('\n', '').split()[2:-1]])
    return charges, np.array(coords)


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


def read_f_ind(file: str):
    charges = []

    with open(file, 'r') as f:
        for line in f:
            if len(line.split()) < 3:
                continue
            charge = line.split()[0].split("-")[1].replace("Z", "")
            charges.append(charge)
    return charges


def main():
    args = parser.parse_args()
    to_zero = args.to_zero
    to_three = args.only_three
    only_pdb = args.only_pdb
    rep_param = args.rep_param
    if rep_param < 0:
        rep_param = None

    improp = None
    charges, coords = read_geom("struct.xyz")
    charge_file = Path.cwd() / "charges.txt"
    if (Path.cwd() / "imps.txt").is_file():
        improp = read_impropers(str(Path.cwd() / "imps.txt"))
    if charge_file.is_file():
        mmcharges = read_charges(str(charge_file))
    else:
        mmcharges = None
    cg = ConnectionGraph()
    cg.add_nodes_from_geometry(charges, coords)
    if not only_pdb:
        cg.set_bonds()

    if mmcharges is not None:
        if len(mmcharges) != len(cg.nodes):
            print("Amount charges should be same to amount atoms")

    if rep_param == None:
        rep_param = len(cg.nodes)
    if args.check_connections:
        cg.sweap_unlogic_bonds()
    force_index = None
    if args.force_index != "":
        force_index = []
        with open(args.force_index, "r") as f:
            for line in f:
                force_index.append(line.split()[0].replace("\n", ""))
    if not only_pdb:
        if force_index != None:
            generate_psf("mol.psf", cg, mmcharges, is_three_ind=to_three, imprs=improp, rep_param=rep_param, force_indexes=force_index)
            generate_par("ff.par", cg, to_zero, is_three_ind=to_three, imprs=improp, rep_param=rep_param, force_indexs=force_index)
            generate_top("ff.top", cg, mmcharges, is_three_ind=to_three, rep_param=rep_param)
        else:
            generate_psf("mol.psf", cg, mmcharges, is_three_ind=to_three, imprs=improp, rep_param=rep_param)
            generate_par("ff.par", cg, to_zero, is_three_ind=to_three, imprs=improp, rep_param=rep_param)
            generate_top("ff.top", cg, mmcharges, is_three_ind=to_three, rep_param=rep_param)
    generate_pdb("mol.pdb", cg, is_three_ind=to_three, rep_param=rep_param, force_indexes=force_index)


if __name__ == '__main__':
    main()
