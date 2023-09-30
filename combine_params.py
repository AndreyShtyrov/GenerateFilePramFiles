from enum import Enum
from typing import List

from constants import VDW_A_EC
from molecule_graph import ConnectionGraph
from constants import atom_alias, alphabet
import sys

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

class MMType(Enum):
    Bond = 1
    Angle = 2
    Dihedral = 3
    VDW = 4
    D1 = 5
    D2 = 6
    D3 = 7
    D4 = 8


class MMParm:

    def __init__(self, line: str):
        self._atom_list = []
        if len(line.split()) == 4:
            self._type = MMType.Bond
            self._atom_list = line.split()[0:2]
            self._num_params = [float(v) for v in line.split()[2:]]
        elif len(line.split()) == 5:
            self._type = MMType.Angle
            self._atom_list = line.split()[:3]
            self._num_params = [float(v) for v in line.split()[3:]]
        elif len(line.split()) == 7:
            if "." in line.split()[1]:
                self._type = MMType.VDW
                self._atom_list = [line.split()[0]]
                self._num_params = [float(v) for v in line.split()[1:]]
            else:
                d_num = int(line.split()[-2])
                if d_num == 1:
                    self._type = MMType.D1
                elif d_num == 2:
                    self._type = MMType.D2
                elif d_num == 3:
                    self._type = MMType.D3
                elif d_num == 4:
                    self._type = MMType.D4
                self._atom_list = line.split()[0:4]
                v1 = float(line.split()[4])
                v2 = int(line.split()[5])
                v3 = float(line.split()[6])
                self._num_params = [v1, v2, v3]
        elif 'VDW' in line.split()[0]:
            self._type = MMType.VDW
            self._atom_list = [line.split()[1]]
            self._num_params = [float(v) for v in line.split()[2:]]

    @staticmethod
    def read_as_vdw(line: str) -> 'MMParm':
        mp = MMParm("VDW " + line)
        return mp

    @property
    def type(self) -> MMType:
        return self._type

    @property
    def atoms(self) -> List[str]:
        return self._atom_list

    @property
    def params(self) -> List[float]:
        return self._num_params

    def print_line(self) -> str:
        line = ''
        for atom in self.atoms:
            line += ' {0:>4}'.format(atom) 
        for param in self.params:
            if type(param) == int:
                line += ' {0:>3}'.format(param)
                continue
            line += ' {0:>9.4f}'.format(param) 
        return line + '\n'

    def print_gauss_line(self, with_out_Z:bool =False) -> str:
        if self.type is MMType.Bond:
            line = "HrmStr1 "
            for atom in self.atoms:
                if with_out_Z:
                    line += ' {0:>4}'.format(atom)
                else:
                    line += ' {0:>4}'.format(atom) + "Z"
            for param in self.params:
                if type(param) == int:
                    line += ' {0:>3}'.format(param)
                    continue
                line += ' {0:>9.4f}'.format(param)
        elif self.type is MMType.Angle:
            line = "HrmBnd1 "
            for atom in self.atoms:
                if with_out_Z:
                    line += ' {0:>4}'.format(atom)
                else:
                    line += ' {0:>4}'.format(atom) + "Z"
            for param in self.params:
                if type(param) == int:
                    line += ' {0:>3}'.format(param)
                    continue
                line += ' {0:>9.4f}'.format(param)
        elif self.type is MMType.VDW:
            line = "VDW "
            for atom in self.atoms:
                if with_out_Z:
                    line += ' {0:>4}'.format(atom)
                else:
                    line += ' {0:>4}'.format(atom) + "Z"
            line += " " + str(self.params[2]) + " "
            line += str(-2 * self.params[1])
        else:
            line = "AmbTrs "
            for atom in self.atoms:
                if with_out_Z:
                    line += ' {0:>4}'.format(atom)
                else:
                    line += ' {0:>4}'.format(atom) + "Z"
            for param in self.params:
                if type(param) == int:
                    line += ' {0:>3}'.format(param)
                    continue
                line += ' {0:>9.4f}'.format(param)
        return line + "\n"

    def get_atom_line(self) -> str:
        line = ""
        for atom in self.atoms:
            line += ' {0:>4}'.format(atom)
        return line

    def __eq__(self, other: 'MMParm'):
        if self.type != other.type:
            return False
        if self.type == MMType.Bond:
            if (self.atoms[0] == other.atoms[-1]
                and self.atoms[-1] == other.atoms[0]) or (self.atoms[0] == other.atoms[0]
                                                         and self.atoms[-1] == other.atoms[-1]):
                return True
        elif self.type == MMType.Angle:
            if self.atoms[1] != other.atoms[1]:
                return False
            if (self.atoms[0] == other.atoms[-1]
                and self.atoms[-1] == other.atoms[0]) or (self.atoms[0] == other.atoms[0]
                                                         and self.atoms[-1] == other.atoms[-1]):
                return True
        elif self._type == MMType.VDW:
            if self.atoms[0] != other.atoms[0]:
                return False
            return True
        elif self.type == MMType.D1 or self.type \
                == MMType.D2 or self.type \
                == MMType.D3 or self.type == MMType.D4:
            t1 = False
            t2 = False
            if (self.atoms[0] == other.atoms[-1]
                and self.atoms[-1] == other.atoms[0]) and (self.atoms[1] == other.atoms[-2]
                                                         and self.atoms[-2] == other.atoms[1]):
                t1 = True
            if (self.atoms[0] == other.atoms[0]
                and self.atoms[-1] == other.atoms[-1]) and (self.atoms[1] == other.atoms[1]
                                                         and self.atoms[-2] == other.atoms[-2]):
                t2 = True
            if t1 or t2:
                return True
        return False


class MMParmContainer:

    def __init__(self):

        self.bonds = []
        self.angles = []
        self.dih = []
        self.vdws = []

    def to_list(self) -> List[MMParm]:
        result = []
        for bond in self.bonds:
            result.append(bond)
        for angle in self.angles:
            result.append(angle)
        for dih in self.dih:
            result.append(dih)
        for vdw in self.vdws:
            result.append(vdw)
        return result

    def rename_all_atoms(self, new_atoms_indexs):
        for bond in self.bonds:
            new_atoms_list = []
            for old_atom in bond.atoms:
                if old_atom in new_atoms_indexs.keys():
                    new_atoms_list.append(new_atoms_indexs[old_atom])
                else:
                    new_atoms_list.append(old_atom)
            bond._atom_list = new_atoms_list
        for angle in self.angles:
            new_atoms_list = []
            for old_atom in angle.atoms:
                if old_atom in new_atoms_indexs.keys():
                    new_atoms_list.append(new_atoms_indexs[old_atom])
                else:
                    new_atoms_list.append(old_atom)
            angle._atom_list = new_atoms_list
        for dih in self.dih:
            new_atoms_list = []
            for old_atom in dih.atoms:
                if old_atom in new_atoms_indexs.keys():
                    new_atoms_list.append(new_atoms_indexs[old_atom])
                else:
                    new_atoms_list.append(old_atom)
            dih._atom_list = new_atoms_list
        for vdw in self.vdws:
            new_atoms_list = []
            for old_atom in vdw.atoms:
                if old_atom in new_atoms_indexs.keys():
                    new_atoms_list.append(new_atoms_indexs[old_atom])
                else:
                    new_atoms_list.append(old_atom)
            vdw._atom_list = new_atoms_list

    def append(self, mm_p: MMParm, check: bool=False):
        if mm_p.type == MMType.Bond:
            if not self.check_if_contain(self.bonds, mm_p):
                self.bonds.append(mm_p)
            else:
                if check:
                    print(mm_p.print_line())
        elif mm_p.type == MMType.Angle:
            if not self.check_if_contain(self.angles, mm_p):
                self.angles.append(mm_p)
            else:
                if check:
                    print(mm_p.print_line())
        elif mm_p.type == MMType.VDW:
            if not self.check_if_contain(self.vdws, mm_p):
                self.vdws.append(mm_p)
            else:
                if check:
                    print(mm_p.print_line())
        elif mm_p.type.value > 4:
            if not self.check_if_contain(self.dih, mm_p):
                self.dih.append(mm_p)
            else:
                if check:
                    print(mm_p.print_line())

    @staticmethod
    def check_if_contain(mms: List[MMParm], mm: MMParm):
        for _mm in mms:
            if _mm == mm:
                return True
        return False

    def get_represent(self) -> List[str]:
        result = []
        for bond in self.bonds:
            result.append(bond.print_line())
        for angle in self.angles:
            result.append(angle.print_line())
        for dih in self.dih:
            result.append(dih.print_line())
        for vdw in self.vdws:
            result.append(vdw.print_line())
        return result

    def get_param(self, line) -> MMParm:
        compare_param = MMParm(line)
        if compare_param.type == MMType.Bond:
            for bond in self.bonds:
                if bond == compare_param:
                    return bond
        elif compare_param.type == MMType.Angle:
            for angle in self.angles:
                if angle == compare_param:
                    return angle
        elif compare_param.type == MMType.VDW:
            for vdw in self.vdws:
                if vdw == compare_param:
                    return vdw
        elif compare_param.type.value > 4:
            for angle in self.dih:
                if angle == compare_param:
                    return angle
        return compare_param

    def get_combined_dh(self, with_out_Z:bool=False):
        result = []
        treated_dihs = []
        for dih in self.dih:
            par_line = "AmbTrs"
            for atom in dih.atoms:
                if with_out_Z:
                    par_line += ' {0:>4}'.format(atom)
                else:
                    par_line += ' {0:>4}'.format(atom) + "Z"
            atom_name = dih.get_atom_line()
            if atom_name in treated_dihs:
                continue
            par_line += "    0 180 0 180  "
            d = self.get_param(atom_name + " 0.0 1 0.0 ")
            par_line += " {0:>5.4f}".format(d.params[0])
            d = self.get_param(atom_name + " 0.0 2 0.0 ")
            par_line += " {0:>5.4f}".format(d.params[0])
            d = self.get_param(atom_name + " 0.0 3 0.0 ")
            par_line += " {0:>5.4f}".format(d.params[0])
            d = self.get_param(atom_name + " 0.0 4 0.0 ")
            par_line += " {0:>5.4f}".format(d.params[0])
            par_line += " 1.0 \n"
            result.append(par_line)
            treated_dihs.append(atom_name)
        return result

    @staticmethod
    def generate_from_psf(cg: ConnectionGraph, to_zero: bool, is_three_ind=False, imprs=None,
                 rep_param=None, force_indexs=None) -> 'MMParmContainer':
        atom_indexs = []
        n_atoms = len(cg.nodes)
        result = MMParmContainer()
        if is_three_ind:
            for i in range(len(cg.nodes)):
                if rep_param is None:
                    atom_indexs.append(generate_un_three_symbol_par(cg.nodes[i].Atom, i, n_atoms))
                else:
                    atom_indexs.append(generate_un_three_symbol_par(cg.nodes[i].Atom, i, rep_param))
        else:
            for i in range(len(cg.nodes)):
                atom_indexs.append(cg.nodes[i].Atom + str(i))
        if force_indexs is not None:
            atom_indexs = force_indexs
        bonds = cg.get_all_pairs()
        for b in bonds:
            line = "{0:>5}".format(atom_indexs[b[0] - 1]) + " "
            line += "{0:>5}".format(atom_indexs[b[1] - 1]) + " "
            if to_zero:
                line += "{0:>7.3f}".format(0.0) + " " + "{0:>6.3f}".format(0.0)
            else:
                line += "{0:>7.3f}".format(110.0) + " " + "{0:>6.3f}".format(1.4)
            result.append(MMParm(line))
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
            result.append(MMParm(line))
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
                result.append(MMParm(line))
        for i in range(len(cg.nodes)):
            line = atom_indexs[i]
            line += " 0.00 "
            line += "{0:>9.6f}".format(0.0) + " " + "{0:>9.6f}".format(0.0)
            line += " 0.00 "
            line += "{0:>9.6f}".format(0.0) + " " + "{0:>9.6f}".format(0.0)
            result.append(MMParm.read_as_vdw(line))
        return result


def main():
    only_include = False
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg.upper() == "ONLY_INCLUDE":
            only_include = True
    rep_list = []
    container = MMParmContainer()
    result = []
    with open("rep.txt", 'r') as f:
        for line in f:
            line = line.replace("\n", "")
            rep_list.append((line.split()[0], line.split()[1]))
    with open("prev_par.par", "r") as f:
        for i, line in enumerate(f):
            if "NONBONDED  NBXMOD 5" in line:
                break
            if "!" in line:
                continue
            len_line = len(line.split())
            prev_line = line
            is_replaced = False
            if len_line > 3:
                for rep_type in rep_list:
                    if rep_type[0] + " " in prev_line:
                        is_replaced = True
                        mod_rep_type = rep_type[1][0] + "_x_" + rep_type[1][1:] + " "
                        line = line.replace(rep_type[0] + " ", mod_rep_type)
                        continue
                    if rep_type[0] + "\t" in prev_line:
                        is_replaced = True
                        mod_rep_type = rep_type[1][0] + "_x_" + rep_type[1][1:] + "\t"
                        line = line.replace(rep_type[0] + "\t", mod_rep_type)
                line = line.replace("_x_", "")
                mm_par = MMParm(line)
                if not only_include:
                    container.append(mm_par, is_replaced)
                else:
                    if is_replaced:
                        container.append(mm_par)

    with open("new_par.par", "r") as f:
        for i, line in enumerate(f):
            if "NONBONDED  NBXMOD 5" in line:
                result.append(line)
                break
            if "!" in line:
                result.append(line)
                continue
            len_line = len(line.split())
            if len_line > 3:
                mm_par = container.get_param(line)
                result.append(mm_par.print_line())
                continue
            result.append(line)
        for line in f:
            result.append(line)
    with open("result.par", 'w') as f:
        f.writelines(result)


if __name__ == "__main__":
    main()


