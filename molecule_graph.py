import numpy as np
from constants import ATOMS_TO_COVALENT_RADII as ATR
from typing import List, Union
import numpy.linalg as lg


class NodeData:
    Atom: str
    x: np.ndarray
    index: int

    def __init__(self, Atom: str, x: np.ndarray, index: int):
        self.Atom = Atom
        self.x = x
        self.index = index

    def get_length(self, other: 'NodeData'):
        return np.linalg.norm(self.x - other.x)


class ConnectionGraph:
    nodes: dict
    connections: dict

    def __init__(self):
        self.nodes = {}
        self.connections = {}
        self._atom_check = AtomConnectionChecker()

    def read_from_psf(self, file: str):
        n_atoms = 0
        with open(file, 'r') as f:
            for line in f:
                if "!NATOM" in line:
                    break
            i = 0
            for line in f:
                if len(line.split()) < 3:
                    break
                charge = line.split()[5][0]
                self.add_nodes(charge, np.array([0, 0, 0]), i)
                i += 1
            for line in f:
                if "!NBOND" in line:
                    break
            for line in f:
                if len(line.split()) < 1:
                    break
                indxs = [int(v) for v in line.split()]
                for i in range(len(indxs) // 2):
                    idx1 = indxs[2 * i]
                    idx2 = indxs[2 * i + 1]
                    self.add_connection(idx1 - 1, idx2 - 1, 1)

    def read_only_connections(self, file: str):
        with open(file, 'r') as f:
            for line in f:
                if "!NBOND" in line:
                    break
            for line in f:
                if len(line.split()) < 1:
                    break
                indxs = [int(v) for v in line.split()]
                for i in range(len(indxs) // 2):
                    idx1 = indxs[2 * i]
                    idx2 = indxs[2 * i + 1]
                    self.add_connection(idx1 - 1, idx2 - 1, 1)

    def set_bound_matrix(self):
        atoms = [i for i in self.nodes.keys()]
        n_atoms = len(atoms)
        self.__distance_matrix = np.zeros((len(self.nodes), len(self.nodes)))
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                bond = lg.norm(self.nodes[atoms[i]].x - self.nodes[atoms[j]].x)


    def add_nodes_from_geometry(self, charges: list, coords: np.ndarray):
        for i in range(len(charges)):
             nd = NodeData(charges[i], coords[3 * i: 3 * i + 3], i)
             self.nodes.update({i: nd})
        self.__distance_matrix = np.zeros((len(self.nodes), len(self.nodes)))

    def add_nodes(self, atom, coord, i):
        nd = NodeData(atom, coord, i)
        self.nodes.update({i: nd})
        self.__distance_matrix = np.zeros((len(self.nodes), len(self.nodes)))

    def remove_node(self, i):
        del(self.nodes[i])

    def has_old_index(self, indx: int) -> bool:
        for node in self.nodes.values():
            if node.index == indx:
                return True
        return False

    def print_struct(self) -> str:
        result = ""
        for key, node in self.nodes.items():
            result += node.Atom + " " + str(node.x[0]) + " " + str(node.x[1]) + " " + str(node.x[2]) + "\n"
        return result

    def add_connection(self, indx1: int, indx2: int, bond: float):
        if indx1 in self.connections.keys():
            con = self.connections[indx1]
            con.append(indx2)
        else:
            self.connections.update({indx1: [indx2]})
        if indx2 in self.connections.keys():
            con = self.connections[indx2]
            con.append(indx1)
        else:
            self.connections.update({indx2: [indx1]})
        self.__distance_matrix[indx1, indx2] = bond
        self.__distance_matrix[indx2, indx1] = bond

    def split(self) -> List['ConnectionGraph']:
        all_atoms = [i for i in self.nodes.keys()]
        bounded_atoms = self.get_bounded_nodes(0)
        bounded_atoms.sort()
        result = [bounded_atoms]
        result_gh = []
        all_bounded_atoms = bounded_atoms.copy()
        while True:
            is_splited = True
            prev = -1
            for atom in all_atoms:
                if atom not in all_bounded_atoms:
                    is_splited = False
                    prev = atom
                    break
            if not is_splited:
                bounded_atoms = self.get_bounded_nodes(prev)
                bounded_atoms.sort()
                result.append(bounded_atoms)
                all_bounded_atoms.extend(bounded_atoms)
            if is_splited:
                break
        for new_gh_atoms in result:
            atoms = [self.nodes[i].Atom for i in new_gh_atoms]
            coords = []
            indexes = [self.nodes[i].index for i in new_gh_atoms]
            for i in new_gh_atoms:
                coords.extend(self.nodes[i].x.tolist())

            coords = np.array(coords)
            gh = ConnectionGraph()
            gh.add_nodes_from_geometry(atoms, coords)
            gh.set_bonds()
            for j in range(len(indexes)):
                gh.nodes[j].index = indexes[j]
            result_gh.append(gh)
        return result_gh

    def get_bounded_nodes(self, prev_pos: Union[int, List[int]]) -> List[int]:
        if type(prev_pos) is int:
            if prev_pos not in self.connections.keys():
                return [prev_pos]
            prev_pos = [prev_pos]
        next_pos = prev_pos
        for cur_pos in self.connections[next_pos[-1]]:
            if cur_pos not in prev_pos:
                prev_pos.append(cur_pos)
                prev_pos = self.get_bounded_nodes(prev_pos)
        return prev_pos

    def iter_nodes(self):
        keys = self.nodes.keys()
        keys = sorted(keys)
        return (self.nodes[i] for i in keys)

    def set_bonds(self):
        atoms = [i for i in self.nodes.keys()]
        n_atoms = len(atoms)
        for i in range(n_atoms):
            if i % 500 == 0:
                print("End " + str(i))
            for j in range(i + 1, n_atoms):
                distance_threshold = (ATR[self.nodes[atoms[i]].Atom] + ATR[self.nodes[atoms[j]].Atom]) * 1.3
                bond = lg.norm(self.nodes[atoms[i]].x - self.nodes[atoms[j]].x)
                if bond < distance_threshold:
                    self.add_connection(i, j, bond)

    def get_indexes(self, numbers: List[int]):
        pass

    def get_all_pairs(self):
        result = []
        prev_sets = []
        for key in range(len(self.nodes)):
            values = self.connections[key]
            for value in values:
                if [key, value] not in prev_sets \
                        and [value, key] not in prev_sets:
                    prev_sets.append([key, value])
                    result.append([key, value])
        aprove_result = []
        for bond in result:
            aprove_result.append([self.nodes[bond[0]].index+1,
                                  self.nodes[bond[1]].index+1])
        return aprove_result

    def get_all_triplets(self):
        prev_sets = []
        result = []

        for pos1 in range(len(self.nodes)):
            values = self.connections[pos1]
            for pos2 in values:
                for pos3 in self.connections[pos2]:
                    if pos3 == pos1:
                        continue
                    if [pos1, pos2, pos3] not in prev_sets \
                            and [pos3, pos2, pos1] not in prev_sets:
                        prev_sets.append([pos1, pos2, pos3])
                        result.append([pos1, pos2, pos3])
        aprove_result = []
        for angle in result:
            pos1 = angle[0]
            pos2 = angle[1]
            pos3 = angle[2]
            if pos3 in self.connections[pos1]:
                if [pos1, pos3, pos2] in aprove_result or [pos2, pos1, pos3] in aprove_result:
                    continue
            aprove_result.append(angle)
        result = aprove_result
        aprove_result = []
        for angle in result:
            aprove_result.append([self.nodes[angle[0]].index+1,
                                  self.nodes[angle[1]].index+1,
                                  self.nodes[angle[2]].index+1])
        return aprove_result

    def get_all_quartet(self):
        prev_sets = []
        result = []
        for pos1 in range(len(self.nodes)):
            values = self.connections[pos1]
            for pos2 in values:
                for pos3 in self.connections[pos2]:
                    if pos3 == pos1:
                        continue
                    for pos4 in self.connections[pos3]:
                        if pos1 == pos4 or pos2 == pos4:
                            continue
                        if [pos1, pos2, pos3, pos4] not in prev_sets \
                                and [pos4, pos3, pos2, pos1] not in prev_sets:
                            prev_sets.append([pos1, pos2, pos3, pos4])
                            result.append([pos1, pos2, pos3, pos4])
        approve_result = []
        for dihedr in result:
            approve_result.append([self.nodes[dihedr[0]].index + 1,
                                  self.nodes[dihedr[1]].index + 1,
                                  self.nodes[dihedr[2]].index + 1,
                                  self.nodes[dihedr[3]].index + 1])

        return approve_result

    def set_all_hydrogen_bonds(self):
        # TODO: set bond threshold as func from vdw radii
        list_acceptable_type_atoms = ["N", "O", "F", "P", "S", "Cl"]
        atoms = [i for i in self.nodes.keys()]
        n_atoms = len(atoms)
        for i in range(n_atoms):
            if self.nodes[i].Atom == "H":
                heavy_atom = self.nodes[self.connections[i][0]].Atom
                if heavy_atom not in list_acceptable_type_atoms:
                    continue
                for j in range(n_atoms):
                    if i == j:
                        continue
                    if j in self.connections[i]:
                        continue
                    if self.nodes[j].Atom not in list_acceptable_type_atoms:
                        continue
                    if self.nodes[j].Atom in ["P", "S"]:
                        bond = lg.norm(self.nodes[j].x - self.nodes[i].x)
                        if bond < 2.5:
                            self.add_connection(i, j, bond)
                    else:
                        bond = lg.norm(self.nodes[j].x - self.nodes[i].x)
                        if bond < 2.0:
                            self.add_connection(i, j, bond)

    def get_shortest_bond(self, other: 'ConnectionGraph'):
        short_distance = lg.norm(self.nodes[0].x - other.nodes[0].x)
        short_pair = [0, 0]
        for i, ivalue in self.nodes.items():
            for j, jvalue in other.nodes.items():
                n_dist = lg.norm(ivalue.x - jvalue.x)
                if short_distance > n_dist:
                    short_pair = [i, j]
                    short_distance = n_dist
        return self.nodes[short_pair[0]].index, other.nodes[short_pair[1]].index, short_distance

    def find_close_to_shortest(self, other: 'ConnectionGraph', bond: list, shortest: float):
        pairs = []
        for i, ivalue in self.nodes.items():
            for j, jvalue in other.nodes.items():
                n_dist = lg.norm(ivalue.x - jvalue.x)
                if abs(n_dist - shortest) < 0.2:
                    if i != bond[0] and j != bond[1]:
                        pairs.append([self.nodes[i].index, other.nodes[j].index, n_dist])
        return pairs

    def set_intermolecular_connections(self):
        list_of_parts = self.split()
        for i in range(len(list_of_parts)):
            for j in range(i+1, len(list_of_parts)):
                idx1, idx2, bond = list_of_parts[i].get_shortest_bond(list_of_parts[j])
                self.add_connection(idx1, idx2, bond)
                pairs = list_of_parts[i].find_close_to_shortest(list_of_parts[j], [idx1, idx2], bond)
                for pair in pairs:
                    self.add_connection(pair[0], pair[1], pair[2])

    def remove_bond(self, indx1: int, indx2: int):
        con = self.connections[indx1]
        con.remove(indx2)
        con = self.connections[indx2]
        con.remove(indx1)
        self.__distance_matrix[indx1, indx2] = 0
        self.__distance_matrix[indx2, indx1] = 0

    def del_node(self, indx: int):
        connections = self.connections[indx]
        for other_indx in connections:
            self.remove_bond(indx, other_indx)
        self.nodes[indx] = None


    def sweap_unlogic_bonds(self):
        for indx, connection in self.connections.items():
            if self._atom_check.need_to_correct(self.nodes[indx], connection):
                to_exclude = self._atom_check.check_connections(self.nodes[indx],
                                                                [self.nodes[i] for i in connection])
                for node_id in to_exclude:
                    self.remove_bond(indx, node_id)


class AtomConnectionChecker:

    def __init__(self):
        self.working_atoms = ["H", "C", "O", "N"]
        self.max_bond_dict = {"H": 1, "C": 4, "O": 3, "N": 4}
        self.average_bond = {"H": 1.08, "C": 1.4, "N": 1.3, "O": 1.2}

    def need_to_correct(self, central_node: NodeData, connected_nodes: list) -> bool:
        if not central_node.Atom in self.working_atoms:
            return False
        if len(connected_nodes) <= self.max_bond_dict[central_node.Atom]:
            return False
        return True

    def check_connections(self, central_node: NodeData, connected_nodes: list) -> List[int]:
        if not central_node.Atom in self.working_atoms:
            return []
        if len(connected_nodes) < self.max_bond_dict[central_node.Atom]:
            return []
        need_to_exclude = len(connected_nodes) - self.max_bond_dict[central_node.Atom]
        bond_lenghts = []
        to_exclude = []
        for node in connected_nodes:
            average = self.average_bond[central_node.Atom]
            if node.Atom == "H":
                average = self.average_bond["H"]
            else:
                if average < self.average_bond[node.Atom]:
                    average = self.average_bond[node.Atom]
            bond_lenghts.append(abs(central_node.get_length(node) - average)/average)
        for _ in range(need_to_exclude):
            max_ind = 0
            max_d = bond_lenghts[0]
            for i in range(len(bond_lenghts)):
                if bond_lenghts[i] > max_d:
                    max_ind = i
                    max_d = bond_lenghts[i]
            bond_lenghts[max_ind] = 0.0
            to_exclude.append(connected_nodes[max_ind].index)
        return to_exclude