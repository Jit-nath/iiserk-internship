import math
from dataclasses import dataclass
from typing import List, Dict

@dataclass(order=True)
class Atom:
    atomID: int
    atomName: str
    residueName: str
    resID: int
    x: float
    y: float
    z: float
    occupancy: float
    tempFactor: float
    chainID: str

def dist(a1: Atom, a2: Atom) -> float:
    dx = a1.x - a2.x
    dy = a1.y - a2.y
    dz = a1.z - a2.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)

def direction(root: Atom, child: Atom) -> List[float]:
    return [child.x - root.x, child.y - root.y, child.z - root.z]

def normalize(vec: List[float]) -> List[float]:
    magnitude = math.sqrt(sum(v * v for v in vec))
    if magnitude == 0:
        return vec
    return [v / magnitude for v in vec]

def dot_product(a: List[float], b: List[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def cross_product(a: List[float], b: List[float]) -> List[float]:
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ]

def create_rotation_matrix(from_vec: List[float], to_vec: List[float]) -> List[List[float]]:
    from_vec = normalize(from_vec)
    to_vec = normalize(to_vec)

    rot_matrix = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ]

    cos_theta = dot_product(from_vec, to_vec)

    if abs(cos_theta - 1.0) < 1e-6:
        return rot_matrix

    if abs(cos_theta + 1.0) < 1e-6:
        if abs(from_vec[0]) < 0.9:
            perp = [1.0, 0.0, 0.0]
        else:
            perp = [0.0, 1.0, 0.0]

        axis = normalize(cross_product(from_vec, perp))

        rot_matrix = [
            [2 * axis[0] * axis[0] - 1, 2 * axis[0] * axis[1],       2 * axis[0] * axis[2]],
            [2 * axis[1] * axis[0],       2 * axis[1] * axis[1] - 1, 2 * axis[1] * axis[2]],
            [2 * axis[2] * axis[0],       2 * axis[2] * axis[1],       2 * axis[2] * axis[2] - 1]
        ]
        return rot_matrix

    axis = normalize(cross_product(from_vec, to_vec))
    sin_theta = math.sqrt(1 - cos_theta * cos_theta)
    one_minus_cos = 1 - cos_theta

    x, y, z = axis
    rot_matrix = [
        [cos_theta + x*x*one_minus_cos,
         x*y*one_minus_cos - z*sin_theta,
         x*z*one_minus_cos + y*sin_theta],

        [y*x*one_minus_cos + z*sin_theta,
         cos_theta + y*y*one_minus_cos,
         y*z*one_minus_cos - x*sin_theta],

        [z*x*one_minus_cos - y*sin_theta,
         z*y*one_minus_cos + x*sin_theta,
         cos_theta + z*z*one_minus_cos]
    ]

    return rot_matrix

def rotate_atom(original: Atom, rot_matrix: List[List[float]], pivot: Atom) -> Atom:
    rotated = Atom(
        atomID=original.atomID,
        atomName=original.atomName,
        residueName=original.residueName,
        resID=original.resID,
        x=0.0, y=0.0, z=0.0,
        occupancy=original.occupancy,
        tempFactor=original.tempFactor,
        chainID=original.chainID
    )

    temp_x = original.x - pivot.x
    temp_y = original.y - pivot.y
    temp_z = original.z - pivot.z

    rotated.x = rot_matrix[0][0] * temp_x + rot_matrix[0][1] * temp_y + rot_matrix[0][2] * temp_z
    rotated.y = rot_matrix[1][0] * temp_x + rot_matrix[1][1] * temp_y + rot_matrix[1][2] * temp_z
    rotated.z = rot_matrix[2][0] * temp_x + rot_matrix[2][1] * temp_y + rot_matrix[2][2] * temp_z

    rotated.x += pivot.x
    rotated.y += pivot.y
    rotated.z += pivot.z

    return rotated

def group_atoms_by_residue(atoms: List[Atom]) -> Dict[int, List[Atom]]:
    residue_groups = {}
    for atom in atoms:
        if atom.resID not in residue_groups:
            residue_groups[atom.resID] = []
        residue_groups[atom.resID].append(atom)
    return residue_groups

def find_root_atom_for_residue(atoms: List[Atom]) -> Atom:

    for atom in atoms:
        if 'P' in atom.atomName:
            return atom

    if atoms:
        return atoms[0]

    return Atom(0, "", "", 0, 0.0, 0.0, 0.0, 0.0, 0.0, "")

def find_tail_ends_for_residue(atoms: List[Atom], root: Atom, num_tails: int = 2) -> List[Atom]:
    distances = [(dist(atom, root), idx) for idx, atom in enumerate(atoms)]
    distances.sort(reverse=True, key=lambda x: x[0])
    tail_ends = [atoms[idx] for _, idx in distances[:min(num_tails, len(distances))]]
    return tail_ends

def calculate_residue_direction(atoms: List[Atom], root: Atom) -> List[float]:
    tail_ends = find_tail_ends_for_residue(atoms, root)
    avg_direction = [0.0, 0.0, 0.0]
    for tail in tail_ends:
        dir_vec = direction(root, tail)
        avg_direction[0] += dir_vec[0]
        avg_direction[1] += dir_vec[1]
        avg_direction[2] += dir_vec[2]

    avg_direction = [d / len(tail_ends) for d in avg_direction]

    return avg_direction

def reposition_residue(atoms: List[Atom],root: Atom,original_direction: List[float],new_direction: List[float]) -> List[Atom]:
    repositioned_atoms = []
    rot_matrix = create_rotation_matrix(original_direction, new_direction)
    for atom in atoms:
        rotated = rotate_atom(atom, rot_matrix, root)
        repositioned_atoms.append(rotated)
    return repositioned_atoms

def reposition_lipids_by_residue( atoms: List[Atom],residue_directions: Dict[int, List[float]]) -> List[Atom]:
    residue_groups = group_atoms_by_residue(atoms)
    all_repositioned_atoms = []
    
    for res_id, residue_atoms in residue_groups.items():
        print(f"\nProcessing residue {res_id} with {len(residue_atoms)} atoms")

        root = find_root_atom_for_residue(residue_atoms)
        original_direction = calculate_residue_direction(residue_atoms, root)

        if res_id in residue_directions:
            new_direction = residue_directions[res_id]
        else:
            new_direction = [0.0, 0.0, 1.0]
            print(f"  Warning: No direction specified for residue {res_id}, using default [0,0,1]")
        
        print(f"  Root atom: {root.atomName} at ({root.x:.2f}, {root.y:.2f}, {root.z:.2f})")
        print(f"  Original direction: [{original_direction[0]:.3f}, {original_direction[1]:.3f}, {original_direction[2]:.3f}]")
        print(f"  New direction: [{new_direction[0]:.3f}, {new_direction[1]:.3f}, {new_direction[2]:.3f}]")

        repositioned_residue = reposition_residue(
            residue_atoms, root, original_direction, new_direction
        )
        
        all_repositioned_atoms.extend(repositioned_residue)
    
    return all_repositioned_atoms

def read_pdb_file(filename: str) -> List[Atom]:
    atoms = []
    try:
        with open(filename, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    try:
                        atom_id = int(line[6:11].strip())
                        atom_name = line[12:16].strip()
                        residue_name = line[17:20].strip()
                        res_id = int(line[22:26].strip())
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        occupancy = float(line[54:60].strip())
                        temp_factor = float(line[60:66].strip())
                        chain_id = line[72:76].strip()

                        single_atom = Atom(
                            atomID=atom_id,
                            atomName=atom_name,
                            residueName=residue_name,
                            resID=res_id,
                            x=x,
                            y=y,
                            z=z,
                            occupancy=occupancy,
                            tempFactor=temp_factor,
                            chainID=chain_id
                        )

                        atoms.append(single_atom)
                    except ValueError as e:
                        print(f"Skipping malformed line: {line.strip()} - Error: {e}")
    except FileNotFoundError:
        print(f"Error: Could not open file {filename}")
        return []

    print(f"Read {len(atoms)} atoms from {filename}")
    return atoms

def write_pdb_file(filename: str, atoms: List[Atom], title: str = "REPOSITIONED LIPIDS"):
    try:
        with open(filename, 'w') as file:
            # Header
            file.write("CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")
            file.write("REMARK Generated by residue-based lipid repositioning program\n")
            file.write(f"REMARK {title}\n")

            # Atoms
            for a in atoms:
                line = (
                    f"ATOM  {a.atomID:5d} {a.atomName:>4} {a.residueName:>3} 1{a.resID:4d}    "
                    f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{a.occupancy:6.2f}{a.tempFactor:6.2f}      {a.chainID}\n"
                )
                file.write(line)

            file.write("END\n")
        print(f"Written {len(atoms)} atoms to {filename}")
    except IOError:
        print(f"Error: Could not create file {filename}")

def main():
    input_file = "lipid.pdb"
    original_atoms = read_pdb_file(input_file)

    residue_directions = {
        1: [1.0, 0.0, 0.0], 
        2: [1.0, 0.0, 0.0], 
        3: [1.0, 0.0, 0.0], 
        4: [1.0, 0.0, 0.0], 
        5: [1.0, 0.0, 0.0], 
        6: [1.0, 0.0, 0.0], 
    }

    repositioned_atoms = reposition_lipids_by_residue(original_atoms, residue_directions)
    write_pdb_file("lipids_by_residue.pdb", repositioned_atoms, "LIPIDS WITH RESIDUE-SPECIFIC ORIENTATIONS")
    

if __name__ == "__main__":
    main()