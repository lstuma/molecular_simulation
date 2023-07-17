from multiprocessing import Pool
import numpy as np

ε: float = -997.1
σ: float = 3.401

atomcount = 0
atoms = None
atom_ids = None
dir = np.zeros((2))

# tweak as necessary (number of processes)
workers = 3

def init_atoms(atoms_list):
    global atom_ids, atomcount, atoms
    
    atomcount = len(atoms_list)
    atoms = np.array((atoms_list))
    atom_ids = np.array([i for i in range(atomcount)])
    
    
def multiupdate(time_passed):
    global atom_ids, atoms
    
    # create multiprocessing pool
    pool = Pool(workers)
    
    # calc forces using multiprocessing
    forces = pool.map(compute_atom, atom_ids)
    
    # update atom positions
    for atom in atoms:
        atom.update_pos(atom.velocity * time_passed)
        
    pool.close()
    pool.join()

    return

def reset_velocity_all(id):
    global atoms
    
    for atom in atoms:
        atom.velocity[:] = 0.0

def compute_atom(id):
    global atomcount, atoms, compute_func
    
    for other_id in range(id+1, atomcount):
        force = compute_func(atoms[id].pos, atoms[other_id].pos)
        atoms[id].velocity += force
        atoms[other_id].velocity -= force
        
def compute_interaction(pos1: np.array, pos2: np.array):
    global ε, σ, dir
    
    offset = pos1-pos2
    dir[:] = offset/(abs(offset)+1e-6)
    distance = np.linalg.norm(pos1 - pos2)
    energy = 4 * ε * ((σ / (distance+1e-20)) ** 12 - (σ / (distance+1e-20)) ** 6)
    return energy * dir

def compute_interaction_lazy(pos1: np.array, pos2: np.array):
    global ε, σ, dir

    offset = pos1-pos2
    dir[:] = offset/(abs(offset)+1e-6)
    distance = np.linalg.norm(pos1 - pos2)
    energy = 4 * ε * ((distance / (σ+1e-20)) ** 12 - (distance / (σ+1e-20)) ** 6)
    return energy * dir


compute_func = compute_interaction_lazy


