from logger import log
from multiprocessing import Pool, Array
import numpy as np

ε: float = -997.1       # emin
σ: float = 3.401        # rmin
mass: float = 39.948    # mass
Θ6 = 48*ε*(σ**6)
Θ12 = 48*ε*(σ**12)

atomcount = 0
atoms = None
shared_atoms = None
atom_ids = None

# list of distances between atoms
atom_pairs = None
shared_atom_pairs = None

# tweak as necessary (number of processes)
workers = 3

total_energy = 0.0
start_total_energy = None
total_energy_diff = None


def init_atoms(atoms_list):
    global atom_ids, atomcount, atoms, atom_pairs, shared_atoms, shared_atom_pairs
    
    atomcount = len(atoms_list)
    shape = (atomcount, 2, 2)
    atom_ids = np.array([i for i in range(atomcount)])
    atoms = Array('d', shape[0]*shape[1]*shape[2])
    shared_atoms = np.frombuffer(atoms.get_obj(), dtype=np.float64).reshape(shape)

    shape = (atomcount, atomcount)
    atom_pairs = Array('d', shape[0]*shape[1])
    shared_atom_pairs = np.frombuffer(atom_pairs.get_obj(), dtype=np.float64).reshape(shape)



def reset_velocity_all():
    global shared_atoms

    # log msg cuz this function is no longer used
    log("resetting global velocity", "error")

    for atom in shared_atoms:
        atom[1][:] = 0.0


def calc_total_energy():
    global total_energy, mass, shared_atoms

    total_energy = 0.0

    for i, atom in enumerate(shared_atoms):
        for other_atom in shared_atoms[i+1:]:
            total_energy += calc_Epot(np.linalg.norm(atom[0] - other_atom[0])+1e-16)
        total_energy += atom[1].sum()*mass

    return total_energy



def calc_atom_velocity(id):
    global atomcount, shared_atoms, shared_atom_pairs

    for other_id in range(id+1, atomcount):
        acceleration, r = calc_acceleration_lazy(shared_atoms[id][0], shared_atoms[other_id][0])
        Δr = shared_atom_pairs[id, other_id] - r
        shared_atom_pairs[id, other_id] = r
        # add force to velocity
        shared_atoms[id][1] += acceleration * Δr
        shared_atoms[other_id][1] -= acceleration  * Δr

def calc_acceleration(pos1: np.array, pos2: np.array):
    global ε, σ

    offset = pos1-pos2
    dir = np.array((offset/(abs(offset)+1e-6)))
    r = np.linalg.norm(pos1 - pos2)+1e-16
    tmp = (σ/r)
    tmp6 = tmp**6
    tmp12 = tmp6*tmp6

    ΔEpot = 4*ε * ( - 12*(tmp12 / r) + 12*(tmp6 / r) )
    acceleration = (-ΔEpot/mass) * dir
    return acceleration, r

def calc_acceleration_lazy(pos1: np.array, pos2: np.array):
    global ε, σ

    offset = pos1-pos2
    dir = np.array((offset/(abs(offset)+1e-6)))
    r = np.linalg.norm(pos1 - pos2)+1e-16
    r6 = r**6
    r7 = r6*r
    r13 = r6*r6

    ΔEpot = (Θ6 / r7) - (Θ12 / r13)
    acceleration = (-ΔEpot/mass) * dir
    return acceleration, r

def calc_Epot(r):
    global ε, σ
    if r == 0.0: r = 1.0
    Epot = -ε * ((σ/r)**12 - 2 * (σ/r)**6)
    return Epot


#
# normal update
#
def update(time_passed):
    global atom_ids, shared_atoms
    
    # update atom positions
    for atom in shared_atoms:
        atom[0] += atom[1] * time_passed


    # calc forces using multiprocessing
    for i in range(atomcount):
        calc_atom_velocity(i)

    return

def slow_update():
    global atom_ids, shared_atoms, start_total_energy, total_energy, total_energy_diff

    if not start_total_energy:
        # calc total energy
        start_total_energy = calc_total_energy()

    # get difference in energy
    total_energy_diff = calc_total_energy() / start_total_energy

    # set total energy back to start total energy
    reversed_total_energy_diff = start_total_energy/total_energy
    for atom in shared_atoms:
        atom[1] *= reversed_total_energy_diff


#
# leapfrog
#
global lf_step
def lf_init(step):
    global lf_step
    lf_step = step

    # velocity half step update
    for id in range(atomcount):
        calc_atom_velocity(id)

def lf_init_parallel(step):
    lf_init(step)

    global pool, workers
    pool = Pool(workers)

def lf_close_parallel():
    global pool
    pool.close()
    pool.join()


def lf_update():
    global atomcount, shared_atoms, lf_step

    # position full step update
    for atom in shared_atoms:
        atom[0] += atom[1] * lf_step

    # velocity full step update
    for id in range(atomcount):
        calc_atom_velocity(id)

def lf_update_parallel():
    global atomcount, shared_atoms, lf_step, workers

    # position full step update
    for atom in shared_atoms:
        atom[0] += atom[1] * lf_step

    # velocity full step update
    global pool
    pool.map(calc_atom_velocity, atom_ids)