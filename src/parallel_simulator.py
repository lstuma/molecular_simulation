from logger import log
from multiprocessing import Pool
import numpy as np

ε: float = -997.1       # emin
σ: float = 3.401        # rmin
mass: float = 39.948    # mass
Θ6 = 48*ε*(σ**6)
Θ12 = 48*ε*(σ**12)

atomcount = 0
atoms = None
atom_ids = None
# list of distances between atoms
atom_pairs = None


# tweak as necessary (number of processes)
workers = 3

total_energy = 0.0
start_total_energy = None
total_energy_diff = None

def init_atoms(atoms_list):
    global atom_ids, atomcount, atoms, atom_pairs
    
    atomcount = len(atoms_list)
    atoms = np.array((atoms_list))
    atom_ids = np.array([i for i in range(atomcount)])

    atom_pairs = np.zeros((atomcount, atomcount))


def reset_velocity_all():
    global atoms

    # log msg cuz this function is no longer used
    log("resetting global velocity", "error")

    for atom in atoms:
        atom.velocity[:] = 0.0


def calc_total_energy():
    global total_energy, mass, atoms

    total_energy = 0.0

    for i, atom in enumerate(atoms):
        for other_atom in atoms[i+1:]:
            total_energy += calc_Epot(np.linalg.norm(atom.pos - other_atom.pos)+1e-16)
        total_energy += atom.velocity.sum()*mass

    return total_energy



def calc_atom_velocity(id):
    global atomcount, atoms, atom_pairs

    for other_id in range(id+1, atomcount):
        acceleration, r = calc_acceleration_lazy(atoms[id].pos, atoms[other_id].pos)
        Δr = atom_pairs[id, other_id] - r
        atom_pairs[id, other_id] = r
        # add force to velocity
        atoms[id].velocity += acceleration * Δr
        atoms[other_id].velocity -= acceleration  * Δr

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
    global atom_ids, atoms
    
    # update atom positions
    for atom in atoms:
        atom.update_pos(atom.velocity * time_passed)


    # calc forces using multiprocessing
    for i in range(atomcount):
        calc_atom_velocity(i)

    return

def slow_update():
    global atom_ids, atoms, start_total_energy, total_energy, total_energy_diff

    if not start_total_energy:
        # calc total energy
        start_total_energy = calc_total_energy()

    # get difference in energy
    total_energy_diff = calc_total_energy() / start_total_energy

    # set total energy back to start total energy
    reversed_total_energy_diff = start_total_energy/total_energy
    for atom in atoms:
        atom.velocity *= reversed_total_energy_diff

#
# leapfrog
#
global lf_step;
def lf_init(step):
    global lf_step
    lf_step = step

    # velocity half step update
    for id in range(atomcount):
        calc_atom_velocity(id)

def lf_update():
    global atomcount, atoms, lf_step

    # position full step update
    for atom in atoms:
        atom.update_pos(atom.velocity * lf_step)

    # velocity full step update
    for id in range(atomcount):
        calc_atom_velocity(id)