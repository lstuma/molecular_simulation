from logger import log
from multiprocessing import Pool
import numpy as np

ε: float = -997.1       # emin
σ: float = 3.401        # rmin
mass: float = 39.948    # mass

atomcount = 0
atoms = None
atom_ids = None


# tweak as necessary (number of processes)
workers = 3

total_energy = 0.0
start_total_energy = None
total_energy_diff = None

def init_atoms(atoms_list):
    global atom_ids, atomcount, atoms
    
    atomcount = len(atoms_list)
    atoms = np.array((atoms_list))
    atom_ids = np.array([i for i in range(atomcount)])


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
        total_energy += atom.velocity.sum()/mass

    return total_energy



def calc_atom_velocity(id, time_passed):
    global atomcount, atoms

    for other_id in range(id+1, atomcount):
        acceleration, distance = calc_acceleration(atoms[id].pos, atoms[other_id].pos)
        # add force to velocity
        atoms[id].velocity += acceleration * time_passed
        atoms[other_id].velocity -= acceleration  * time_passed

def calc_acceleration(pos1: np.array, pos2: np.array):
    global ε, σ

    offset = pos1-pos2
    dir = np.array((offset/(abs(offset)+1e-6)))
    distance = np.linalg.norm(pos1 - pos2)+1e-16
    tmp = (σ/distance)
    tmp6 = tmp**6
    tmp12 = tmp6*tmp6

    ΔEpot = 4*ε * ( - 12*(tmp12 / distance) + 12*(tmp6 / distance) )
    acceleration = (-ΔEpot/mass) * dir
    return acceleration, distance

def calc_acceleration_lazy(pos1: np.array, pos2: np.array):
    global ε, σ

    offset = pos1-pos2
    dir = np.array((offset/(abs(offset)+1e-6)))
    distance = np.linalg.norm(pos1 - pos2)+1e-16
    tmp = (σ/distance)
    tmp6 = tmp**6
    tmp12 = tmp6*tmp6

    ΔEpot = 4*ε * ( - 12*(tmp12 / distance) + 12*(tmp6 / distance) )
    acceleration = (-ΔE1pot/mass) * dir
    return acceleration, distance

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
    
    # calc forces using multiprocessing
    for i in range(atomcount):
        calc_atom_velocity(i, time_passed)

    # update atom positions
    for atom in atoms:
        atom.update_pos(atom.velocity * time_passed)
    return

def slow_update():
    global atom_ids, atoms, start_total_energy, total_energy, total_energy_diff

    if not start_total_energy:
        # calc total energy
        start_total_energy = calc_total_energy()

    # get difference in energy
    total_energy_diff = start_total_energy / calc_total_energy()

    # set total energy back to start total energy
    for atom in atoms:
        atom.velocity *= total_energy_diff

#
# leapfrog
#
global lf_step;
def lf_init(step):
    global lf_step
    lf_step = step

    # velocity half step update
    for id in range(atomcount):
        calc_atom_velocity(id, step/2.0)

def lf_update(step):
    global atomcount, atoms

    # velocity full step update
    for id in range(atomcount):
        calc_atom_velocity(id, step/2.0)


    # position full step update
    for atom in atoms:
        atom.update_pos(atom.velocity * time_passed)