import visualizer as gfx
from random import randrange, choice
import math
from logger import log, log_up, log_down
import logger
import argparse
import time, datetime
import numpy as np

import parallel_simulator as fast_simulator

# extra output?
verbose = False

# Wether to allow lag when fps is low or compromise accuracy instead
_fixed_deltatime = False

molecule_amount = 20
molecule_sim_radius = 300.0
scale_multiplier = 15.0

_interval, _slow_interval = None, None
simulation_speed, cps, fps = 1.0, 0.0, 0.0

# how fast the simulation should run
simulation_speed_factor = 1
simulation_fps_save = 100

def init(randomize=True, initialize_gfx=True):
    global molecules
    molecules = []
    global r
    #r = 3.401 # Ångström
    r = 3.401
    global ε
    ε = 1.654 # KJ/mol
    ε = -997.1
    global molecular_mass
    molecular_mass = 39.948 # u
    global kinetic_energy, total_velocity, last_kinetic_energy, last_total_velocity
    kinetic_energy, total_velocity, last_kinetic_energy, last_total_velocity = None, None, 0.0, 0.0
    global molecule_sim_radius
    global total_time_passed
    total_time_passed = 0.
    global scale_multiplier
    global molecule_amount
    global sim_width, sim_height
    sim_width, sim_height = 1080, 720

    # Adjust molecule_sim_radius
    sim_width /= scale_multiplier
    sim_height /= scale_multiplier
    
    if initialize_gfx: gfx.init()
    
    if randomize:
        for i in range(molecule_amount):
            molecules.append(Molecule(randomize=randomize))
    else:
        side_length = round(math.sqrt(molecule_amount))
        offset = 5
        step_x, step_y = (1070/scale_multiplier) // side_length, (710/scale_multiplier) // side_length
        for i in range(side_length):
            for j in range(side_length):
                molecules.append(Molecule(x=offset + step_x*i, y = offset + step_y*j))
    
class Molecule:
    def __init__(self, x=0.0, y=0.0, velocity=(0.0,0.0), randomize=False):
        if randomize:
            self.pos = np.array((float(randrange(1,int(1070//scale_multiplier))), float(randrange(1,int(710//scale_multiplier)))))
        else:
            self.pos = np.array((float(x), float(y)))
        
        self.velocity = np.array((float(velocity[0]), float(velocity[1])))
        self.velocity_x, self.velocity_y = self.velocity[0], self.velocity[1]
        
        # Create graphical counterpiece
        self.graphic = gfx.Circle(pos=(self.pos * scale_multiplier), color=gfx.color('rainbow'), radius=10*(scale_multiplier/15))
    
    @property
    def energy(self):
        ...
    
    @property
    def x(self):
        return self.pos[0]

    def update_pos(self, var):
        self.pos += var
        self.graphic.shape.position += var
    
    @x.setter
    def x(self, var):
        global scale_multiplier, sim_width
        self.pos[0] = var % sim_width
        self.graphic.shape.x = self.pos[0] * scale_multiplier
    
    @property
    def y(self):
        return self.pos[1]
    
    @y.setter
    def y(self, var):
        global scale_multiplier, sim_height
        self.pos[1] = var % sim_height
        self.graphic.shape.y = self.pos[1] * scale_multiplier

    def calc_direction(self, offset):
        total_offset = sum(abs(i) for i in offset)
        if total_offset != 0: return -offset[0]/total_offset, -offset[1]/total_offset
        else: return randrange(-1, 1)/2, randrange(-1, 1)/2
        
    def calc_offset(self, other):
        offset = self.x - other.x, self.y - other.y
        return offset

    def calc_force(self, σ):
        global r
        global ε
        if σ != 0:
            force = -ε * ((σ/r)**12 - 2 * (σ/r)**6)
            return force
        # molecules are on same position, this shouldn't happen but cannot be avoided with time leaps
        else: return 0
    
    def calc_distance(self, offset):
        return math.hypot(offset[0], offset[1])
    
    def calc_velocity(self):
        global molecules, molecule_sim_radius
        self.velocity[:] = 0
        for molecule in molecules:
            if not molecule == self:
                offset = self.calc_offset(molecule)
                distance = self.calc_distance(offset)
                if distance <= molecule_sim_radius:
                    direction = self.calc_direction(offset)
                    force = self.calc_force(distance)
                    self.velocity[0] += force * direction[0]
                    self.velocity[1] += force * direction[1]
    
    def __sub__(self, other):
        return self.calc_distance(other=other)
    

# update will allways be called at about the interval, but might slightly differ depending on performance. setting fixed_deltatime to True will slow down simulation speed to call the update at fixed intervals relative to the simulation
def update(time_passed, mute=False):
    global total_time_passed, _fixed_deltatime, _interval, cps, simulation_speed
    if time_passed >= 0.1: simulation_speed = _interval/time_passed
    if _fixed_deltatime: time_passed = _interval
    cps = 1/time_passed
    total_time_passed += time_passed
    
    for molecule in molecules:
        molecule.calc_velocity()
        molecule.x += molecule.velocity_x * time_passed
        molecule.y += molecule.velocity_y * time_passed
    
    if not mute:
        log("total velocity: " + "%.6f" % calc_total_velocity())
        log("kinetic energy: " + "%.6f" % calc_kinetic_energy())
        log("time passed: " + "%.2f" % total_time_passed)
        log_up(2)
    
def slow_update(time_passed):
    global last_total_velocity, total_velocity, _fixed_deltatime, _slow_interval, cps, simulation_speed
    if _fixed_deltatime: time_passed = _slow_interval
    velocity_loss = last_total_velocity - total_velocity
    last_total_velocity = total_velocity
    log_down(3)
    log("calculations per second: " + "%.6f" % cps, "info" if cps >= (1/_interval)-10 else "warning")
    log("simulation speed: " + "%.6f" % simulation_speed, "info" if simulation_speed >= 0.9 else "warning")
    log("kinetic difference: " + "%.6f" % velocity_loss, "info" if abs(velocity_loss) <= 50 else "warning")
    log_up(6)
    

def calc_total_velocity():
    global total_velocity
    total_velocity = 0
    for molecule in molecules:
        total_velocity += abs(molecule.velocity_x) + abs(molecule.velocity_y)
    return total_velocity


def calc_kinetic_energy():
    global kinetic_energy, total_velocity, molecular_mass
    kinetic_energy = total_velocity * molecular_mass
    return kinetic_energy

def load_molecule_positions(fp, current_pos, next_pos):
    distance = next_pos - current_pos - 1
    for i in range(distance): fp.readline()
    return fp.readline().split(",")[:-1]

def from_file_update(time_passed):
    global molecules, total_time_passed, _interval, intervals_length, _duration, last_intervals_passed, intervals_passed, fps, cps, fp, simulation_speed, simulation_speed_factor, verbose
    time_passed *= simulation_speed_factor
    if total_time_passed <= _duration: total_time_passed += time_passed
    fps = 1/time_passed
    if intervals_passed % 500 == 50:
        simulation_speed = intervals_passed*_interval/total_time_passed
    log("    " + "total time passed: " + "%.4f" % total_time_passed + "s")
    log("%.4f" % fps + " FPS\t" + "%.4f" % cps + " CPS\t")
    log("~%.4f" % simulation_speed + " simulation speed")
    log(str(intervals_passed) + " intervals passed")
    log_up(3)
    if intervals_passed < intervals_length:
        # reading
        log(logger.colors['info']+"[r]"+logger.colors['reset'], end="\r")
        positions = load_molecule_positions(fp, last_intervals_passed, intervals_passed)
        last_intervals_passed = intervals_passed
        intervals_passed = int(total_time_passed // _interval)
        for molecule in molecules:
            position = [float(i) for i in positions.pop(0).split("x")]
            # writing
            log(logger.colors['success']+"[w]"+logger.colors['reset'], end="\r")
            molecule.x = position[0]
            molecule.y = position[1]

def simulate_from_file(filepath):
    global fp
    fp = open(filepath, "r")

    global molecule_amount, molecule_sim_radius, scale_multiplier, _duration, intervals_length, simulation_speed_factor, verbose
    log("loading simulation at \033[;32m" + filepath, end="\n\n")
    molecule_amount = float(fp.readline()[:-1])
    interval = float(fp.readline()[:-1])
    _duration = float(fp.readline()[:-1])
    molecule_sim_radius = float(fp.readline()[:-1])
    scale_multiplier = float(fp.readline()[:-1])
    intervals_length = int(_duration//interval)

    log("molecule amount: " + str(molecule_amount))
    log("interval: " + str(interval) + "s " + "(%.1f Hz)" % (1/interval))
    log("duration of simulation: " + str(_duration) + "s")
    log("simulation speed: " + "%.2f" % simulation_speed_factor, end="\n\n")

    global _interval, _slow_interval, _fixed_deltatime, cps, last_intervals_passed
    last_intervals_passed = 0
    _interval = interval
    _fixed_deltatime = True
    cps = 1/_interval

    init(False)
    gfx.run(fixed_callback=from_file_update, interval=_interval)


def simulate_to_file(filepath, _duration=100, fixed_callback=update, interval=0.2, slow_fixed_callback=None, slow_interval=None):
    if slow_fixed_callback != None or slow_interval != None:
        log("slow_fixed_callback and slow_interval are unused in simulations written to files!", "warning")

    global _interval, _slow_interval, _slow_interval, _fixed_deltatime, verbose
    _interval, _slow_interval = interval, slow_interval
    _fixed_deltatime = True
    intervals_length = int(_duration/interval)

    real_start_time = time.time()

    with open(filepath, "w") as fp:
        global molecules, molecule_amount, molecule_sim_radius, scale_multiplier, simulation_fps_save

        # calculating properties to get real_ftps as well as intervals_per_frame
        seconds_per_frame = 1.0/simulation_fps_save
        intervals_per_frame, _ = divmod(seconds_per_frame, interval)
        intervals_per_frame = int(intervals_per_frame)
        # fps that can be hit with set interval
        real_fps = 1.0/(intervals_per_frame*interval)

        # hours:minutes:seconds to display how long the simulation will be in log
        _duration_minutes, _duration_seconds = divmod(_duration, 60)
        _duration_hours, _duration_minutes = divmod(_duration_minutes, 60)

        # details about simulation
        log("simulating " + str(molecule_amount) + " molecules for " + "%d:%02d:%02d" % (_duration_hours, _duration_minutes, _duration_seconds) + " (" + str(intervals_length) + " intervals)")
        log("simulating at " + str(real_fps) + " FPS")
        log("creating (max) " + str(_duration*real_fps) + " x " + str(molecule_amount) + " datapoints", end="\n\n")

        # speeeeeeeed
        fast_simulator.init_atoms(molecules)

        # leapfrog
        fast_simulator.lf_init_parallel(interval)

        intervals_passed = 1
        time_passed = 1e-10
        fp.write(str(molecule_amount)+"\n"+str(seconds_per_frame)+"\n"+str(_duration)+"\n"+str(molecule_sim_radius)+"\n"+str(scale_multiplier)+"\n")

        # estimate
        seconds_left, minutes_left, hours_left, real_time_passed = 0, 0, 0, 0

        while time_passed <= _duration:
            if intervals_passed%100 == 1:
                progress = (intervals_passed/intervals_length)*100
                _progress = progress/100

                # estimate
                real_time_passed = time.time() - real_start_time
                seconds_left = int((real_time_passed/(_progress))*(1-(_progress)))
                minutes_left, seconds_left = divmod(seconds_left, 60)
                hours_left, minutes_left = divmod(minutes_left, 60)

            if intervals_passed%intervals_per_frame == 1:
                for atom in fast_simulator.shared_atoms:
                    fp.write(str(atom[0][0])+"x"+str(atom[0][1])+",")
                fp.write("\n")


            # calculate new pos + velocity
            fast_simulator.lf_update_parallel()


            # every 100 intervals slow_update is called
            if intervals_passed%100 == 0:
                fast_simulator.slow_update()
                if verbose:
                    # display progress
                    log(str(intervals_passed) + "/" + str(intervals_length))
                    log("[PROGRESS] \033[;34m" + ("▰"*int(progress/2)) + "\033[;32m" + ("▱"*(50-int(progress/2))) + "\033[0;0m " + "%.1f" % progress + " %")
                    log("%d:%02d:%02d" % (hours_left, minutes_left, seconds_left) + " left     ")
                    log('[START] %.4f KJ/mol' % fast_simulator.start_total_energy, level='warning')
                    log('[NOW] %.4f KJ/mol' % fast_simulator.total_energy, level='warning')
                    log('energy diff: ' + '%.4f' % fast_simulator.total_energy_diff + ' (' + choice(['/', '|', '\\', '-']) + ')')
                    log_up(5)

            time_passed += interval
            intervals_passed += 1

        fast_simulator.lf_close_parallel()
        log("DONE!", "success", begin="\n\n")
        log("simulation written to \033[;32m" + filepath, "success")

def run(fixed_callback=update, interval=0.2, slow_fixed_callback=slow_update, slow_interval=1, fixed_deltatime=False):
    global _interval, _slow_interval, _fixed_deltatime
    _interval, _slow_interval = interval, slow_interval
    _fixed_deltatime = fixed_deltatime
    gfx.run(fixed_callback=fixed_callback, interval=interval, slow_fixed_callback=slow_fixed_callback, slow_interval=slow_interval)

def main(interval=0.001, slow_interval=3, randomize=True, fixed_deltatime=False):
    log("starting with deltatime "+ ("fixed" if fixed_deltatime else "loose"))
    log("initializig...")
    init(randomize)
    log("finished initializing", "success", end="\n\n")

    run(fixed_callback=update, interval=interval, slow_fixed_callback=slow_update, slow_interval=slow_interval, fixed_deltatime=fixed_deltatime)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="\033[;32mMOLECULAR SIMULATOR \033[;37m|\t\033[0;0m generate and simulate Van Der Waals Forces between particles")


    parser.add_argument("--generate", "-g", action="store_true", help="generate a new simulation")
    parser.add_argument("--load", "-l", action="store_true", help="load a generated simulation")

    parser.add_argument("--filepath", "-f", dest="filepath", help="filepath to save the generated simulation to or to load the simulation from")
    parser.add_argument("--duration", "-d", dest="duration", help="duration [in seconds] of the simulation \033[;37m[\033[;33mgenerate only\033[;37m]\033[0;0m", default=20.0, type=float)
    parser.add_argument("--interval", "-i", dest="interval", help="time between recalculations [in femtoseconds] of velocity for molecules \033[;37m[\033[;33mgenerate only\033[;37m]\033[0;0m", default=0.01, type=float)
    parser.add_argument("--molecule_amount", "-a", dest="molecule_amount", help="amount of molecules \033[;37m[\033[;33mgenerate only\033[;37m]\033[0;0m", default= 50, type=int)
    parser.add_argument("--sim_radius", "-sr", dest="sim_radius", help="max simulation radius regarding each molecule \033[;37m[\033[;33mgenerate only\033[;37m]\033[0;0m", default=300, type=float)
    parser.add_argument("--scale_factor", "-sf", dest="scale_factor", help="higher values reduce simulation space \033[;37m[\033[;33mgenerate only\033[;37m]\033[0;0m", default=15.0, type=float)
    parser.add_argument("--random_start", "-r", dest="random_start", help="use random starting positions for molecules \033[;37m[\033[;33mgenerate only\033[;37m]\033[0;0m", action="store_true")
    parser.add_argument("--simulation_speed", "-sp", dest="simulation_speed", help="how fast the simulation should be displayed \033[;37m[\033[;33mload only\033[;37m]\033[0;0m", default=1.0, type=float)
    parser.add_argument("--frames_per_second", "-fps", dest="fps", help="in how many fps the simulation should be saved \033[;37m[\033[;33mload only\033[;37m]\033[0;0m", default=100, type=int)
    parser.add_argument("--verbose", "-v", dest="verbose", help="output extra info in console", action="store_true")


    args = parser.parse_args()

    if args.generate and args.load:
        log("please select either \033[;32mgenerate\033[0;0m or \033[;32mload\033[0;0m", "error")
    if not args.filepath:
        log("filepath required", "error")

    verbose = args.verbose

    if args.generate:
        molecule_amount = args.molecule_amount
        molecule_sim_radius = args.sim_radius
        scale_multiplier = args.scale_factor
        simulation_fps_save = args.fps
        init(args.random_start, initialize_gfx=False)
        simulate_to_file(filepath=args.filepath, _duration=args.duration, interval=args.interval*1e-15)
    elif args.load:
        intervals_passed = 0
        simulation_speed_factor = args.simulation_speed
        simulate_from_file(args.filepath)

