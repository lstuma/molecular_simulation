import visualizer as gfx
from random import randrange
import math
from logger import log

molecule_amount = 100
molecule_sim_radius = 300.0
scale_multiplier = 15.0


def init(randomize=True):
    global molecules
    molecules = []
    global r_min
    r_min = 3.401 # KJ/mol
    global E_min
    E_min = -1.654 # Ångström
    global molecular_mass
    molecular_mass = 39.948 # u
    global kinetic_energy, total_velocity, last_kinetic_energy, last_total_velocity
    kinetic_energy, total_velocity, last_kinetic_energy, last_total_velocity = None, None, 0.0, 0.0
    global molecule_sim_radius
    global total_time_passed
    total_time_passed = 0.0
    global scale_multiplier
    global molecule_amount
    global sim_width, sim_height
    sim_width, sim_height = 1080, 720

    # Adjust molecule_sim_radius
    sim_width /= scale_multiplier
    sim_height /= scale_multiplier
    
    gfx.init()
    
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
    def __init__(self, x=0, y=0, velocity=(0,0), randomize=False):
        if randomize:
            self._x, self._y = randrange(1,1070//scale_multiplier), randrange(1,710//scale_multiplier)
        else:
            self._x, self._y = x, y
        
        self.velocity_x, self.velocity_y = velocity
        
        # Create graphical counterpiece
        self.graphic = gfx.Circle(pos=(self._x * scale_multiplier, self._y * scale_multiplier), color=(randrange(50, 255), randrange(50, 255), randrange(50, 255)), radius=10*(scale_multiplier/15))
        
    
    @property
    def energy(self):
        ...
    
    @property
    def x(self):
        return self._x
    
    @x.setter
    def x(self, var):
        global scale_multiplier, sim_width
        self._x = var % sim_width
        self.graphic.shape.x = self._x * scale_multiplier
    
    @property
    def y(self):
        return self._y
    
    @y.setter
    def y(self, var):
        global scale_multiplier, sim_height
        self._y = var % sim_height
        self.graphic.shape.y = self._y * scale_multiplier

    def calc_direction(self, offset):
        total_offset = sum(abs(i) for i in offset)
        if total_offset != 0: return -offset[0]/total_offset, -offset[1]/total_offset
        else: return randrange(-1, 1)/2, randrange(-1, 1)/2
        
    def calc_offset(self, other):
        offset = self.x - other.x, self.y - other.y
        return offset

    def calc_force(self, distance):
        global r_min
        global E_min
        if distance != 0:
            force = E_min * ((r_min/distance)**12 - 2 * (r_min/distance)**6)
            return force
        # molecules are on same position, this shouldn't happen but cannot be avoided with time leaps
        else: return 0
    
    def calc_distance(self, offset):
        return math.hypot(offset[0], offset[1])
    
    def calc_velocity(self):
        global molecules, molecule_sim_radius
        self.velocity_x, self.velocity_y = 0, 0
        for molecule in molecules:
            if not molecule == self:
                offset = self.calc_offset(molecule)
                distance = self.calc_distance(offset)
                if distance <= molecule_sim_radius:
                    direction = self.calc_direction(offset)
                    force = self.calc_force(distance)
                    #log("calc force: " + str("%.6f" % force))
                    self.velocity_x += force * direction[0]
                    self.velocity_y += force * direction[1]
    
    def __sub__(self, other):
        return self.calc_distance(other=other)
    

# fixed update will allways be called at about the interval, but might slightly differ in ranges up to 1e-4
def fixed_update(time_passed, mute=False):
    global total_time_passed
    total_time_passed += time_passed
    
    for molecule in molecules:
        molecule.calc_velocity()
        molecule.x += molecule.velocity_x * time_passed
        molecule.y += molecule.velocity_y * time_passed
    
    if not mute: log("total velocity: " + str("%.6f" % calc_total_velocity()) + "\t|\tkinetic energy: " + str("%.6f" % calc_kinetic_energy()) + "\t|\ttime passed: " + str("%.2f" % total_time_passed), end="", begin="\r")
    
def slow_fixed_update(time_passed):
    global last_total_velocity, total_velocity
    velocity_loss = last_total_velocity - total_velocity
    last_total_velocity = total_velocity
    if abs(velocity_loss) >= 10: log("kinetic difference: " + str("%.6f" % velocity_loss), "warning" if abs(velocity_loss) <= 50 else "warning2", begin="\n")
    

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
    

def run(fixed_callback=fixed_update, interval=0.2, slow_fixed_callback=slow_fixed_update, slow_interval=1):
    gfx.run(fixed_callback=fixed_callback, interval=interval, slow_fixed_callback=slow_fixed_callback, slow_interval=slow_interval)

def main(interval=0.001, slow_interval=1, randomize=True):
    log("initializig...")
    init(randomize)
    log("finished initializing", "success")

    run(fixed_callback=fixed_update, interval=interval, slow_fixed_callback=slow_fixed_update, slow_interval=slow_interval)


if __name__ == '__main__':
    main()