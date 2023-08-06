# Van der Wall (VdW) Forces Simulation
## Physics

We can calculate the potential (for VdW forces) energy Epot between two Argon atoms using the following Formula:

    Epot = 4 * ε ((σ/r)^12 - (σ/r)^6)

where ε and σ are constants depending on the atom (for Argon: ε=0.238eV, σ=3.405Å=340.5pm) and r is the distance between the atoms

    Epot = 4 * ε ((σ/r)^12 - (σ/r)^6)
    <=> Epot = (4 * ε * σ^12) / r^12 - (4 * ε * σ^6) / r^6

since ε and σ we can calculate the values for (4 * ε * σ^12) and (4 * ε * σ^6), let's call them ψ12 and ψ6

	ψ12 = (4 * ε * σ^12)
	ψ6 = (4 * ε * σ^6)

So now we have the following formula to calculate Epot between two Argon atoms:
    
    Epot = ψ12/r^12 - ψ6/r^6)

Using the following formula:
	
	1/2 * mass * v^2 = Ekin                | * 2 / mass
    <=> v^2 = Ekin*2 / mass                | sqrt()
    <=> v = sqrt(abs(Ekin)*2 / mass)

we can calculate the change in velocity Δv between two argon atoms using this formula:

    Δv = sqrt(abs(ΔEkin)*2 / mass)

and since ΔEkin = -ΔEpot, we can calculate it using:

    ΔEkin = -(Epot2 - Epot1)

Where Epot1 is the last calculated Epot between the two atoms and Epot2 is the now newly calculated Epot between the two atoms.
Therefore we get the formula:

    Δv = sqrt(abs(Epot2-Epot1)*2 / mass)

Now we only need to recalculate the recalculate the potential Energy for every timestep, update the velocity be Δv and the position by the time_passed*velocity

## Implementation

Currently the simulation only simulates the VdW forces on the x-axis.
I am focusing on fixing on axis before switching on the y-axis as well.

### Constants
const int pos_scale = 1e+2;
 - a factor applied to all position values for the atoms to work against floating point underflow
const long double ε = 0.238;
 - defined previously
const long double σ = 3.405*1e+2*pos_scale;
 - defined previously
const long double mass = 6.63352088e-7;
 - the mass given in 1e+19 * Kg

long double Θ6;
 - unused (was used previously for calcluation of ΔEpot)
long double Θ12;
 - unused (was used previously for calcluation of ΔEpot)
long double ψ6;
 - used for calculation of Epot (defined previously)
long double ψ12;
 - used for calculation of Epot (defined previously)

typedef struct {
    long double x;
    long double y;
} Vec2;
typedef struct {
    Vec2 pos;
    Vec2 velocity;
} Atom;

int atomcount = 0;
 - amount of atoms (will be set on initialization of atoms)

Atom* atoms = NULL;
 - array with all the atoms

long double** atom_pairs = NULL;
 - array containing the distance r between atom i and j

int workers = 3;
 - amount of threads for multithreading (unused)
long long int intervals_duration;
 - amount of intervals that will be simulated (interval=timestep)
long long int ipf;
 - intervals per frame (one frame marks the point where positions of all atoms are captured and written to the output file)
float progress;
 - progress from 0.0 to 1.0 of simulation generation
long double simulation_speed = 1e-12;
 - would not recommend changing this; time factor -> one second replay time is one picosecond (ps) in the simulation

bool mute = true;
 - should progress be displayed; change this at will

float replay_sim_scale = 3.0;
 - only used for simulation replay, impacts the scale of the simulation in replay 
int sim_width = 1080;
 - fixed values for replay (do not change these)
int sim_height = 720;
 - fixed values for replay (do not change these)
    
### Methods
void init_atoms(int _atomcount)
 - int _atomcount: amount of atoms to be spawned
Will spawn _atomcount atoms at the middle of the screen (when in replay) with a set distance between each other.
The velocity for the first two atoms will be set as well.
The function works as intended. It will be called by the init(...) function.
However, it may be useful to change some values to check the behaviour of the atoms.
Additionally, I would recomend starting with two atoms when trying to figure out why the simulation destabilizes.

void init(int _atomcount, float duration, long double interval, int fps)
 - int _atomcount: amount of atoms to be spawned
 - float duration: seconds to be simulated in (simulation_speed * duration) seconds
 - long double interval: timestep between each calculation in (simulation_speed * interval) seconds (recommended: simulation_speed*interval = 1fs)
 - int fps: how often per second (not factoring in simulation_speed) the positions of the atoms should be captured (values of 60-100 are recommended, if you don't plan on slowing the simulation down in replay)
Will set all required parameters for simulation.
The function works as intended. It will be called by the generate(...) function.
    
long double calc_Epot(long double r)
 - long double r: distance between the two atoms
returns: the potential energy Epot for the distance r
Will calculate Epot for the given distance r

long double calc_Δv(long double Epot1, long double Epot2, long double Δr, long double r)
 - long double Epot1: Epot calculated last timestep
 - lon  g double Epot2: Epot calculated this timestep
 - long dobule Δr: change in distance between the atoms
 - long double r: distance between the atoms
returns: Δv in pm/ns (i hope)
Will calculate Δv using the values given as exlained in the previous section.
If I had to guess, either this function or the calc_Epot(...) function must be behaving wrongly.
    

Vec2 compute_interaction(int id1, int id2)
 - int id1: index for atom 1 in array atoms
 - int id2: index for atom 2 in array atoms
returns: a vector describing Δv (Δv_vec)
    
void compute_atom(int i)
 - int i: index of atom in array atoms
Calculates the change in velocity for the atom and updates the velocity respectively.

void update(long double time_passed)
 - long double time_passed: time since last calculation (will usually just be the interval set if you don't do something real freaky)
will call compute_atom(i) for each atom and update the positions according to the time_passed and the velocity of the atom

void simulate(char* filepath, long double duration, long double interval, int atomcount, bool random_start, int fps)
 - char* filepath: where to save the simulation
 - float duration: seconds to be simulated in (simulation_speed * duration) seconds
 - long double interval: timestep between each calculation in (simulation_speed * interval) seconds (recommended: simulation_speed*interval = 1fs)
 - int atomcount: amount of atoms to simulate (i would recommend getting two atoms to work correctly before scaling up)
 - bool random_start: unused
 - int fps: how often per second (not factoring in simulation_speed) the positions of the atoms should be captured (values of 60-100 are recommended, if you don't plan on slowing the simulation down in replay)
generates a simulation 
    
void generate_test_data(char* filepath, int type)
 - char* filepath: where to save the tests results
 - int type: which test to run
generate different test data (usually just maps out a function for it to be viewed as a graph later) which can be viewed by running fastsim_testvis.py (might need to change what tests will be evaluated in the code)

### run.sh
If you have the gcc compiler installed, running the shell script run.sh should automatically recompile and run fastsim.c

### To play a simulation

To relay a simulation you can use the simulator.py script in the parent folder.

```bash
python3 ../simulator.py -l -f /simulations/Csim5
```

Run the script with the arguments -l (load) -f (filepath) [filepath]

You can find more info on the simulator.py script in the README.md (it is broken apart from replay)

#### Help Page
```
$ python simulator.py --help
```
```
usage: simulator.py [-h] [-generate] [-load] [--filepath FILEPATH] [--duration DURATION] [--interval INTERVAL] [--molecule_amount MOLECULE_AMOUNT] [--sim_radius SIM_RADIUS] [--scale_factor SCALE_FACTOR] [--random_start]

MOLECULAR SIMULATOR |  generate and simulate Van Der Waals Forces between particles

options:
-h, --help            show this help message and exit
-generate, -g         generate a new simulation
-load, -l             load a generated simulation
--filepath FILEPATH, -f FILEPATH
filepath to save the generated simulation to or to load the simulation from
--duration DURATION, -d DURATION
duration [in seconds] of the simulation [generate only]
--interval INTERVAL, -i INTERVAL
time between recalculations [in seconds] of velocity for molecules [generate only]
--molecule_amount MOLECULE_AMOUNT, -a MOLECULE_AMOUNT
amount of molecules [generate only]
--sim_radius SIM_RADIUS, -sr SIM_RADIUS
max simulation radius regarding each molecule [generate only]
--scale_factor SCALE_FACTOR, -sf SCALE_FACTOR
higher values reduce simulation space [generate only]
--random_start, -r    use random starting positions for molecules [generate only]
--simulation_speed SIMULATION_SPEED, -sp SIMULATION_SPEED
how fast the simulation should be displayed [load only]
```
