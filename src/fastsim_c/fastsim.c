#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

// constants needed for energy calculation
const double ε = -997.1;
const double σ = 3.401;

const double mass = 39.948    // mass
const double Θ6 = 48*ε*powf(σ,6);
const double Θ12 = 48*ε*powf(σ,12):


typedef struct {
    float x;
    float y;
} Vec2;

typedef struct {
    Vec2 pos;
    Vec2 velocity;
} Atom;

// amount of atoms -> will be set on init_atoms()
int atomcount = 0;

// array containing all atoms
Atom* atoms = NULL;

double** atom_pairs = NULL;

// amount of threads
int workers = 3;


void init_atoms(int _atomcount)
{
    // set atomcount to actual value
    atomcount = _atomcount;
    
    // initialize array with atoms
    atoms = (Atom*)malloc(_atomcount * sizeof(Atom));
    
    // set positions etc. for atoms
    for(int i = 0; i < _atomcount; i++) {
        atoms[i].pos.x = 0.0;
        atoms[i].pos.y = 0.0;
        atoms[i].velocity.x = 0.0;
        atoms[i].velocity.y = 0.0;
    }

    // initialize array with atom distances
    atom_pairs = (double**)malloc(_atomcount*sizeof(double*)+_atomcount*sizeof(double));

    for(int i = 0; i < _atomcount; i++)
        for(int j = 0; j < _atomcount; j++)
            atom_pairs[i][j] = 0.0;


    return;
}

Vec2 compute_interaction(int id1, int id2)
{
    // calc position offset between atoms
    Vec2 offset;
    offset.x = atoms[id1].pos.x-atoms[id2].pos.x;
    offset.y = atoms[id1].pos.y-atoms[id2].pos.y;

    // calc direction
    Vec2 dir;
    dir.x = offset.x/(fabsf(offset.x)+1e-6);
    dir.x = offset.x/(fabsf(offset.x)+1e-6);
    
    // calc distance
    double r = sqrt(offset.x*offset.x + offset.y*offset.y);
    atom_pairs[id1][id2] = r;
    double r6 = powf(r, 6);


    // calc energy
    double ΔEpot = (Θ6 / (r6*r)) - (Θ12 / (r6*r6));
    double acceleration = (-ΔEpot/mass);

    Vec2 acceleration_vec;
    acceleration_vec.x = dir.x * ΔEpot;
    acceleration_vec.y = dir.y * ΔEpot;
    return acceleration_vec;
}


void compute_atom(int i)
{
    for(int j = i+1; j < atomcount; j++)
    {
        double Δr = atom_pairs[i][j];

        Vec2 Δv = compute_interaction(i, j);

        Δr -= atom_pairs[i][j];

        // calculate velocity change
        Δv.x *= Δr;
        Δv.y *= Δr;

        atoms[i].velocity.x += Δv.x;
        atoms[i].velocity.y += Δv.y;

        atoms[j].velocity.x -= Δv.x;
        atoms[j].velocity.y -= Δv.y;
    }
}


void update(float time_passed)
{
    reset_velocity_all();

    for(int i = 0; i < atomcount; i++)
    {
        compute_atom(i);
        atoms[i].pos.x += atoms[i].velocity.x * time_passed;
        atoms[i].pos.y += atoms[i].velocity.y * time_passed;
    }
}


void generate(char* filepath, float duration, float interval, int atomcount, bool random_start, int fps)
{
    // initialize atom array
    init_atoms(atomcount);

    // amount of intervals that will be run
    int intervals_duration = (int)(duration/interval);
    // intervals per frame
    int ipf = (int)((1.0/fps)/interval);
    if(ipf<=10) printf("\n\033[31mERROR \033[0m| <= 1 IPF (HIGHLY UNSTABLE) -> PLEASE CHOOSE LOWER INTERVAL!");

    // progress (0% to 100%)
    float progress = 0.0;

    for(int intervals_passed = 0; intervals_passed < intervals_duration; intervals_passed++)
    {
        if(intervals_passed%ipf==0)
        {
            //write(atoms);
        }

        update(interval);

        progress = intervals_passed/intervals_duration;

        intervals_passed += 1;
        printf("%d, ", intervals_passed);
    }
    printf("\nDONE!");

    return;
}

int main(int argc, char** kwargs)
{
    printf("\nHello World!");
    generate("\n/home/sdoxl/programming/projects/molecular_simulation/simuations/C_sim1", 10.0, 0.001, 20, true, 60);
    return 0;
}