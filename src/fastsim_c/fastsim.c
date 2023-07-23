#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

// constants needed for energy calculation
const int ε = -997.1;
const int σ = 3.401;

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

    return;
}

void reset_velocity_all()
{
    // go through all atoms and set velocity to zero for later recalculation of forces
    for(int i = 0; i < atomcount; i++)
    {
        atoms[i].velocity.x = 0;
        atoms[i].velocity.y = 0;
    }

    return;
}


Vec2 compute_interaction(Vec2 pos1, Vec2 pos2)
{
    // calc position offset between atoms
    Vec2 offset;
    offset.x = pos1.x-pos2.x;
    offset.y = pos1.y-pos2.y;

    // calc direction
    Vec2 dir;
    dir.x = offset.x/(fabsf(offset.x)+1e-6);
    dir.x = offset.x/(fabsf(offset.x)+1e-6);
    
    // calc distance
    float distance = sqrt(offset.x*offset.x + offset.y*offset.y);
    
    // calc energy
    float energy = 4 * ε * (powf(σ / (distance+1e-20), 12) - powf(σ / (distance+1e-20), 6));
    
    Vec2 force;
    force.x = dir.x * energy;
    force.y = dir.y * energy;
    return force;
}


void compute_atom(int i)
{
    for(int j = i+1; j < atomcount; j++)
    {
    Vec2 force = compute_interaction(atoms[i].pos, atoms[j].pos);

        atoms[i].velocity.x += force.x;
        atoms[i].velocity.y += force.y;

        atoms[j].velocity.x -= force.x;
        atoms[j].velocity.y -= force.y;
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