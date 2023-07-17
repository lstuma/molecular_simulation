#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

// constants needed for energy calculation
ε = -997.1;
σ = 3.401;

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


void init_atoms(_atomcount)
{
    // set atomcount to actual value
    atomcount = _atomcount
    
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

void compute_atom(i)
{
    for(int j = i+1; j < atomcount; j++)
    {
        Vec2 force = compute_interaction(atoms[i].pos, atoms[j].pos)  

        atoms[i].velocity.x += force.x
        atoms[i].velocity.y += force.y
        
        atoms[j].velocity.x -= force.x
        atoms[j].velocity.y -= force.y
    }
}


Vec2 compute_interaction(Vec2 pos1, Vec2 pos2)
{
    // calc position offset between atoms
    Vec2 offset;
    offset.x = pos1.x-pos2.x;
    offset.y = pos1.y-pos2.y;

    // calc direction
    Vec2 dir;
    dir.x = offset-x/(fabsf(offset.x)+1e-6);
    dir.x = offset-x/(fabsf(offset.x)+1e-6);
    
    // calc distance
    float distance = sqrt(offset.x*offset.x + offset.y*offset.y);
    
    // calc energy
    float energy = 4 * ε * powf(σ / (distance+1e-20), 12) - powf(σ / (distance+1e-20), 6));
    
    Vec2 force;
    force.x = dir.x * energy;
    force.y = dir.y * energy;
    return force;
}

void generate(char* filepath, float duration, float interval, int molecule_amount, float sim_radius, bool random_start)
{
    printf("%s\n", filepath);
    return;
}

int main(int argc, char** kwargs)
{
    printf("Hello World!\n");
    generate("/home/sdoxl/programming/projects/molecular_simulation/simuations/C_sim1", 10.0, 0.001, 10, 20, true);
    return 0;
}