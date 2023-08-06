#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>

/*
*   CLARIFICATIONS:
*    - simulation speed will be stretched so that 1 second (s) equals one picosecond (ps)
*    - distance will be given in milli Ångström (mÅ)
*
 */

// constants needed for energy calculation
// epsilon and sigma found at https://personal.utdallas.edu/~bxw109120/exercises/CHEM6V19_E_LJ17_MD.pdf
const int pos_scale = 1e+2;
const long double ε = 0.238;         // unit in eV
const long double σ = 3.405*1e+2*pos_scale;         // in centi Ångström (cÅ) or 1e-12m (pm)
const long double mass = 6.63352088e-7;     // mass given in 1e+19 * Kg

long double Θ6;
long double Θ12;
long double ψ6;
long double ψ12;

typedef struct {
    long double x;
    long double y;
} Vec2;

typedef struct {
    Vec2 pos;       // should be given in cÅ (pm)
    Vec2 velocity;  // should be given in cÅ/s (pm/s)
} Atom;

// amount of atoms -> will be set on init_atoms()
int atomcount = 0;

// array containing all atoms
Atom* atoms = NULL;

long double** atom_pairs = NULL;

// amount of threads
int workers = 3;
// amount of intervals that will be simulated and interalvs per frame (ipf)
long long int intervals_duration; long long int ipf;
float progress;
long double simulation_speed = 1e-12;

// for logging
bool mute = true;

// for simulation replay
float replay_sim_scale = 3.0;
int sim_width = 1080;
int sim_height = 720;

void init_atoms(int _atomcount)
{
    // set atomcount to actual value
    atomcount = _atomcount;
    
    // initialize array with atoms
    atoms = (Atom*)malloc(_atomcount * sizeof(Atom));
    
    // set positions etc. for atoms
    for(int i = 0; i < _atomcount; i++) {
        atoms[i].pos.x = 0.5 * sim_width/replay_sim_scale + i * 350;
        atoms[i].pos.y = 0.5 * sim_height/replay_sim_scale;
        atoms[i].velocity.x = 0.0;
        atoms[i].velocity.y = 0.0;
    }
    atoms[0].velocity.x = 100;
    atoms[1].velocity.x = -atoms[0].velocity.x;

    // initialize array with atom distances
    atom_pairs = (long double**)malloc(_atomcount*sizeof(long double*));

    for(int i = 0; i < _atomcount; i++) {
        atom_pairs[i] = (long double*)malloc(_atomcount*sizeof(long double));
        for(int j = 0; j < _atomcount; j++) {
            // calc position offset between atoms
            Vec2 offset;
            offset.x = (atoms[i].pos.x-atoms[j].pos.x)/pos_scale;
            offset.y = (atoms[i].pos.y-atoms[j].pos.y)/pos_scale;

            // calc distance
            atom_pairs[i][j] = sqrt(offset.x*offset.x + offset.y*offset.y);
          }
    }

    return;
}

void init(int _atomcount, float duration, long double interval, int fps)
{
    // amount of intervals that will be run
    intervals_duration = (long long int)(duration/interval);
    printf("\033[34mINFO  \033[37m|\033[0;30m SIMULATING %lld INTERVALS\n", intervals_duration);

    // intervals per frame
    ipf = (long long  int)((1.0/fps)/interval);
    if(ipf<=10) {
        printf("\033[31mERROR \033[37m|\033[0;30m >=10 IPF (=%lld) (HIGHLY UNSTABLE) -> PLEASE CHOOSE LOWER INTERVAL!\n", ipf);
        printf("\033[35mEXIT  \033[37m|\033[0;30m EXITING!\n"); exit(1);
    } else {
        printf("\033[34mINFO  \033[37m|\033[0;30m <=10 IPF (=%lld) STABLE :)\n", ipf);
    }
    // progress (0% to 100%)
    progress = 0.0;

    // set beautiful constants
    Θ6 = 48*ε*powl(σ,6);
    Θ12 = 48*ε*powl(σ,12);
    ψ6 = 4*ε*powl(σ,6);
    ψ12 = 4*ε*powl(σ,12);

    // initialize atoms
    init_atoms(_atomcount);

    return;
}

long double calc_Epot(long double r) {
    long double r6 = powl(r, 6);
    long double r12 = r6*r6;
    return ψ12/r12 - ψ6/r6;
}

long double calc_Δv(long double Epot1, long double Epot2, long double Δr, long double r) {
    // calc energy difference to previous position
    long double ΔEpot = Epot1 - Epot2;
    ΔEpot *= 0.5; // we only want ΔEpot for one atom not for both
    long double ΔEkin = 1.60218 * ΔEpot; // ΔEkin given in 1e+19 * J

    // 1 if we are getting closer to σ and -1 if we are getting further away from σ
    int dir = (r<=σ && Δr>=0 || r>=σ && Δr<=0)?1:-1;

    long double Δv = sqrtf(2.0*fabsl(ΔEkin)/mass)*1e+12*dir; // Δv given in pm/s
    if(Δv!=0.0)    printf("%d : %.20Lf, %.20Lf, %.20Lf, ", dir, Δv, Δr, r);
    return Δv * pos_scale;
}

Vec2 compute_interaction(int id1, int id2)
{
    // calc position offset between atoms
    Vec2 offset;
  offset.x = (atoms[id2].pos.x-atoms[id1].pos.x)/pos_scale;
    offset.y = (atoms[id2].pos.y-atoms[id1].pos.y)/pos_scale;

    // change in distance
    long double Δr = atom_pairs[id1][id2];

    // calc direction
    Vec2 dir;
    dir.x = offset.x/(fabsl(offset.x)+1e-6);
    dir.y = offset.y/(fabsl(offset.y)+1e-6);
    
    // calc distance
    long double Epot1 = calc_Epot(atom_pairs[id1][id2]);
    atom_pairs[id1][id2] = sqrtl(offset.x*offset.x + offset.y*offset.y);
    long double Epot2 = calc_Epot(atom_pairs[id1][id2]);

    // calc Δr
    Δr = atom_pairs[id1][id2] - Δr;

    long double Δv = calc_Δv(Epot1, Epot2, Δr, atom_pairs[id1][id2]);


    Vec2 Δv_vec;
    Δv_vec.x = dir.x * Δv;
    Δv_vec.y = dir.y * Δv;

    return Δv_vec;
}


void compute_atom(int i)
{
    for(int j = i+1; j < atomcount; j++)
    {
        Vec2 Δv = compute_interaction(i, j);

        atoms[i].velocity.x += Δv.x;
        atoms[i].velocity.y += Δv.y;

        atoms[j].velocity.x -= Δv.x;
        atoms[j].velocity.y -= Δv.y;
    }
}


void update(long double time_passed)
{
    for(int i = 0; i < atomcount; i++)
    {
        compute_atom(i);

        atoms[i].pos.x -= atoms[i].velocity.x * time_passed * pos_scale;
        //atoms[i].pos.y += atoms[i].velocity.y;
    }
}

void simulate(char* filepath, long double duration, long double interval, int atomcount, bool random_start, int fps)
{
    // initialize needed variables for simulation
    init(atomcount, duration, interval, fps);

    // open file
    FILE* fp = fopen(filepath, "w");

    interval *= simulation_speed;
    printf("\033[34mINFO  \033[37m|\033[0;30m SIMULATION SPEED AT: %.20Lf AND THEREFORE INTERVAL AT %.20Lf\n", simulation_speed, interval);

    // write metadata
    fprintf(fp, "%d\n%.20lf\n%.20Lf\n%.20lf\n%.20f\n",
        atomcount,
        1.0/fps,
        duration,
        1.0,
        replay_sim_scale);

    printf("\033[34mINFO  \033[37m|\033[0;30m SIMULATING NOW ...\n");
    for(int intervals_passed = 0; intervals_passed < intervals; intervals_passed++)
    {
        if(intervals_passed%ipf==0)
        {
            if(!mute) {
                // log progress
                progress = (long double)intervals_passed/(long double)intervals_duration;
                int log_progress = (int)(progress*25);
                printf("\033[34mINFO  \033[37m|\033[0;30m %d/%lld ", intervals_passed, intervals_duration);
                printf("\033[33m"); for(int i = 0; i < log_progress; i++) printf("#");
                printf("\033[37m"); for(int i = log_progress; i < 25; i++) printf(".");
                printf("\033[0;30m %.2f %c \r", progress, '%');
            }

            // write to file
            for(int i = 0; i < atomcount; i++)
                fprintf(fp, "%.20Lfx%.20Lf,", atoms[i].pos.x, atoms[i].pos.y);
            fprintf(fp, "\n");
        }

        update(interval);

        intervals_passed += 1;
    }
    fclose(fp);
    printf("\nDONE!");

    return;
}

void generate_test_data(char* filepath, int type)
{
    // set beautiful constants
    Θ6 = 48*ε*powl(σ,6);
    Θ12 = 48*ε*powl(σ,12);
    ψ6 = 4*ε*powl(σ,6);
    ψ12 = 4*ε*powl(σ,12);

    // open file
    FILE* fp = fopen(filepath, "w");

    switch(type) {
        case 0:
            // graph of calc_Epot(r)
            fprintf(fp, "Epot(r)\n");
            for(long double r = -10; r <= 10; r+=1e-2) {
                long double Epot = calc_Epot(r);
                fprintf(fp, "%.20Lf,%.20Lf\n", r, Epot>100?100:Epot);
            }
            break;

        case 1:
            // graph of calc_Δv(Epot(3.405), Epot(r)) -> will cut of extreme values
            fprintf(fp, "Δv(3.405,i)\n");
            for(long double r = -7; r <= 7; r+=1e-2) {
                long double Δv = calc_Δv(calc_Epot(3.405), calc_Epot(r), r-3.405, r);
                if(Δv > -2e+15)
                    fprintf(fp, "%.20Lf,%.20Lf\n", r, Δv);
            }
            break;

        case 2:
            // graph of calc_Δv(Epot(2.405), Epot(r)) -> will cut of extreme values
            fprintf(fp, "Δv(2.405,i)\n");
            for(long double r = -7; r <= 7; r+=1e-2) {
                long double Δv = calc_Δv(calc_Epot(2.405), calc_Epot(r), r-2.405, r);
                if(Δv > -2e+15)
                    fprintf(fp, "%.20Lf,%.20Lf\n", r, Δv);
            }
            break;

        case 3:
            // graph of calc_Δv(Epot(4.405), Epot(r)) -> will cut of extreme values
            fprintf(fp, "Δv(4.405,i)\n");
            for(long double r = -7; r <= 7; r+=1e-2) {
                long double Δv = calc_Δv(calc_Epot(4.405), calc_Epot(r), r-4.405, r);
                if(Δv > -2e+15)
                    fprintf(fp, "%.20Lf,%.20Lf\n", r, Δv);
            }
            break;

        case 4:
            // graph of calc_Δv(Epot(0), Epot(r)) -> will cut of extreme values
            fprintf(fp, "Δv(0,i)\n");
            for(long double r = -7; r <= 7; r+=1e-2) {
            long double Δv = calc_Δv(calc_Epot(0), calc_Epot(r), r, 0);
            if(Δv > -2e+15)
                fprintf(fp, "%.20Lf,%.20Lf\n", r, Δv);
            }
            break;

        case 5:
            // heatmap of calc_Δv -> will cut of extreme values
            for(long double i = 2; i <= 5; i+=1e-3) {
                for(long double j = 2; j <= 5; j+=1e-3) {
                    long double Δv = calc_Δv(calc_Epot(i), calc_Epot(j), j-i, j);
                    fprintf(fp, "%.20Lf,", Δv);
                }
                fprintf(fp, "\n");
            }
            break;
    }
    fclose(fp);
}

int main(int argc, char** kwargs)
{
    // test0
    generate_test_data("/home/lstuma/programming/projects/molecular_simulation/simulations/test0", 0);
    // test1.0
    generate_test_data("/home/lstuma/programming/projects/molecular_simulation/simulations/test1_0", 1);
    // test1.1
    generate_test_data("/home/lstuma/programming/projects/molecular_simulation/simulations/test1_1", 2);
    // test1.2
    generate_test_data("/home/lstuma/programming/projects/molecular_simulation/simulations/test1_2", 3);
    // test1.3
    generate_test_data("/home/lstuma/programming/projects/molecular_simulation/simulations/test1_3", 3);
    // test2
    //generate_test_data("/home/lstuma/programming/projects/molecular_simulation/simulations/test2", 5);

    // generate simulation
    simulate("/home/lstuma/programming/projects/molecular_simulation/simulations/Csim5", 10, 0.0001, 2, true, 60);
    return 0;
}
