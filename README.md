# Molecular Simulation

## Usage
### Using GUI
```
$ python gui.py
```
[![https://imgur.com/EBflVbH.png](https://imgur.com/EBflVbH.png)](https://imgur.com/EBflVbH.png)
[![https://imgur.com/kYgUC8Y.png](https://imgur.com/kYgUC8Y.png)](https://imgur.com/kYgUC8Y.png)
### Shell Usage
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
##### Generating A Simulation:
```
$ python3 simulator.py --generate --random_start --filepath /tmp/simulation --duration 10 --interval 0.0001 --molecule_amount 150
```
```
INFO		| simulating 150 molecules for 10.0 seconds

INFO        | c     100000/100000        [PROGRESS] ▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰ 100.0 %      0:00:01 left
SUCCESS		| DONE!
SUCCESS		| simulation written to /tmp/simulation
```

##### Loading A Simulation
```
$ python simulator.py --load --filepath /tmp/simulation
```
```
INFO		| loading simulation at /tmp/sim2

INFO		| molecule amount: 500.0
INFO		| interval: 0.001s
INFO		| length of simulation: 20.0s
INFO        | w     total time passed: 4.3645s      203.7920 FPS    1000.0000 CPS   0.9215 simulation speed         |       4359 intervals passed
```
[![https://imgur.com/i6OOEyM.png](https://imgur.com/i6OOEyM.png)](https://imgur.com/i6OOEyM.png)
