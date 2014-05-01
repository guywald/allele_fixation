Allele Fixation - Matlab Project

Run Start.m file.

Description:
Matlab implementation for a university project on allele fixation simulation.
This project simulates the effect of a usefull allele on its fixation in individuals in a population.

In each turn of the simulation individuals move and reproduce. Each individual can move an amount of steps accordings to it's genotype. The simulation shows how these alleles fixate according to different settings.

With different set-ups the simulation can show fixation models such as Genetic Drift.

The simulation is built from two parts:
1) Individual Based Model: Controlled by the Simulation Control panel
2) Analytical model: Controlled by the Analytical Model Control panel.

How to operate:

Simulation control  - 
 - Number of steps: (The Phenotype for each genotype) defines the amount of movments an individual can move in each turn of the game. Example: If AA has 3 it can move each turn 3 steps.

 - Number of creatures: The amount of creates in a simulation. The number stays fixed during a simulation (When Parents reproduce their offsprings (one male, one female) replace them.

 - Initial Genotypes: The Initial allele distribution between the individuals.

 - Dimension of board: The size of the board. For example: 5 means a 5x5 board.

 - Rep. Age: Reproduction age of the individuals. Each individual starts with the age of 1 and each turn it increments. Only on Rep. Age they can reproduce.

 - Draw every: Frequency of when to draw the board. 1 - Every turn, 2 - Every two turns. This affects the speed performance.

 - Draw Board: Check if you want to draw the board or not at all. Affects simulation speed performance.

 - Run: Start / Continue the simulation.

 - Step: Pauses the simulation, each additional click advances the game in one step (not a turn. A turn can be multiple steps according to the value in Number of steps).

 - Reset: Clears the simulation.


Analytical Model Control:

The analytical model calculates the Frequency of an allele spanning generations according to the following formula:

Delta(q) = q_(t+1)-q_t = (pq[p(w_(Aa) - w_(AA)) + q(w_(aa) - w_(Aa))])/(p^2 x w_(AA)+ 2pq x w_(Ab) + q^2 x w_(aa))

Delta(q) describes the change of allele 'a' frequency along the generations.
The code uses numeric integraion Ode23 to calculate this.

Allele 'a' is q in the formula, 'A' is p which is (1-q).

 - Fitness: Fitness values for the analytical model.

 - 'a' frequency: The initial frequency of the 'a' allele.

 - Number of generations: The number of generations for the formula.

Enjoy,
Guy
