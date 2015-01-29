# Design of the BilliardModels package

## File structure
The simulator, which calculates billiard trajectories, is in BilliardModels.jl. The visualization is in BilliardVisualization.jl.

## Data structures
Firstly a billiard table should be constructed. This consists of a collection of Obstacles.

A particle may then be constructed on the billiard table. The particle and billiard table are passed to the billiard_dynamics function, which generates a trajectory by following consecutive collisions with obstacles.

