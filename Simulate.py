from Multi import Multi
from Perm import Perm
import time
import numpy as np
import numba as jit
    
# Define the functions which do the worker

Var = np.array([
        250,   # Max polymer size
        1,     # Temperature
        0.25,  # Epsilon
        0.8,   # Sigma
        0,     # Force
        6,     # Theta
        10000, # Polymer count
        ]);
Add = np.zeros(np.size(Var));  Add[4]  = 0.01; # Add 0.01 per step to force
Stop = np.zeros(np.size(Var)); Stop[4] = 0.01;  # Stop when we reach a force of 0.5

# Now we have to create the worker function
def Worker(qm, qw):
    Working = 0;
    while(1):
        Data = qm.get();
        if (Data == "Done"):
            return;
        Var = Data[0];
    
        testSim = Perm()

        testSim.initialise(10000,10000, Data[4])
        %time testSim.run_simulation()
        
        qw.put(int(Data[1]));
    
if __name__ == '__main__':
    Multi(Worker, Var, Add, Stop, True);