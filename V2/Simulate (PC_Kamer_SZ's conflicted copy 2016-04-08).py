from Multi import Multi
from PERMSimulation import PERMSimulation
import time
import os
import sys
import numpy as np
import numba as jit
import warnings as warnings
np.errstate(divide="ignore")
np.errstate(over="ignore")
np.errstate(invalid="ignore")
warnings.filterwarnings("ignore")
    
# Define the functions which do the worker

Var = np.array([
        250,   # Max polymer size
        0.2,   # Temperature
        0.25,  # Epsilon
        0.8,   # Sigma
        0.0,   # Force
        6,     # Theta
        10000, # Polymer count
        ]);
Add = np.zeros(np.size(Var));  Add[4]  = 0.2;
Stop = np.zeros(np.size(Var)); Stop[4] = 1;

# Now we have to create the worker function
def Worker(qm, qw, i):
    UnitId = i;
    print("Unit Id = "+str(i));
    Working = 0;
    while(1):
        Data = qm.get();
        if (not Data[1] == UnitId):
            qm.put((Data[0], Data[1]));
            time.sleep(0.1);
            continue;
        if (Data == "Done"):
            return;
        Var = Data[0];
    
        testSim = PERMSimulation()
        testSim.initialise(int(Var[6]),int(Var[6]))
        testSim.T = Var[1];
        testSim.F = Var[4];
        print("Unit #"+str(UnitId)+"starting simulation for T = "+str(testSim.T)+" & F = "+str(testSim.F));
        testSim.run_simulation()
        
        i = 1;
        file = str(os.getcwd())+"/Data/Save_"+str(i)+".npz"
        while (os.path.isfile(file)):
            file = str(os.getcwd())+"/Data/Save_"+str(i)+".npz"
            i+=1
        testSim.save_results(file)
        qw.put(UnitId);
        
    
if __name__ == '__main__':
    Multi(Worker, Var, Add, Stop, True);