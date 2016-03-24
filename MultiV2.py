# The goal of this script is to do the parallel processing of multiple polymers

# Import required libraries
from multiprocessing import Process, Queue, cpu_count
from numba import jit
from copy import deepcopy
import numpy as np
import math
import time
import random
import os
import msvcrt
import sys
import psutil

# Assign default variables
bead_position = np.zeros([2,2])       # Initialize array
bead_position[1,:] = [1,0]            # Set initial condition
polymer_size = 2                      # Amount of beads in the polymer
max_polymer_size = 250                # Maximum polymer size
polymer_count = 1E5                   # Amount of polymers that we are simulating

T = 5                                 # Temperature in kelvin
epsilon = 0.25                        # Epsilon voor de lennard Jones
sigma = 0.8                           # Sigma voor de Lennard Jones

theta = 6                             # Number of angles
theta_weight = np.zeros([theta,1])    # Weight of angles
theta_weight_sum = 0                  # Sum of the weights
polymer_weight = 1                    # Weight

multi_step_size = 3                   # Multiprocessor data size communication

# Create standard Polymer
BasePolymer = {
        'Weight':polymer_weight,      # Current weight of the polymer
        'Weight3':1,                  # Weight at the time of the first addition
        'Size':polymer_size,          # Current size of the polymer
        'Alive':True,                 # Is the polymer alive or not?
        'BeadPosition':bead_position  # Positions of the beads
    }
polarray = []

# Define the queue's
q  = Queue();                         # Queue to communicate the jobs
q2 = Queue();                         # Queue to communicate the results

# Define functions

def restart_program():
    """Restarts the current program
    """
    os.execv(sys.executable, [sys.executable] + sys.argv)

def workerFunction(q, q2): # This function controls the workers and defines what they should do
    while (1):
        dicti = q.get();
        if (dicti == 'Exit'):
            break;
        listi = dicti['List']
        
        try:
            len(listi)
        except:
            print("ListError - Self");
            continue;
            
        if (dicti['Task'] == 'addBead'):
            for i in range(0,len(listi)-1):
                listi[i] = addBead(listi[i]);
            q2.put(dicti);
        time.sleep(0.001);

@jit
def addBead(var): # This function adds a bead to a polymer
    # We start by checking the energy and thus the change of all theta directions
    theta_weight_sum = 0
    while (theta_weight_sum == 0):
        theta_offset = random.uniform(0,1)
        theta_weight = np.zeros([theta,1])
        for i in range(0,theta):
            E = 0
            new_bead_x = var['BeadPosition'][var['Size']-1, 0] + math.cos(2 * math.pi * (i+theta_offset) / theta);
            new_bead_y = var['BeadPosition'][var['Size']-1, 1] + math.sin(2 * math.pi * (i+theta_offset) / theta);
            for j in range(0,var['Size']):
                dx = var['BeadPosition'][j, 0] - new_bead_x
                dy = var['BeadPosition'][j, 1] - new_bead_y
                if (dx < 3 and dy < 3): # Cut-off distance
                    d = math.sqrt(dx**2+dy**2)
                    E += 4 * epsilon *((sigma / d)**12 - (sigma/d)**6)
            theta_weight[i] = math.exp(-E/T)
            theta_weight_sum += theta_weight[i]
    
    # Next we use a roulette weel algorith to determine which direction we go
    rand = random.uniform(0,1)
    var['Weight'] = var['Weight'] * theta_weight_sum
    theta_weight = np.divide(theta_weight, theta_weight_sum)
    for i in range(0,theta):
        rand -= theta_weight[i]
        if (rand < 0):
            direction = i
            new_bead = np.zeros([1,2])
            new_bead[0,0] = var['BeadPosition'][var['Size']-1,0] + math.cos(2 * math.pi * (i+theta_offset) / theta)
            new_bead[0,1] = var['BeadPosition'][var['Size']-1,1] + math.sin(2 * math.pi * (i+theta_offset) / theta)
            var['BeadPosition'] = np.append(var['BeadPosition'], new_bead,axis=0)
            var['Size'] += 1
            break
    if (var['Size'] == 3):
        var['Weight3'] = var['Weight']
    return var


# Master function - This function initializes everything and starts the workers
if __name__ == '__main__':
    # Create the workers
    cpuList = []
    for i in range(0,cpu_count()): #Cpu count utilizes your computer. use cpu_count()-1 to do other things.
        cpuList.append(Process(target=workerFunction, args=(q,q2,)))
        cpuList[i].start();
    while(1):
        polymer_size = 2;
        upLimInput = int(sys.argv[1]);
        print("Start");
        start_time = time.time();

        # Create a subfolder in Data were we are going to save everything
        i = 1;
        while (os.path.isdir(str(os.getcwd())+"/Data/"+str(sys.argv[2])+"_"+str(i))):
            i+=1;
        os.makedirs(str(os.getcwd())+"/Data/"+str(sys.argv[2])+"_"+str(i))
        folder = str(os.getcwd())+"/Data/"+str(sys.argv[2])+"_"+str(i)+"/";
        print("Folder = "+str(folder));

        # Create a list of all polymers
        polarray = []

        # Fill the list with N polymers
        N = 2
        for i in range(0,966):
            polarray.append(BasePolymer)

        # Now we start the actual process
        steptime = 0;
        polymer_start_size = polymer_size
        while (polymer_size < max_polymer_size):
            if (polymer_size > polymer_start_size):
                print("Iteration took: "+str((time.time() - steptime))+"\t\t Speed = "+(str(Amount / (time.time() - steptime))));
            steptime = time.time();
            if msvcrt.kbhit():
                if (msvcrt.getch() == 'q'):
                    break
                print ("Key pressed: "+ str(msvcrt.getch()))
            # Add a bead to every polymer. We send this task to the workers.
            c = 0;
            jobCount = 0;
            # Submit the tasks to the workers
            while( c < len(polarray)):
                if (len(polarray) - c < 0):
                    q.put({'List':polarray[c:], 'Id':c, 'Task':'addBead'})
                else:
                    q.put({'List':polarray[c:c+multi_step_size], 'Id':c, 'Task':'addBead'})
                jobCount += 1;
                c = c + multi_step_size;
            # Retrieve the tasks from the workers
            while (jobCount > 0):
                dicti = q2.get();
                if (len(polarray)-dicti['Id'] < 0):
                    polarray[dicti['Id']:] = dicti['List']
                else:
                        polarray[dicti['Id']:dicti['Id']+multi_step_size] = dicti['List']
                jobCount -= 1;

            # Next we need to prune the polymers
            # We do this in the main thread because it shouldn't take much time and saves overhead
            TotWeight = 0;
            Amount = len(polarray);
            for i in range(0,Amount):
                TotWeight += polarray[i]['Weight'];
            AvgWeight = TotWeight / Amount
            print("AVgWeight = "+str(AvgWeight)+" Amomunt = "+str(Amount))
            upLim = 2 ** (math.log(Amount,10)-upLimInput)
            lowLim = 1.5 ** (math.log(Amount,10)-1)
            for i in range(0,Amount):
                if (polarray[i]['Weight'] > upLim * AvgWeight / polarray[i]['Weight3']): # Branch
                    polarray[i]['Weight'] /= 2;
                    polarray.append(deepcopy(polarray[i]))
                if (polarray[i]['Weight'] < lowLim * AvgWeight / polarray[i]['Weight3']): # Prune
                    rand = random.uniform(0,1)
                    if (rand < 0.5):
                        polarray[i]['Weight'] *= 2
                    else:
                        polarray[i]['Alive'] = False;

            # Now we don't want the dead polymers to crowd the list,
            # so we are going to remove them and save them in a file.
            DeadList = [];
            i = Amount;
            while (i > 0):
                i -= 1;
                if (polarray[i]['Alive'] == False):
                    DeadList.append(polarray.pop(i))
            # Save the dead polymers to a file
            np.save(folder+"/Dead_"+str(polymer_size), DeadList);
            polymer_size += 1;
            print("PolymerSize = "+str(polymer_size)+" Living = "+str(len(polarray))+" Just died = "+str(len(DeadList)));
        print("Exiting");

        # Now we have to save the polymers to a data file
        np.save(folder+'/Survivors', polarray);
        print("Finished in: "+str(time.time() - start_time));
        print("So now we are done. Congratulations!");
        time.sleep(0.5);
        print(""); time.sleep(0.5);
        print(""); time.sleep(0.5);
        print(""); time.sleep(0.5);
        print(""); time.sleep(0.5);
        print(""); time.sleep(0.5);
        print(""); time.sleep(0.5);
        print(""); time.sleep(0.5);
        print("Just kidding");
        print(""); time.sleep(0.5);
        print("Restarting!");
        # Since we can never have enough data, we are now restarting the programm and running it again
    for i in range(0,100): # Spam the exit message to the queue to make sure all workers finish
        q.put('Exit');

        # Stop the workers when their work is done - Good job guys
    for i in range(0, len(cpuList)):
        cpuList[i].join()