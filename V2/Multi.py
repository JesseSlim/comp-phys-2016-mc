import numpy as np
import time;
from multiprocessing import Process, Queue, cpu_count

# This class facilitates multiprocesser 
class Multi:
    
    ###
    # You can initialize the class with 5 variables:
    # Main is a reference to the main function which in adition to controlling the workers can do actions
    # Worker is the specific function which the workers shoudl complete
    # Var are variables which are passed to the Main and the workers
    # Add is a change in these variables. If you want them to stay the same put them to zero.
    # The Add variable should have the same size as the Var varable.
    # Stop is the moment 
    # UseAll can be set to false to keep 1 cpu free
    ###
    def __init__(self, Worker, Var, Add, Stop, UseAll=True):
        # First we check if the input is correct
        if (not hasattr(Worker, '__call__')):
            print("Error: Worker is not a function");
            return;
        if (np.size(Var) != np.size(Add) or np.size(Var) != np.size(Stop)):
            print("Error: The Var, Add and Stop variables don't have the same size");
            return;
        
        # Now everything is correct so we start up the function.
        # First of all we start the communication queue's
        qm = Queue(); # Queue which communicates from the master to the workers
        qw = Queue(); # Queue which communicates from the workers to the master
        qs = Queue(); # Queue which communicates settings
        # Start the workers
        cpuList = []; # List of the workers
        cpuCount = cpu_count() - (not UseAll);
        cpuCount = 1;
        cpuStatus = np.zeros([cpuCount]);
        for i in range(0,cpuCount): #Cpu count utilizes your computer. use cpu_count()-1 to do other things.
            cpuList.append(Process(target=Worker, args=(qm,qw)))
            cpuList[i].start();
            qm.put((Var, i));
            Var = Var + Add;
            cpuStatus[i] = 1;
        
        # Increment the used values and run till it is complete;
        while (not np.all(np.greater_equal(Var,Stop))):
            if (not qw.empty()):
                i = qw.get();
                qm.put((Var, i));
                Var = Var + Add
            time.sleep(0.01);
        
        # Tell everyone this is the last job:
        for i in range(0,cpuCount):
            qm.put("Done");
        
        # Wait for all jobs ot finish
        while (np.sum(cpuStatus) > 0):
            i = qw.get();
            cpuStatus[int(i)] = 0;
            time.sleep(0.5);
        
        # All jobs have been completed. Stop the workers.s
        for i in range(0, len(cpuList)):
            cpuList[i].join();
            
        # Done
        print("Done");