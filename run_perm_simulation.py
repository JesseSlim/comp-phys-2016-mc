from PERMSimulation import PERMSimulation
import numpy as np
import sys
import time

sim = PERMSimulation()

Ts = np.linspace(0.2, 6.0, 30)
Fs = np.concatenate((np.linspace(0.0, 0.1, 11), np.linspace(0.2, 1.0, 9)))

sim_id = str(sys.argv[1])

if (sys.argv[2] == "fjc"):
    do_saw = False
elif(sys.argv[2] == "saw"):
    do_saw = True
else:
    sys.exit("unsupported simulation mode")

num_pol = int(sys.argv[3])

Fblock_start = int(sys.argv[4])
Fblock_end = int(sys.argv[5])
Fs = Fs[Fblock_start:Fblock_end]

print("Fs for this simulation batch: " + str(Fs))

for F in Fs:
    for T in Ts:
        print("Starting simulation with parameters (sim id " + sim_id + "): ")
        print("T:   " + str(T))
        print("F:   " + str(F))
        filename = "results/perm-sim-id-%s-numpol-%d-T-%.2f-F-%.2f.npz" % (sim_id, num_pol, T, F)
        sim.T = T
        sim.F = F
        sim.enable_LJ_interaction = do_saw
        sim.enable_pruning_enriching = do_saw
        sim.initialise(num_pol, num_pol)
        sim.run_simulation()
        sim.save_results(filename)