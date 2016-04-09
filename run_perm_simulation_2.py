from PERMSimulation import PERMSimulation
import numpy as np
import sys
import time

sim = PERMSimulation()

# Ts = np.linspace(0.2, 6.0, 30)
Ts = np.linspace(0.02, 3.0, 150)

sim_id = str(sys.argv[1])

if (sys.argv[2] == "fjc"):
    do_saw = False
elif(sys.argv[2] == "saw"):
    do_saw = True
else:
    sys.exit("unsupported simulation mode")

num_pol = int(sys.argv[3])

Tblock_start = int(sys.argv[4])
Tblock_end = int(sys.argv[5])
Ts = Ts[Tblock_start:Tblock_end]

print("Fs for this simulation batch: " + str(Fs))
print("Ts for this simulation batch: " + str(Ts))
print("Performing a %s simulation" % ("saw" if do_saw else "fjc"))

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