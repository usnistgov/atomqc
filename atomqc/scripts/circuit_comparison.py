from jarvis.db.figshare import (
    get_wann_phonon,
    get_hk_tb,
    get_wann_electron,
)
from jarvis.io.qiskit.inputs import HermitianSolver
from collections import defaultdict
from matplotlib.gridspec import GridSpec
%matplotlib inline
import matplotlib.pyplot as plt
from qiskit.aqua.components.optimizers import COBYLA
the_grid = GridSpec(1, 3)
plt.rcParams.update({'font.size': 18})
plt.figure(figsize=(16,6))

info=defaultdict()
w, ef,atoms = get_wann_electron("JVASP-816")
kps = [[0.0, 0.0, 0.0], [0.5, 0., 0.5]]
for kk, k in enumerate(kps):
    print ()
    plt.subplot(the_grid[kk])

    hk = get_hk_tb(w=w, k=k)
    HS = HermitianSolver(hk)
    n_qubits = HS.n_qubits()
    Q = QuantumCircuitLibrary(n_qubits=n_qubits,reps=2)
    circs = [
        Q.circuit1(),
        Q.circuit2(),
        Q.circuit3(),
        Q.circuit4(),
        Q.circuit5(),
        Q.circuit6(),
    ]
    for ii, i in enumerate(circs):
        en, vqe_result, vqe = HS.run_vqe(var_form=i)

        vals, vecs = HS.run_numpy()
        diff = abs(en - vals[0])
        plt.axhline(vals[0]-ef,c='g')
        plt.plot(ii+1,en-ef,'s',c='b')
        print("VQE,numpy", ii, kk, en-ef, vals[0]-ef,diff)
        if kk==1:
          plt.ylim([-4,-1])
          plt.title('(b) Al X-point')
        else:
           plt.title('(a) Al $\Gamma$-point')
           plt.ylabel('Energy (eV)')
           plt.ylim([-13,-9])
        plt.xticks([1,2,3,4,5,6])
        plt.xlabel('Circuit number')
        #info[kk].append()
        #print(i)
    plt.tight_layout()
plt.subplot(the_grid[2])
wtbh, Ef, atoms = get_wann_electron("JVASP-35680") 
kpt = [0.5, 0., 0.5] # X-point
hk = get_hk_tb(w=wtbh, k=kpt)
HS = HermitianSolver(hk)
n_qubits = HS.n_qubits()
vals,vecs = HS.run_numpy()
plt.axhline(vals[0]-Ef,c='g')
#circ = QuantumCircuitLibrary(n_qubits=n_qubits).circuit6()
for i in range(1,6):
  print ('reps=',i)
  en, vqe_result, vqe = HS.run_vqe(reps=i)#var_form=circ)
  plt.plot(i,en-Ef,'s',c='b')
plt.xticks([1,2,3,4,5])
plt.title('(c) PbS X-point')
plt.xlabel('#Reps (Circuit-6 )')
plt.ylim([-10,-7])
plt.tight_layout()
