%matplotlib inline
import matplotlib.pyplot as plt
%matplotlib inline
import matplotlib.pyplot as plt
from jarvis.db.figshare import get_wann_electron,get_wann_phonon
from jarvis.io.qiskit.inputs import get_bandstruct
from jarvis.core.kpoints import Kpoints3D as Kpoints
from jarvis.db.figshare import get_wann_electron, get_wann_phonon, get_hk_tb
from jarvis.io.qiskit.inputs import HermitianSolver
import numpy as np
from jarvis.db.figshare import get_hk_tb
import numpy as np
from qiskit import BasicAer
from qiskit.aqua.operators import WeightedPauliOperator, MatrixOperator, op_converter
from qiskit import Aer
from qiskit.aqua import QuantumInstance, aqua_globals
from qiskit.aqua.algorithms import VQE, NumPyMinimumEigensolver
from jarvis.db.figshare import get_wann_electron,get_wann_phonon
from qiskit.aqua.components.optimizers import COBYLA, L_BFGS_B, SLSQP,SPSA,POWELL,CG, ADAM
from jarvis.core.kpoints import Kpoints3D as Kpoints
from matplotlib.gridspec import GridSpec

the_grid = GridSpec(1, 3)
plt.rcParams['figure.figsize'] = (24, 8)
plt.rcParams.update({'font.size': 18})
w,ef,atoms=get_wann_electron(jid="JVASP-816")
#hk=get_hk_tb(w=w,k=[0.,0.,0.])
hk=get_hk_tb(w=w,k=[0.5,0.5,0.])
Hamil_Mat=MatrixOperator(hk)
# backend= Aer.get_backend("statevector_simulator")
Hamil_Qop = op_converter.to_weighted_pauli_operator(Hamil_Mat)
H2_op=Hamil_Qop
#H2_op=HermitianSolver(hk)
optimizers = [COBYLA(maxiter=500), L_BFGS_B(maxiter=500), SLSQP(maxiter=500),CG(maxiter=500),SPSA(maxiter=500)]
converge_cnts = np.empty([len(optimizers)], dtype=object)
converge_vals = np.empty([len(optimizers)], dtype=object)

for i, optimizer in enumerate(optimizers):
    print('\rOptimizer: {}        '.format(type(optimizer).__name__), end='')
    aqua_globals.random_seed = 50
    var_form = EfficientSU2(4)

    counts = []
    values = []
    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)

    vqe = VQE(H2_op, var_form, optimizer, callback=store_intermediate_result,
              quantum_instance=QuantumInstance(backend=BasicAer.get_backend('statevector_simulator')))
    result = vqe.compute_minimum_eigenvalue(operator=H2_op)
    converge_cnts[i] = np.asarray(counts)
    converge_vals[i] = np.asarray(values)
print('\rOptimization complete      ');
plt.subplot(the_grid[0])

plt.rcParams['figure.figsize'] = (24, 8)
plt.rcParams.update({'font.size': 18})
for i, optimizer in enumerate(optimizers):
    plt.plot(converge_cnts[i], converge_vals[i], label=type(optimizer).__name__)
plt.xlabel('Eval count')
plt.ylabel('Energy (eV)')
#pylab.title('Energy convergence for various optimizers')
plt.legend(loc='upper right');
plt.xlim([0,400])
plt.title('(a)')


reps=5
plt.subplot(the_grid[2])
plt.title('(c)')
#we,ef,atoms=get_wann_electron(jid="JVASP-816")
wp,atoms=get_wann_phonon(jid="JVASP-816",factor=30)
info = {}
line_density=5
kpoints = Kpoints().kpath(atoms, line_density=line_density)
labels = kpoints.to_dict()["labels"]
kpts = kpoints.to_dict()["kpoints"]
new_kp = []
new_labels = []
count = 0
kp = np.arange(len(kpts))
for i, j in zip(kp, labels):
    if j != "":
        if count > 1 and count < len(labels) - 1:
            if labels[count] != labels[count + 1]:
                new_kp.append(i)
                new_labels.append("$" + str(j) + "$")
        else:
            new_kp.append(i)
            new_labels.append("$" + str(j) + "$")
    count += 1
print("kpts", len(kpts))

eigvals_q = []
eigvals_np = []
for ii, i in enumerate(kpts):
  #if ii<5:
    try:
        print("kp=", ii, i)
        hk = get_hk_tb(w=wp, k=i)
        HS = HermitianSolver(hk)
        vqe_vals, _ = HS.run_vqd(reps=reps)
        np_vals, _ = HS.run_numpy()
        print("np_vals", np_vals)
        print("vqe_vals", vqe_vals)
        eigvals_q.append(vqe_vals)
        eigvals_np.append(np_vals)
        # break
    except Exception as exp:
        print(exp)
        pass
eigvals_q =np.array(eigvals_q)
eigvals_np =np.array(eigvals_np)


for ii, i in enumerate(eigvals_q.T):
    if ii == 0:
        plt.plot(i, '*',c="b", label="VQD")
    else:
        plt.plot(i, '*',c="b")

 kpoints = Kpoints().kpath(atoms,line_density=line_density)
labels = kpoints.to_dict()["labels"]
kpts = kpoints.to_dict()["kpoints"]
eigs = wp.band_structure_eigs(kpath=kpts).T
for i, ii in enumerate(eigvals_np.T):
     if i == 0:
            plt.plot(ii, color="g", label="Exact")
     else:
            plt.plot(ii, color="g")
tol=6
plt.ylim([tol, np.max(eigvals_q)])
plt.xticks(new_kp, new_labels)
plt.ylabel('Freq. (cm$^{-1}$)')
plt.legend()
#"""





plt.subplot(the_grid[1])
plt.title('(b)')
wp,ef,atoms=get_wann_electron(jid="JVASP-816")
#wp,atoms=get_wann_phonon(jid="JVASP-816")
info = {}
line_density=5
kpoints = Kpoints().kpath(atoms, line_density=line_density)
labels = kpoints.to_dict()["labels"]
kpts = kpoints.to_dict()["kpoints"]
new_kp = []
new_labels = []
count = 0
kp = np.arange(len(kpts))
for i, j in zip(kp, labels):
    if j != "":
        if count > 1 and count < len(labels) - 1:
            if labels[count] != labels[count + 1]:
                new_kp.append(i)
                new_labels.append("$" + str(j) + "$")
        else:
            new_kp.append(i)
            new_labels.append("$" + str(j) + "$")
    count += 1
print("kpts", len(kpts))


eigvals_q = []
eigvals_np = []
for ii, i in enumerate(kpts):
  #if ii<5:
    try:
        print("kp=", ii, i)
        hk = get_hk_tb(w=wp, k=i)
        HS = HermitianSolver(hk)
        vqe_vals, _ = HS.run_vqd(reps=reps)
        np_vals, _ = HS.run_numpy()
        print("np_vals", np_vals)
        print("vqe_vals", vqe_vals)
        eigvals_q.append(vqe_vals)
        eigvals_np.append(np_vals)
        # break
    except Exception as exp:
        print(exp)
        pass
eigvals_q =np.array(eigvals_q)-ef
eigvals_np =np.array(eigvals_np)-ef


for ii, i in enumerate(eigvals_q.T):
    if ii == 0:
        plt.plot(i-ef, '*',c="b", label="VQD")
    else:
        plt.plot(i-ef, '*',c="b")

kpoints = Kpoints().kpath(atoms,line_density=line_density)
labels = kpoints.to_dict()["labels"]
kpts = kpoints.to_dict()["kpoints"]
eigs = wp.band_structure_eigs(kpath=kpts).T
for i, ii in enumerate(eigvals_np.T):
#for i, ii in enumerate(eigvals_np):
     if i == 0:
            plt.plot(ii-ef, color="g", label="Exact")
     else:
            plt.plot(ii-ef, color="g")
tol=0
#plt.ylim([tol, np.max(eigvals_q)])
plt.xticks(new_kp, new_labels)
plt.ylabel('Freq. (cm$^{-1}$)')
plt.legend()
