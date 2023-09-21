import numpy as np
from scipy import linalg
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import random_statevector, Operator
from qiskit.extensions import UnitaryGate
import quimb as qu
import quimb.tensor as qtn

def pauli_y(qubit,system_size):
    """
    Parameters: 
    integer: the qubit on which is applied the Pauli Matrix Y (from 0 to system_size-1)
    integer: number of qubits of the whole system 
    
    Return: numpy array: the corresponding matrix operator
    """
    if qubit==0:
        S=np.array([[0,-1j],[1j,0]])
    else:
        S=np.array([[1,0],[0,1]])
    for j in range(1,system_size):
        if j==qubit:
            S=np.kron(np.array([[0,-1j],[1j,0]]),S)
        else:
            S=np.kron(np.array([[1,0],[0,1]]),S)
    return S
    

def pauli_z(qubit,system_size):
    """
    Parameters: 
    integer: the qubit on which is applied the Pauli Matrix Z (from 0 to system_size-1)
    integer: number of qubits of the whole system 
    
    Return: numpy array: the corresponding matrix operator
    """
    if qubit==0:
        S=np.array([[1,0],[0,-1]])
    else:
        S=np.array([[1,0],[0,1]])
    for j in range(1,system_size):
        if j==qubit:
            S=np.kron(np.array([[1,0],[0,-1]]),S)
        else:
            S=np.kron(np.array([[1,0],[0,1]]),S)
    return S

def pauli_p(qubit,system_size):
    """
    Parameters: 
    integer: the qubit on which is applied the Pauli Matrix X (from 0 to system_size-1)
    integer: number of qubits of the whole system 
    
    Return: numpy array: the corresponding matrix operator
    """
    if qubit==0:
        S=np.array([[0,1],[0,0]])
    else:
        S=np.array([[1,0],[0,1]])
    for j in range(1,system_size):
        if j==qubit:
            S=np.kron(np.array([[0,1],[0,0]]),S)
        else:
            S=np.kron(np.array([[1,0],[0,1]]),S)
    return S

def pauli_m(qubit,system_size):
    """
    Parameters: 
    integer: the qubit on which is applied the Pauli Matrix X (from 0 to system_size-1)
    integer: number of qubits of the whole system 
    
    Return: numpy array: the corresponding matrix operator
    """
    if qubit==0:
        S=np.array([[0,0],[1,0]])
    else:
        S=np.array([[1,0],[0,1]])
    for j in range(1,system_size):
        if j==qubit:
            S=np.kron(np.array([[0,0],[1,0]]),S)
        else:
            S=np.kron(np.array([[1,0],[0,1]]),S)
    return S



def random_gaussian_interactions(J,omegas):
    """
    Parameters: 
     numpy array: hermitian matrix of coupling strength
     list of float: qubits frequencies
    
    Return: numpy array: the hamiltonian operator defined using raising and lowering operators
    """
    nb_qubits=len(J)
    hamiltonian=np.zeros([2**nb_qubits,2**nb_qubits],dtype=np.complex128)
    for i in range(nb_qubits):
        for j in range(i+1,nb_qubits):
            z_chain=[pauli_z(k,nb_qubits) for k in range(i+1,j)]
            hamiltonian-=J[i,j]*np.linalg.multi_dot([pauli_p(i,nb_qubits)]+z_chain+[pauli_m(j,nb_qubits)])
            hamiltonian-=J[j,i]*np.linalg.multi_dot([pauli_p(j,nb_qubits)]+z_chain+[pauli_m(i,nb_qubits)])
        hamiltonian+=omegas[i]/2*pauli_y(i,nb_qubits)
    return hamiltonian


def circuit_from_tensor_circuit(n_qubits, tensor_circ):
    qc_final = QuantumCircuit(n_qubits)
    for gate in tensor_circ.gates:
        gate_qubits = [abs(n_qubits-1-a) for a in gate.qubits]
        # print(gate_qubits)
        if gate.label == 'X':
            qc_final.x(*gate_qubits)
        elif gate.label == 'Y':
            qc_final.y(*gate_qubits)
        elif gate.label == 'Z':
            qc_final.z(*gate_qubits)
        elif gate.label == 'U3':
            qc_final.u(*gate.params, *gate_qubits)
        elif gate.label == 'CX':
            # gate_qubits = [(2-a) for a in gate.qubits]
            qc_final.cx(*gate_qubits)
        else:
            print('Abnormal gate, be careful')
    return qc_final 

def tensor_fidelity(unitary1, unitary2):
    fid_tab =[]
    n_qubits = int(np.log(len(unitary1[0,:]))/np.log(2))
    for _ in range(5000):
        psi0 = qu.rand_ket(2**n_qubits)
        psif_exact = unitary1 @ psi0
        psif_apprx = unitary2 @ psi0
        fid_tab.append(100 * qu.fidelity(psif_apprx, psif_exact))
    
    return np.mean(fid_tab)

def compute_fidelity2(qc1,qc2):
    n_qubits=qc1.num_qubits
    unitary_op1 = Operator(qc1).data
    unitary_op2 = Operator(qc2).data
    final_fid = [] 
    for _ in range(5000):
        state = random_statevector(2**n_qubits, seed=None).data
        vec1 = unitary_op1 @ state
        vec2 = unitary_op2 @ state
        final_fid.append(abs(np.transpose(np.conjugate(vec1)) @ vec2))
    # print(final_fid)
    final_fid=np.mean(final_fid)
    return final_fid

def single_qubit_layer(circ, gate_round=None):
    """Apply a parametrizable layer of single qubit ``U3`` gates.
    """
    for i in range(circ.N):
        # initialize with random parameters
        params = qu.randn(3, dist='uniform')
        circ.apply_gate(
            'U3', *params, i, 
            gate_round=gate_round, parametrize=True)
        
def two_qubit_layer(circ, two_qubit_gate, reverse=False, gate_round=None): #WARNING qubits in opposite odering when going from quimb to qiskit and vice versa
    """Apply a layer of constant entangling gates.
    """
    regs = range(0, circ.N - 1)
    if reverse:
        regs = reversed(regs)
    
    for i in regs:
            circ.apply_gate(two_qubit_gate, i, i+1, gate_round=gate_round)
            
def ansatz_circuit(n, depth, two_qubit_gate, **kwargs): #WARNING qubits in opposite odering when going from quimb to qiskit and vice versa
    """Construct a circuit of single qubit and entangling layers.
    """
    circ = qtn.Circuit(n, **kwargs)
    
    for r in range(depth):
        # single qubit gate layer
        single_qubit_layer(circ, gate_round=r)
        
        # alternate between forward and backward two_qubit_gate layers
        two_qubit_layer(
            circ, two_qubit_gate=two_qubit_gate, gate_round=r, reverse= (r % 2 == 0) )
        
    # add a final single qubit layer
    single_qubit_layer(circ, gate_round=r + 1)
    
    return circ

def tensor_network_recompilation(target, depth, ansatz_fct = ansatz_circuit, two_qubit_gate = 'CX', n_iter = 5000, n_hop = 20, temperature = 10, print_bool = True):

    if type(target) == np.ndarray:
        
        unitary = target
        n_qubits = int(np.log(len(unitary[0,:]))/np.log(2))
        
        qc_original = QuantumCircuit(n_qubits)
        qc_original.append(UnitaryGate(unitary),list(range(n_qubits)))
        coupling_map=[[i,i+1] for i in range(n_qubits-1)]+[[i+1,i] for i in range(n_qubits-1)]
        qc_original_trans = transpile(qc_original, optimization_level = 0,  initial_layout = list(np.arange(n_qubits)), coupling_map=coupling_map, basis_gates=['id', 'rz', 'sx', 'x', 'cx'])
        if print_bool: print('Original number of CNOT : '+str(qc_original_trans.count_ops()['cx']))
        
    if type(target) == QuantumCircuit:
        
        unitary = Operator(target).data
        n_qubits =  int(np.log(len(unitary[0,:]))/np.log(2))
        
        qc_original = target
        coupling_map=[[i,i+1] for i in range(n_qubits-1)]+[[i+1,i] for i in range(n_qubits-1)]
        qc_original_trans = transpile(qc_original, optimization_level = 0,  initial_layout = list(np.arange(n_qubits)), coupling_map=coupling_map, basis_gates=['id', 'rz', 'sx', 'x', 'cx'])
        if print_bool: print('Original number of CNOT : '+str(qc_original_trans.count_ops()['cx']))
        
    #Target tensor from unitary
    U = qtn.Tensor(
        data = unitary.reshape([2] * (2 * n_qubits)),
        inds = [f'k{i}' for i in range(n_qubits)]+[f'b{i}' for i in range(n_qubits)],
        tags = {'U_TARGET'}
        )
    #Random occurence of the Ansatz circuit
    
    qc_rand = ansatz_fct(n_qubits, depth, two_qubit_gate = two_qubit_gate)
    if print_bool: print('Ansatz number of CNOT : ' + str(circuit_from_tensor_circuit(n_qubits, qc_rand).count_ops()['cx']))
    
    #Tensor representation of the Ansatz circuit
    
    V = qc_rand.get_uni()
    
    #Loss function and optimizer definition
    def loss(V, U):
        return 1 - abs((V.H & U).contract(optimize='auto-hq')) / 2**n_qubits
    
    tnopt = qtn.TNOptimizer(
        V,                        # the tensor network we want to optimize
        loss,                     # the function we want to minimize
        loss_constants={'U': U},  # supply U to the loss function as a constant TN
        tags=['U3'],              # only optimize U3 tensors
        autodiff_backend='jax',   # use 'autograd' for non-compiled optimization
        optimizer='L-BFGS-B',     # the optimization algorithm
    )
    
    #Optimization
    V_opt = tnopt.optimize_basinhopping(n = n_iter, nhop = n_hop, temperature = temperature)
    V_opt_dense = V_opt.to_dense()
    tensor_fid = tensor_fidelity(unitary, V_opt_dense)
    if print_bool: print(f"Fidelity (tensor): {tensor_fid:.10f} %")
    
    #Updating quantum circuit and translate it to qiskit
    qc_rand.update_params_from(V_opt)
    qc_final = circuit_from_tensor_circuit(n_qubits, qc_rand)
    if print_bool: print(f"Fidelity (circuit): {100*compute_fidelity2(qc_original, qc_final):.10f} %")
    
    return qc_final

#Add final measurement in the y direction
def add_measure(qc):
    nb_qubits=qc.num_qubits
    qc2=qc.copy()
    for i in range(nb_qubits):
        qc2.rx(np.pi/2,i)
    qc2.measure_all()
    qc_meas=transpile(qc2,optimization_level =0,basis_gates=['cx','x','sx','rz','id'],approximation_degree=1)
    return qc_meas