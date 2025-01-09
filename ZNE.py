from qiskit import QuantumCircuit,transpile

def bruit(qc,n):
    qct=qc.copy()
    qct=transpile(qct,optimization_level=0,basis_gates=['cx','x','sx','rz','id'],approximation_degree=1)
    qc2=QuantumCircuit(qc.num_qubits,qc.num_clbits)
    for gate in qct:
        if gate[0].name=='cx':
            for k in range(n):
                qc2.cx(gate[1][0].index,gate[1][1].index)
        if gate[0].name=='x':
            qc2.x(gate[1][0].index)
        if gate[0].name=='sx':
            qc2.sx(gate[1][0].index)
        if gate[0].name=='rz':
            theta=gate[0].params[0]
            qc2.rz(theta,gate[1][0].index)
        if gate[0].name=='measure':
            i=gate[1][0].index
            j=gate[2][0].index 
            qc2.measure(i,j)
        if gate[0].name=='barrier':
            qc2.barrier(list(gate[1][k].index for k in range(len(gate[1]))))
    qc2=transpile(qc2,optimization_level=0,basis_gates=['cx','x','sx','rz','id'],approximation_degree=1)
    return qc2