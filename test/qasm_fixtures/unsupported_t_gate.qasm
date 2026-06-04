OPENQASM 3.0;

qubit[1] q;
bit[1]   c;

// T gate is non-Clifford and should raise an informative error
t q[0];

c[0] = measure q[0];
