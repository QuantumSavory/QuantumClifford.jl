OPENQASM 3.0;

// Test all supported single-qubit Clifford gates
qubit[6] q;
bit[6]   c;

id   q[0];
x    q[0];
y    q[1];
z    q[2];
h    q[3];
s    q[4];
sdg  q[5];

// Two-qubit Clifford gates
cx   q[0], q[1];
cz   q[2], q[3];
swap q[4], q[5];

// Measurement
c[0] = measure q[0];
c[1] = measure q[1];
c[2] = measure q[2];
c[3] = measure q[3];
c[4] = measure q[4];
c[5] = measure q[5];
