OPENQASM 3.0;

// Declarations
qubit[4] q1;
qubit[4] q2;
bit[4] c1;
bit c2;

// Single qubit tests
h q1[0];
h q2;
h q1[{3,2}];
h q1[0:2];
h q1[-1];

// Two qubit tests
cx q1[0], q1[1];
cx q1, q2;
cx q2[0], q1;
cx q1, q2[0];
cx q1[{0,2}], q1[{1,3}];
cx q1[0:2], q1[1:3];
cx q1[-1], q2[-2];

// Test all gates
id q1[0];
x q1[1];
y q1[2];
z q1[3];
s q2[0];
sdg q2[1];

cz q1[2], q1[1];
swap q2[2], q2[1];

// Test measurements and resets
measure q1[3];
c2 = measure q1[0];
c1 = measure q1;
c1[1:3] = measure q1[1:3];
c1[{0,2}] = measure q2[{3,2}];
c1[-2] = measure q1[-3];

reset q1[0];
reset q1;
reset q1[1:3];
reset q2[{3,2}];