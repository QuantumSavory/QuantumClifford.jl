OPENQASM 3.0;

// Two separate qubit registers to test offset arithmetic
qubit[2] a;
qubit[2] b;
bit[2]   ca;
bit[2]   cb;

h  a[0];
cx a[0], b[0];
x  b[1];

ca[0] = measure a[0];
ca[1] = measure a[1];
cb[0] = measure b[0];
cb[1] = measure b[1];
