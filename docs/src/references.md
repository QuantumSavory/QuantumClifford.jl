# Suggested reading

For the basis of the tableaux methods first read [gottesman1998heisenberg](@cite) followed by the more efficient approach described in [aaronson2004improved](@cite).

The tableaux can be canonicalized (i.e. Gaussian elimination can be performed on them) in a number of different ways, and considering the different approaches provides useful insight. The following methods are implemented in this library:

- The default one: [garcia2012efficient](@cite)
- Useful when in need of tracing out a set of qubits: [audenaert2005entanglement](@cite)
- Useful when defining logical operators of codes: [gottesman1997stabilizer](@cite)

For the use of these methods in error correction and the subtle overlap between the two fields consider these resources. They are also useful in defining some of the specific constraints in commutation between rows in the tableaux:

- [steane2007tutorial](@cite)
- [calderbank1998quantum](@cite)
- [mackay2004sparse](@cite)
- [wilde2009logical](@cite)

These publications describe the uniform sampling of random stabilizer states:

- [koenig2014efficiently](@cite)
- [bravyi2020hadamard](@cite)
- [berg2020simple](@cite)
- [li2019measurement](@cite)

For circuit construction routines (for stabilizer measurements for a given code):
- [cleve1997efficient](@cite)
- [gottesman1997stabilizer](@cite) (and its erratum)
- [grassl2002algorithmic](@cite)
- [grassl2011variations](@cite)

Quantum Error Corecting Codes:
- 2D Color Code [Fowler_2011](@cite)
- 3D Color Code [bombin2015gauge](@cite)
- 2D Surface Code [bravyi1998quantum](@cite)
- 3D Surface Code[Vasmer_2019](@cite)
- Floquent Code [fahimniya2023faulttolerant](@cite) 
- Quantum Low-Density Parity Check (QLDPC) Code [Breuckmann_2021](@cite)
- Hypergraph Product Code [Krishna_2021](@cite)
- Bias-tailored Quantum LDPC Codes [Roffe_2023](@cite)

# References

```@bibliography
```
