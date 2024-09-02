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

For quantum code construction routines:
- [cleve1997efficient](@cite)
- [gottesman1996class](@cite)
- [gottesman1997stabilizer](@cite)
- [yu2013all](@cite)
- [chao2018quantum](@cite)
- [kitaev2003fault](@cite)
- [fowler2012surface](@cite)
- [knill1996concatenated](@cite)
- [bravyi2024high](@cite)

For classical code construction routines:
- [muller1954application](@cite)
- [reed1954class](@cite)
- [raaphorst2003reed](@cite)
- [abbe2020reed](@cite)
- [djordjevic2021quantum](@cite)
- [hocquenghem1959codes](@cite)
- [bose1960class](@cite)
- [bose1960further](@cite)
- [error2024lin](@cite)

# References

```@bibliography
```