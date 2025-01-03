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
- [steane1999quantum](@cite)
- [campbell2012magic](@cite)
- [anderson2014fault](@cite)
- [wang2024coprime](@cite)
- [voss2024multivariatebicyclecodes](@cite)
- [lin2024quantum](@cite)
- [bravyi2024high](@cite)
- [haah2011local](@cite)

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
- [golay1949notes](@cite)
- [huffman2010fundamentals](@cite)
- [bhatia2018mceliece](@cite)

For brute-force simulators:
- [fatima2021faster](@cite)
- [de2019massively](@cite)
- [markov2008simulating](@cite)
- [de2007massively](@cite)
- [markov2018quantum](@cite)

For efficient classical simulators:
- [gottesman1998heisenberg](@cite)
- [aaronson2004improved](@cite)
- [terhal2002classical](@cite)
- [bartlett2002efficient](@cite)
- [jozsa2008matchgate](@cite)

For Born rule probability estimators:
- [seddon2021quantifying](@cite)
- [pashayan2020estimation](@cite)
- [pashayan2019classical](@cite)
- [pashayan2015estimating](@cite)
- [rall2019simulation](@cite)
- [howard2017application](@cite)
- [veitch2012negative](@cite)
- [mari2012positive](@cite)

For pure-state sampling simulators:
- [bravyi2016improved](@cite)
- [garcia2017geometry](@cite)
- [bravyi2016trading](@cite)
- [bravyi2019simulation](@cite)
- [garcia2012efficient](@cite)
- [pashayan2022fast](@cite)

# References

```@bibliography
```