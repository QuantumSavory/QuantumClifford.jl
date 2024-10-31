# [Generalized Stabilizer Representation](@id Generalized-Stabilizer-Overview)

Gottesman's introduction of stabilizer formalism in 1997 greatly impacted quantum complexity and coding
theory. The key insight of the Gottesman-Knill theorem lies in utilizing the Heisenberg representation[^1] for
quantum states, allowing classical simulations to work with only `n` Pauli operators, rather than processing
an exponentially large complex vector with approximately `2ⁿ` entries for an `n`-qubit state. However, this
approach is limited to stabilizer circuits with Clifford gates and measurements. While effective, the theorem
has a narrow scope, making it essential to generalize it for broader quantum circuit simulations. Theodore
Yoder[^2] introduced the generalized stabilizer representation to address this challenge.

# Advances in Stabilizer Formalism

Since its inception, the stabilizer formalism has undergone several improvements. Notable enhancements include:

```@raw html
<div class="mermaid">
timeline
    title Related Work in Generalization of the Gottesman-Knill Theorem
    1997 : Gottesman introduces stabilizer formalism and the Gottesman-Knill theorem.
    2002 : Bartlett et al. expand to continuous variable quantum computation.
    2004 : Aaronson and Gottesman improve measurement time complexity to 𝒪(n²).
    2006 : Anders and Briegel achieve 𝒪(n log n) speedup in time complexity with graph states.
    2012 : Bermejo-Vega and Van den Nest generalize to any finite Abelian group from n-qubits ℤ₂ⁿ.
    2012 : Yoder develops the Generalized Stabilizer with a novel state representation.
</div>
```

# Generalized Stabilizer Representation

The generalized stabilizer representation provides a flexible framework for simulating quantum circuits by:

- Enabling the representation of any quantum state, pure or mixed.
- Allowing simulations of arbitrary quantum circuits, including unitary operations, measurements, and
quantum channels.

This representation expands on the stabilizer formalism by incorporating non-stabilizer states and circuits,
enabling the simulation of non-Clifford gates and broader quantum channels for diverse quantum computations.

Unlike previous methods that may use a superposition of stabilizer states to represent arbitrary states,
this approach employs the tableau construction developed by Aaronson and Gottesman[^3]. This method implicitly
represents a set of orthogonal stabilizer states, forming a stabilizer basis capable of representing arbitrary
quantum states. Updating the tableau takes only twice as long as updating a single stabilizer, enabling
efficient updates of the entire stabilizer basis with minimal computational overhead.

# Strong Simulation vs Weak Simulation

Yoder[^2] investigated the classical simulation of quantum circuits in the context of `strong simulation`, which
focuses on calculating exact probabilities for specific measurement outcomes. A notable gap exists between weak
and strong simulation problems, with some circuits being `#P`-complete for strong simulation while weak simulation
is in `BPP`.The key differences between these two approaches are as follows:

```@raw html
<div class="mermaid">
mindmap
  root((Classical Simulation of Quantum Circuits))
    Strong Simulation
      Goal: Exact probability of specific outcome
      Example: Calculate probability of outcome 01
      Precision: High, requires precise probabilities
      Complexity: Can be #P-complete
    Weak Simulation
      Goal: Sample outcome close to quantum distribution
      Example: Sample outcomes like 00, 01, 10, 11
      Precision: Lower, approximate sampling
      Complexity: Generally in BPP
</div>
```

# Simulation of Quantum Channels

The generalized stabilizer representation enables the simulation of arbitrary quantum channels, beyond just
unitary gates and measurements. It does this by decomposing the Kraus operators of a channel into Pauli
operators from the state’s tableau, allowing for a broader range of quantum operations.

# Advantages of the Generalized Stabilizer

The proposed representation combines the rapid update capabilities of stabilizer states with the generality of
density matrices. Key features include:

- High update efficiency for unitary gates, measurements, and quantum channels, influenced by the sparsity of
the density matrix, `Λ(χ)`, which indicates the count of non-zero elements in `χ`.

-  Simulations maintain linear complexity with respect to the number of measurements, and the representation
remains straightforward, reflecting the principle that measurements simplify quantum states through collapse.

# Implications for Classical and Quantum Computation

Investigating stabilizer circuits enhances our understanding of classical and quantum computation. Simulating these
circuits is a complete problem in the classical complexity class `⊕L`, a subset of `P`, indicating that stabilizer
circuits may not be universal in classical computation contexts. Surprisingly, adding just one non-Clifford gate to
circuits with Clifford gates and measurements generally enables universal quantum computation—a contrast that highlights
intriguing questions about the computational boundaries between classical and quantum systems.

[^1]: [gottesman1998heisenberg](@cite)

[^2]: [yoder2012generalization](@cite)

[^3]: [gottesman1997stabilizer](@cite)
