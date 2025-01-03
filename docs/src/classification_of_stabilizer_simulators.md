# [Classification of Stabilizer Simulators](@id classification_of_stabilizer_simulators)

Brute-force simulation algorithms, such as *Schrödinger-style* ([fatima2021faster](@cite)), *Feynman-style*
([de2019massively](@cite), [markov2008simulating](@cite), [de2007massively](@cite)), and *hybrid* simulators
([markov2018quantum](@cite)) provide high precision for simulating *universal* quantum circuits, but
they can become highly resource-intensive for circuits with moderate width (around *40* qubits) or depth.
Alternatively, efficiently classically simulable quantum circuits, like stabilizer circuits, can be
simulated using the *Gottesman-Knill* theorem, allowing the simulation of thousands of qubits with hundreds
of thousands of gates. Research to overcome the limitations of these methods falls into two main categories:
*Born rule probability estimators* that use a quasi-probabilistic representation of the density matrix, and
*pure-state sampling simulators*.

### Aaronson and Gottesman's Simulator

- Introduced in *Improved Simulation of Stabilizer Circuits* [aaronson2004improved](@cite)
- **Efficient for Stabilizer Circuits**: Introduced a classical simulation algorithm efficient
for stabilizer circuits.
- **Simulates Non-Stabilizer Circuits**: Handles non-stabilizer circuits with **exponential**
run-time cost depending on the number of non-stabilizer gates.
- **Limitation**: The run-time does not depend on the specific properties of the additional
non-stabilizer gates, incurring a heavy penalty even for small deviations from stabilizer gates.

### Research Categories to Overcome Limitation

#### **Born Rule Probability Estimators**

- **Quasi-Probabilistic Representation**
  - *Quantifying quantum speedups: improved classical simulation from tighter magic monotones* ([seddon2021quantifying](@cite))
  - *From estimation of quantum probabilities to simulation of quantum circuits* ([pashayan2020estimation](@cite))
  - *On the classical simulability of quantum circuits* ([pashayan2019classical](@cite))
  - *Estimating outcome probabilities of quantum circuits using quasiprobabilities* ([pashayan2015estimating](@cite))
  - *Simulation of Qubit Quantum Circuits via Pauli Propagation* ([rall2019simulation](@cite))
  - *Application of a resource theory for magic states to fault-tolerant quantum computing* ([howard2017application](@cite))
  - *Negative Quasi-Probability as a Resource for Quantum Computation* ([veitch2012negative](@cite))
  - *Positive Wigner functions render classical simulation of quantum computation efficient* ([mari2012positive](@cite))

#### **Pure-State Sampling Simulators**

- **Bravyi and Gosset Algorithms**
    - *Improved classical simulation of quantum circuits dominated by Clifford gates* ([bravyi2016improved](@cite))
    - *Efficient Inner-product Algorithm for Stabilizer States* ([garcia2012efficient](@cite))
    - *On the geometry of stabilizer states* ([garcia2017geometry](@cite))
    - *Trading classical and quantum computational resource* ([bravyi2016trading](@cite))
    - *Simulation of quantum circuits by low-rank stabilizer decompositions* ([bravyi2019simulation](@cite))
    - *Fast Estimation of Outcome Probabilities for Quantum Circuits* ([pashayan2022fast](@cite))

### Quasi-Probabilistic Simulators

- **Purpose**: Produce additive precision estimates of Born rule probabilities.
- **Representation**: Density matrices are expressed as a linear combination of a preferred set of operators (frame).

#### **Frame Choices**:

- **Examples**:
  - **Weyl-Heisenberg displacement operators**
     - *From estimation of quantum probabilities to simulation of quantum circuits* ([pashayan2020estimation](@cite))
     - *On the classical simulability of quantum circuits* ([pashayan2019classical](@cite))
     - *Simulation of Qubit Quantum Circuits via Pauli Propagation* ([rall2019simulation](@cite))
  - **Stabilizer states**
     - *Quantifying quantum speedups: improved classical simulation from tighter magic monotones* ([seddon2021quantifying](@cite))
     - *Application of a resource theory for magic states to fault-tolerant quantum computing* ([howard2017application](@cite))
  - **Phase-point operators**
      - *Estimating outcome probabilities of quantum circuits using quasiprobabilities* ([pashayan2015estimating](@cite))
      - *Negative Quasi-Probability as a Resource for Quantum Computation* ([veitch2012negative](@cite)
      - *Positive Wigner functions render classical simulation of quantum computation efficient* ([mari2012positive](@cite))

#### **Special Mention**:

- **Dyadic Frame Simulator**:
  - Introduced in *Quantifying quantum speedups: improved classical simulation from tighter magic monotones* ([seddon2021quantifying](@cite))
  - Decomposes density matrices into stabilizer *dyads* ``|L\rangle\langle R|``.
  - Circuits promoted to universality using magic states.

#### **Run-Time**:

- Depends quadratically on the dyadic negativity.
- Dyadic negativity measures deviation from convex combinations of stabilizer dyads.

### Bravyi and Gosset Algorithms

- **BG-Estimation Algorithm**: Produces multiplicative precision estimates of Born rule probabilities.
- **BG-Sampling Algorithm**: Samples approximately from the quantum circuit outcome distribution.

#### **Methodology**:

- **State Representation**: Initial states are expressed as a linear combination of stabilizer states.
- **Efficiently Simulable Circuits**:
  - Superposition of polynomially many stabilizer states.
  - Clifford gates and computational basis measurements.
- **Promoting Universality**: Allowing magic states in initial conditions.

#### **Run-Time Dependence**:

- Linear in the stabilizer rank of the state.
- **Stabilizer Rank**: Minimal number of stabilizer states required to represent the state as a linear combination.
- **Approximation and Stabilizer Extent**: Approximate stabilizer rank can be bounded by the *stabilizer extent*
``\xi`` divided by ``\epsilon^2``, where ``\epsilon`` is the approximation error.

### Sum Over Cliffords Algorithm

 - Introduced in *Simulation of quantum circuits by low-rank stabilizer decompositions* ([bravyi2019simulation](@cite))
- **Variant**: Simulates non-Clifford gates using a linear combination of Clifford gates.
- **Run-Time**: Scales linearly with the stabilizer extent of states like ``|T_\phi^\dagger \rangle``.

### Mixed-State Stabilizer Rank Simulator

  - Introduced in *Quantifying quantum speedups: improved classical simulation from tighter magic monotones* ([seddon2021quantifying](@cite))
- **Improvements**:
  - Generalized BG-sampling algorithm to include mixed states.
  - Improved run-time dependence on error tolerance for approximate sampling.

#### **Mixed-State Extent**:

- **Definition**: Quantity governing run-time for mixed-state simulators.
- **Comparison**: For $n$-qubit product states, dyadic negativity, stabilizer extent, and mixed-state extent are equivalent.
