# [ECC example with Perturbative Expansions](@id noisycircuits_pf_ecc_example)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.Experimental.NoisyCircuits
end
``` 

!!! warning "Unstable"
    This is experimental functionality with an unstable API.

This example demonstrates the simulation of a noisy quantum circuit using perturbative expansions, where error correction is incorporated through a DecoderCorrectionGate. After encoding the Steane code and introducing noise, the decoder correction gate applies the most likely correction based on the measured syndrome, leading to a noticeable improvement in the logical success rate. This illustrates how perturbative expansion methods can be combined with error-correcting circuits to study performance under noise.

Some of the relevant structures and functions used in this code are:
 - Register[@ref]
 - TableDecoder[@ref]
 - DecoderCorrectionGate[@ref] 
 - petrajectories[@ref]

Import with
`using QuantumClifford.Experimental.NoisyCircuits` and `using QuantumClifford.ECC: DecoderCorrectionGate`

These imports enable the simulation of noisy circuits using the perturbative expansion method and error correction using the DecoderCorrectionGate

Below is an example of this functionality in use:


```@example 1
using QuantumClifford
using QuantumClifford.ECC: Steane7, parity_checks, naive_encoding_circuit, naive_syndrome_circuit, TableDecoder, DecoderCorrectionGate
code = parity_checks(Steane7()) #creates the stabilizer tableau for the Steane code
register = Register(one(MixedDestabilizer, 13), 6) # A register with 13 qubits, 7 data qubits and 6 ancillas
decoder = TableDecoder(code) # creates a decoder based on the Steane code
noiseless_state = one(MixedDestabilizer, 7)
encoding_circuit = naive_encoding_circuit(code) #encodes the circuit as a state within the Steane codespace
mctrajectory!(noiseless_state, encoding_circuit) # applies each gate in ecirc to noiseless_state
syndrome_circuit,_,_ = naive_syndrome_circuit(code) # caputures the operations needed to extract the error syndrome
verify = VerifyOp(noiseless_state, indices = 1:7) # verifies the first 7 qubits of the register with the prepared state
epsilon = 0.1
noise = NoiseOp(UnbiasedUncorrelatedNoise(epsilon), (1:7)) # introduces noise into the circuit

circuit = vcat(encoding_circuit, noise, syndrome_circuit, verify) 
```
we can run perturbative expansions with [`petrajectories`](@ref).
```@example 1
pet = petrajectories(register, circuit; max_order=3)
pet
```
We see that the true_success rate is pretty terrible.
In order to improve this rate with the help of error correction, we can use the DecoderCorrectionGate:

```@example 1
correction_gate = DecoderCorrectionGate(decoder, data_qubits = 1:7, syndrome_bits = 1:6) # guesses the most probable error correcting Pauli Operation based on the recorded syndrome
```
Let's include the correction gate in our list of circuit operations and place it after syndrome collection.
```@example 1
circuit = vcat(encoding_circuit, noise, syndrome_circuit, correction_gate, verify) 
pet = petrajectories(register, circuit; max_order=3)
pet
```
We can see that the true_success rate has gone up significantly.