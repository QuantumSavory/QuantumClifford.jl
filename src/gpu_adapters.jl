"""
copies the memory content of the object to CPU

You can only use this function if CUDA.jl is imported

For more advanced users `to_cpu(data, element_type)` will reinterpret elements of data and converts them to `element_type`.
For example based on your CPU architecture, if working with matrices of UInt32 is faster than UInt64, you can use `to_cpu(data, UInt32)`




```jldoctest
julia> using QuantumClifford: to_cpu, to_gpu

julia> using CUDA # without this import, to_cpu, to_gpu are just function

julia> stab = random_stabilizer(3)
- X_Z
+ _ZZ
+ __Z

julia> stab_gpu = to_gpu(stab);

julia> apply!(stab_gpu, sHadamard(1));

julia> stab_result_cpu = to_cpu(stab_gpu)
- Z_Z
+ _ZZ
+ __Z
```

```jldoctest
julia> using QuantumClifford: to_cpu, to_gpu

julia> using CUDA # without this import, to_cpu, to_gpu are just function

julia> pf_gpu = to_gpu(PauliFrame(1000, 2, 2));
julia> circuit = [sMZ(1, 1), sHadamard(2), sMZ(2, 2)];
julia> pftrajectories(pf_gpu, circuit);
julia> measurements = to_cpu(pf_gpu.measurements);
julia> sum(measurements, dims=1)
1×2 Matrix{Int64}:
 0  492
```

See also: [`to_gpu`](@ref)
"""
function to_cpu end


"""
copies the memory content of the object to GPU

You can only use this function if CUDA.jl is imported

For more advanced users `to_gpu(data, element_type)` will reinterpret elements of data and converts them to `element_type`.
For example based on your GPU architecture, if working with matrices of UInt64 is faster than UInt32, you can use `to_gpu(data, UInt64)`




```jldoctest
julia> using QuantumClifford: to_cpu, to_gpu

julia> using CUDA # without this import, to_cpu, to_gpu are just function

julia> stab = random_stabilizer(3)
- X_Z
+ _ZZ
+ __Z

julia> stab_gpu = to_gpu(stab);

julia> apply!(stab_gpu, sHadamard(1));

julia> stab_result_cpu = to_cpu(stab_gpu)
- Z_Z
+ _ZZ
+ __Z
```

```jldoctest
julia> using QuantumClifford: to_cpu, to_gpu

julia> using CUDA # without this import, to_cpu, to_gpu are just function

julia> pf_gpu = to_gpu(PauliFrame(1000, 2, 2));
julia> circuit = [sMZ(1, 1), sHadamard(2), sMZ(2, 2)];
julia> pftrajectories(pf_gpu, circuit);
julia> measurements = to_cpu(pf_gpu.measurements);
julia> sum(measurements, dims=1)
1×2 Matrix{Int64}:
 0  492
```

See also: [`to_cpu`](@ref)
"""
function to_gpu end
