using Test
using CUDA
using QuantumClifford
using Random
include("../src/gpu_canon.jl")    
using .gpu_canon


const MIN_BLOCKS = 36


function test_copy_kernel()
    println("Running test_copy_kernel...")
    arr = rand(Float32, 100)
    d_arr = custom_gpu_linear_copy(arr)
    copied_back = Array(d_arr)
    @test arr ≈ copied_back 
    println("Passed test_copy_kernel.")
end

function test_swap_columns()
    println("Running test_swap_columns...")
    mat = CUDA.fill(UInt64(0), 4, 4)
    mat[:, 1] .= UInt64(2)
    mat[:, 2] .= UInt64(3)
    gpu_swap_columns!(mat, 1, 2)
    result = Array(mat)
    @test all(result[:, 1] .== UInt64(3)) && all(result[:, 2] .== UInt64(2))
    println("Passed test_swap_columns.")
end

function test_swap_columns_and_phases()
    println("Running test_swap_columns_and_phases...")
    mat = CUDA.fill(UInt64(0), 4, 4)
    phases = CUDA.fill(UInt8(0), 4)
    CUDA.@allowscalar begin
        mat[:, 1] .= UInt64(2)
        mat[:, 2] .= UInt64(3)
        phases[1] = UInt8(1)
        phases[2] = UInt8(2)
    end
    
    gpu_swap_columns_and_phases!(mat, phases, 1, 2)
    result_mat = Array(mat)
    result_phases = Array(phases)
    @test all(result_mat[:, 1] .== UInt64(3)) && all(result_mat[:, 2] .== UInt64(2))
    @test result_phases[1] == UInt8(2) && result_phases[2] == UInt8(1)
    println("Passed test_swap_columns_and_phases.")
end

function test_pivot_search()
    println("Running test_pivot_search...")
    mat = CUDA.fill(UInt64(0), 4, 4)
    
    CUDA.@allowscalar mat[1, 2] = UInt64(1)
    
    pivot = gpu_pivot_search!(mat, 1, UInt64(1), 1, 4)
    @test pivot == 2
    println("Passed test_pivot_search.")
end

function test_update_kernel()
    println("Running test_update_kernel...")
    mat = CUDA.fill(UInt64(0), 4, 4)
    phases = CUDA.fill(UInt8(0), 4)
    
    CUDA.@allowscalar begin
        mat[1, 1] = UInt64(1)
        phases[1] = UInt8(1)
    end

    gpu_update_kernel!(mat, phases, 1, 1, UInt64(1), 4, 4, 1, true)
    
    CUDA.synchronize()

    result_mat = Array(mat)
    result_phases = Array(phases)
    
    @test result_mat[1, 1] == UInt64(1) 
    println("Passed test_update_kernel.")
end


function test_canonicalization()
    println("Running test_canonicalization...")
    cpu_stab = random_stabilizer(4)

    CUDA.@allowscalar begin
        gpu_stab = cpu_to_gpu_stabilizer(cpu_stab)
        gpu_canon.canonicalize!(gpu_stab; phases=true)
        cpu_canon = gpu_to_cpu_stabilizer(gpu_stab)
    end

    @test cpu_canon != cpu_stab 
    println("Passed test_canonicalization.")
end


function test_conversion_functions()
    println("Running test_conversion_functions...")
    cpu_stab = random_stabilizer(4)
    gpu_stab = cpu_to_gpu_stabilizer(cpu_stab)
    cpu_back = gpu_to_cpu_stabilizer(gpu_stab)
    @test cpu_stab == cpu_back
    println("Passed test_conversion_functions.")
end

function test_cpu_gpu_stabilizer_equivalence()
    println("Running test_cpu_gpu_stabilizer_equivalence...")
    cpu_stab = random_stabilizer(4)
    gpu_stab = cpu_to_gpu_stabilizer(cpu_stab)
    cpu_stab_converted = gpu_to_cpu_stabilizer(gpu_stab)
    @test cpu_stab == cpu_stab_converted
    println("Passed test_cpu_gpu_stabilizer_equivalence.")
end

function test_canonicalize!()
    println("Running test_canonicalize!...")

    CUDA.@allowscalar begin
        cpu_stab = random_stabilizer(4)
        gpu_stab = cpu_to_gpu_stabilizer(cpu_stab)
        gpu_canon.canonicalize!(gpu_stab; phases=true)
        cpu_canon = gpu_to_cpu_stabilizer(gpu_stab)
    end

    @test cpu_canon != cpu_stab  
    println("Passed test_canonicalize!.")
end


function test_random_stabilizer()
    println("Running test_random_stabilizer...")

    CUDA.@allowscalar begin
        cpu_stab = random_stabilizer(4)
        gpu_stab = cpu_to_gpu_stabilizer(cpu_stab)
        gpu_canon.canonicalize!(gpu_stab; phases=true)
        cpu_canon = gpu_to_cpu_stabilizer(gpu_stab)
    end

    @test cpu_canon != cpu_stab  
    println("Passed test_random_stabilizer.")
end
function run_all_tests()
    println("\nStarting all tests...\n")

    @testset "GpuCanon Tests" begin
        @testset "Copy Kernel Test" begin
            test_copy_kernel()
        end
        @testset "Swap Columns Test" begin
            test_swap_columns()
        end
        @testset "Swap Columns and Phases Test" begin
            test_swap_columns_and_phases()
        end
        @testset "Pivot Search Test" begin
            test_pivot_search()
        end
        @testset "Update Kernel Test" begin
            test_update_kernel()
        end
        @testset "Canonicalization Test" begin
            test_canonicalization()
        end
        @testset "Conversion Functions Test" begin
            test_conversion_functions()
        end
        @testset "CPU-GPU Stabilizer Equivalence Test" begin
            test_cpu_gpu_stabilizer_equivalence()
        end
        @testset "Canonicalize! Test" begin
            test_canonicalize!()
        end
        @testset "Random Stabilizer Test" begin
            test_random_stabilizer()
        end
    end

    println("\nAll tests completed!")
end


run_all_tests()

