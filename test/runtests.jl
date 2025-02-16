using QuantumClifford
using TestItemRunner

if get(ENV, "GPU_TESTS", "") != "true"
    println("skipping gpu tests (set GPU_TESTS=true to test gpu)")
end

# filter for the test
testfilter = ti -> begin
  exclude = Symbol[]
  if get(ENV,"JET_TEST","")!="true"
    push!(exclude, :jet)
  end
  if !(VERSION >= v"1.10")
    push!(exclude, :doctests)
    push!(exclude, :aqua)
  end

  if get(ENV, "GPU_TESTS", "")!="true"
    push!(exclude, :gpu)
  end

  if !(Base.Sys.islinux() & (Int===Int64))
    push!(exclude, :bitpack)
  end

  return all(in([:alloccc]), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter
