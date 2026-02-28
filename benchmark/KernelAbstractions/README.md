# Usage Instructions

1. Modify the values listed in `implementation/definitions/[benchmark_configuration, tuning_parameters].jl` as desired, keeping in mind the limitations of the available device memory and the extent of the benchmark runtime.
2. Ensure that all the packages listed in `implementation/imports.jl` are installed as this is not handled automatically.
3. Ensure that the backend package(s) listed in the pertinent `benchmark_platform_*.jl` script are properly setup and configured.
4. Pass said script as an argument to julia (optionally, also set the number of executing host threads) and await for the benchmark to conclude.
5. Navigate to the newly created directory (hierarchy) and find the matching platform/timestamp pair in order to inspect the results.

# Noteworthy Details

- Even though OpenCL should in principle have cross-vendor compatibility, the JuliaGPU toolchain requires SPIR-V support in order to compile native julia code. At the present moment, this benchmark runs on the host and relies on `pocl_jll` to provide this intermediate representation. In the event that this is not necessary and/or desired, the line importing it may simply be removed from the OpenCL script.
- Modern hardware is able to dynamically change its clock speed depending on its current workload and operating conditions, possibly boosting performance in certain situations. This is *antithetical* to the reliability and reproducibility of any benchmark results since it sullies the measured performance metrics. As such, it is advised to consult with the relevant hardware manufacturer and/or operating system documentation for instructions on how to disable such features, or at the very least ensuring that they display minimal variation.
