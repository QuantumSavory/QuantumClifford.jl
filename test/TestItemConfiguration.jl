module TestItemConfiguration

export JULIATI_CI_FILTER, JULIATI_PLATFORM_FILTER, platform_supported

# Bit-packing is supported only on 64-bit Linux.
function platform_supported(tags)
    return (!(:bitpack in tags) || (Sys.islinux() && Int === Int64)) &&
           (!(:julia_1_14_required in tags) || VERSION >= v"1.14") &&
           (!(:oscar_required in tags) || (!Sys.iswindows() && Sys.ARCH == :x86_64)) &&
           (!(:tesseract_required in tags) || !Sys.iswindows())
end

const JULIATI_PLATFORM_FILTER = """
    (!(:bitpack in tags) || (Sys.islinux() && Int === Int64)) &&
    (!(:julia_1_14_required in tags) || VERSION >= v\"1.14\") &&
    (!(:oscar_required in tags) || (!Sys.iswindows() && Sys.ARCH == :x86_64)) &&
    (!(:tesseract_required in tags) || !Sys.iswindows())
    """

const JULIATI_CI_FILTER = "($JULIATI_PLATFORM_FILTER) && !(:too_slow_for_ci in tags)"

end
