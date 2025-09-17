
#=============================================================================#
# The output should be both stylish and informative.
const plot_style = Dict(
    :xticks => sizes_MiB,
    :xscale => :log2,
    :shape => :circle,
    :background_color => :transparent
    )

# (La)TeX hates SVG but the Plots package has issues with transparent PDFs.
const file_format = "svg"
#=============================================================================#
