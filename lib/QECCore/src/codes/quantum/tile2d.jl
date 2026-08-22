"""
`Tile2D` is a generalization of surface codes that offers flexibility in terms of *locality* and stabilizer *check weight*
without compromising on the 2D locality of the 2D surface code [steffan2025tilecodeshighefficiencyquantum](@cite). It encodes more
logical qubits than surface code, and provides ``\\mathcal{O}(1)``-locality.

The `horiz` and `vert` arguments give zero-based edge coordinates in the X-stabilizer tile, with both coordinates in `0:B-1`.
The Z-stabilizer tile is derived automatically by reflection and exchange of the horizontal and vertical supports. This
implementation constructs the unrotated rectangular layout with `Lx` by `Ly` tiles.

Here is an example of weight-6 `[[288, 8, 12]]` 2D Tile code from Table I of [steffan2025tilecodeshighefficiencyquantum](@cite).

!!! note
    The `[[288, 8, 12]]` code was first discovered by Liang *et al.* [liang2025planar](@cite) in the context of constructing planar
    quantum low-density parity-check codes with *open* boundary conditions.

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> B = 3;

julia> horizX = [(0,0), (2,1), (2,2)];

julia> vertX = [(0,2), (1,2), (2,0)];

julia> Lx, Ly = 10, 10;

julia> c = Tile2D(B, horizX, vertX, Lx, Ly);

julia> code_n(c), code_k(c)
(288, 8)
```
"""
struct Tile2D <: AbstractCSSCode
    """Size of the tile box ``(B \\times B)`` determining the support of a stabilizer."""
    B::Int
    """Zero-based positions of horizontal edges within the tile box."""
    horiz::Vector{Tuple{Int,Int}}
    """Zero-based positions of vertical edges within the tile box."""
    vert::Vector{Tuple{Int,Int}}
    """Number of tiles along the x-direction."""
    Lx::Int
    """Number of tiles along the y-direction."""
    Ly::Int

    function Tile2D(B::Int, horiz::Vector{Tuple{Int,Int}}, vert::Vector{Tuple{Int,Int}}, Lx::Int, Ly::Int)
        (B > 0 && Lx > 0 && Ly > 0) || throw(ArgumentError("B, Lx, Ly must be positive"))
        for (x, y) in Iterators.flatten((horiz, vert))
            (0 <= x < B && 0 <= y < B) || throw(ArgumentError("Tile coordinates must satisfy 0 <= x, y < B"))
        end
        new(B, copy(horiz), copy(vert), Lx, Ly)
    end
end

function _rectangular_layout(tile::Tile2D)
    black = Tuple{Int,Int}[]
    red = Tuple{Int,Int}[]
    blue = Tuple{Int,Int}[]
    B, Lx, Ly = tile.B, tile.Lx, tile.Ly
    for x in 0:Lx-1, y in 0:Ly-1
        push!(black, (x,y))
    end
    for x in 0:Lx-1, t in 1:B-1
        push!(red, (x, -t))
        push!(red, (x, Ly-1+t))
    end
    for y in 0:Ly-1, t in 1:B-1
        push!(blue, (-t, y))
        push!(blue, (Lx-1+t, y))
    end
    return black, red, blue
end

function _complement_tile(tile::Tile2D)
    B = tile.B
    horiz_z = [(B-1-x, B-1-y) for (x,y) in tile.vert]
    vert_z  = [(B-1-x, B-1-y) for (x,y) in tile.horiz]
    Tile2D(B, horiz_z, vert_z, tile.Lx, tile.Ly)
end

function _physical_qubits(tile::Tile2D)
    qubits = Set{Tuple{Symbol,Int,Int}}()
    black, _, _ = _rectangular_layout(tile)
    for (vx,vy) in black
        for x in 0:tile.B-1, y in 0:tile.B-1
            push!(qubits, (:h, vx+x, vy+y))
            push!(qubits, (:v, vx+x, vy+y))
        end
    end
    return qubits
end


function _edges((vx,vy)::Tuple{Int,Int}, tile::Tile2D)
    edges = Tuple{Symbol,Int,Int}[]
    for (x,y) in tile.horiz
        push!(edges, (:h, vx+x, vy+y))
    end
    for (x,y) in tile.vert
        push!(edges, (:v, vx+x, vy+y))
    end
    return edges
end

function parity_matrix_xz(tile::Tile2D)
    tileZ   = _complement_tile(tile)
    # "We will always restrict ourselves to (rotated) rectangular shapes" [steffan2025tilecodeshighefficiencyquantum](@cite).
    black, red, blue = _rectangular_layout(tile)
    physical = _physical_qubits(tile)
    Xrows = Vector{Vector{Tuple{Symbol,Int,Int}}}()
    Zrows = Vector{Vector{Tuple{Symbol,Int,Int}}}()
    # "We fine-tune the layout to the specific support of the stabilizers. First remove
    # all qubits that are not supported in any X-type stabilizer or are not supported
    # in any Z-type stabilizer" [steffan2025tilecodeshighefficiencyquantum](@cite).
    for v in black
        push!(Xrows, filter(in(physical), _edges(v, tile)))
        push!(Zrows, filter(in(physical), _edges(v, tileZ)))
    end
    for v in red
        push!(Xrows, filter(in(physical), _edges(v, tile)))
    end
    for v in blue
        push!(Zrows, filter(in(physical), _edges(v, tileZ)))
    end
    xqubits = Set(q for row in Xrows for q in row)
    zqubits = Set(q for row in Zrows for q in row)
    qubits = intersect(xqubits, zqubits)
    for rows in (Xrows, Zrows)
        for row in rows
            filter!(in(qubits), row)
        end
        # "Finally, we remove all stabilizers whose support has become empty because of aforementioned procedure" [steffan2025tilecodeshighefficiencyquantum](@cite).
        filter!(!isempty, rows)
    end
    qubits = sort!(collect(qubits); by=q -> (q[1] == :v, q[2], q[3]))
    qindex = Dict(q => i for (i,q) in enumerate(qubits))
    Hx = spzeros(Bool, length(Xrows), length(qubits))
    Hz = spzeros(Bool, length(Zrows), length(qubits))
    for (i,row) in enumerate(Xrows), q in row
        Hx[i, qindex[q]] = true
    end
    for (i,row) in enumerate(Zrows), q in row
        Hz[i, qindex[q]] = true
    end
    return Hx, Hz
end

parity_matrix_x(tile::Tile2D) = parity_matrix_xz(tile)[1]

parity_matrix_z(tile::Tile2D) = parity_matrix_xz(tile)[2]
