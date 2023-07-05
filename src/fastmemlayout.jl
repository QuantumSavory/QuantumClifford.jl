"""Convert a tableau to a memory layout that is fast for row operations (e.g. canonicalization or multiplication with a dense `CliffordOperator`).

```jldoctest
julia> s, c = random_stabilizer(300), random_clifford(300);

julia> sc, sr = fastcolumn(copy(s)), fastrow(copy(s));

julia> tc, tr = (@benchmark apply!(sc,c)), (@benchmark apply!(sr,c));

julia> minimum(tc).time > minimum(tr).time
true
```

See also: [`fastrow`](@ref)"""
function fastrow end

"""Convert a tableau to a memory layout that is fast for column operations (e.g. sparse gate application).

See also: [`fastrow`](@ref)"""
function fastcolumn end

fastrow(t::Tableau{Tzv,Tm}) where {Tzv, Tm} = t
fastrow(t::Tableau{Tzv,Tm}) where {Tzv, Tm<:Adjoint} = Tableau(t.phases, t.nqubits, collect(t.xzs))
fastcolumn(t::Tableau{Tzv,Tm}) where {Tzv, Tm} = Tableau(t.phases, t.nqubits, collect(t.xzs')')
fastcolumn(t::Tableau{Tzv,Tm}) where {Tzv, Tm<:Adjoint} = t

fastrow(s::Stabilizer) = Stabilizer(fastrow(s.tab))
fastcolumn(s::Stabilizer) = Stabilizer(fastcolumn(s.tab))

fastrow(s::Destabilizer) = Destabilizer(fastrow(s.tab))
fastcolumn(s::Destabilizer) = Destabilizer(fastcolumn(s.tab))

fastrow(s::MixedStabilizer) = MixedStabilizer(fastrow(s.tab), s.rank)
fastcolumn(s::MixedStabilizer) = MixedStabilizer(fastcolumn(s.tab), s.rank)

fastrow(s::MixedDestabilizer) = MixedDestabilizer(fastrow(s.tab), s.rank)
fastcolumn(s::MixedDestabilizer) = MixedDestabilizer(fastcolumn(s.tab), s.rank)

fastrow(s::PauliFrame) = PauliFrame(fastrow(s.frame), s.measurements)
fastcolumn(s::PauliFrame) = PauliFrame(fastcolumn(s.frame), s.measurements)
