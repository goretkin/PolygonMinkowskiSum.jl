module PolygonMinkowskiSum
export minkowski_sum

using GeometryTypes: Point
using LinearAlgebra: normalize, dot

"""
Polygon with no holes
"""
struct SimplePolygon{VP}
    points::VP
end

function diff(p::SimplePolygon)
    [p.points[mod1(n + 1, end)] - p.points[n] for n=1:length(p.points)]
end

function crosses(p::SimplePolygon)
    d = diff(p)
    [cross(d[i], d[mod1(i+1, end)]) for i=1:length(d)]
end

function find_ccw_violations(p::SimplePolygon)
    # TODO have option for strict or not (allow co-linear points?)
    cs = crosses(p)
    bad_i = findall(cs .< 0) # TODO generic zero
    return bad_i
end

is_ccw(p1, p2, p3) = cross(p1, p2) ≥ 0 && cross(p2, p3) ≥ 0

# @recipe macro requires RecipesBase to be defined in scope of macro use
using RecipesBase: RecipesBase, @recipe

@recipe function f(p::SimplePolygon)
    cycle = p.points[[1:end...,1]]
    seriestype --> :shape
    (x, y) = ([p[1] for p in cycle], [p[2] for p in cycle])
    x, y
end

# TODO doesn't work for `Array`. What is a generic constructor?
rot_quarter(v::P) where {P} = P(-v[2], v[1]) # TODO constrain P

# asin(angle from v1 to v2) * |v1| * |v2|
cross(v1, v2) = v1[1] * v2[2] - v1[2] * v2[1]

"""
find the vector that points in the same direction (positive dot product),
is ccw from v (positive cross product), that is next (smallest cross product)
"""
function arg_next_direction(vs, v)
    I = eachindex(vs)
    i_best = nothing
    c_best = typemax(Float64) # TODO generic
    for i in I
        v_i = normalize(vs[i])
        if dot(v_i, v) > 0
            c = cross(v_i, v)
            if c > 0 && c < c_best
                c_best = c
                i_best = i
            end
        end
    end
    return i_best
end

"""
Assume all point sequences are arranged counter-clockwise
Assume both shapes are convex


Effectively https://en.wikipedia.org/wiki/Minkowski_addition#Two_convex_polygons_in_the_plane
"""
function minkowski_sum(this::SimplePolygon, alongthis::SimplePolygon, plotdebug=nothing)
    IJ = find_minkowski_sum(this, alongthis)
    r = SimplePolygon([this.points[i] + alongthis.points[j] for (i,j) in IJ])
    if plotdebug !== nothing
        ss = [string((k,i,j)) for (k, (i, j)) in enumerate(IJ)]
        plotdebug.ns.scatter!(plotdebug.plots["result"], r.points, series_annotation=ss)
        plotdebug.ns.scatter!(plotdebug.plots["poly1"], this.points, series_annotation=map(string,1:length(this.points)))
        plotdebug.ns.scatter!(plotdebug.plots["poly2"], alongthis.points, series_annotation=map(string,1:length(alongthis.points)))
    end
    return r
end

"""
only works if there are no co-linear edges
"""
function sort_find_sort_polar(p1::SimplePolygon, p2::SimplePolygon)
    polar_i(p) = [(angle=atan(y,x), index=i) for (i,(x,y)) in enumerate(diff(p))]

    (s1, s2) = ([(vertex=pi, index_polygon=1) for pi in polar_i(p)] for p in (p1, p2))
    # TODO these sorts need to be stable with repsect to cyclic order, and they
    s1 = sort(s1)
    s2 = sort(s2)
    s1, s2, sort(vcat(s1, s2))
end

function find_minkowski_sum2(p1::SimplePolygon, p2::SimplePolygon)
    ps = (p1, p2)
    s1, s2, sboth = sort_find_sort_polar(ps...)

    if first(sboth) == first(s2)
        return [(i2, i1) for (i1, i2) in find_minkowski_sum(p2, p1)]
    end
    @assert first(sboth) == first(s1)

    r = Tuple{Int64, Int64}[]

    is = [s1[1].vertex.index, s2[end].vertex.index]

    j = 2
    while true
        push!(r, tuple(is...))
        if !(j <= length(sboth)); break; end
        s = sboth[j]
        is[s.index_polygon] = s.vertex.index
        j += 1
    end
    is = [s1[end].vertex.index, s2[1].vertex.index]
    push!(r, tuple(is...))
    return r
end

function find_minkowski_sum(this::SimplePolygon, alongthis::SimplePolygon)
    result = Tuple{Int64,Int64}[]
    n = 1

    d = alongthis.points[n + 1] - alongthis.points[n]
    # TODO don't create an array for these differences.
    this_diff = diff(this)
    @assert all(cross(this_diff[i], this_diff[mod1(i + 1, end)])>0 for i = 1:length(this_diff))

    dots = [dot(d, d1) for d1 in this_diff]
    flips = [dots[mod1(i+1, end)] * dots[i] for i in 1:length(dots)]
    # generic zero (if units)
    # TODO avoid collect
    m = findfirst(((d,f),)->(d>=zero(d) && f<=zero(f)), collect(zip(dots, flips)))
    #@show m
    push!(result, (m, n))

    while length(result) < length(this.points) + length(alongthis.points)
        @show n, m
        # TODO avoid recomputing
        # TODO check for repeated points (causes div-by-zero in `normalize`)
        d_path = normalize(alongthis.points[mod1(n+1, end)] - alongthis.points[n])
        d_poly = normalize(this.points[mod1(m+1, end)] - this.points[mod1(m, end)])
        # TODO these assertions are being violated, understand why they are wrong or why the algorithm is wrong
        #@assert cross(d, d_path) >= 0
        #@show cross(d, d_poly)
        #@assert cross(d, d_poly) >= 0
        if cross(d, d_path) < cross(d, d_poly)
            n += 1
            d = d_path
        else
            m += 1
            d = d_poly
        end
        push!(result,
            (mod1(m, length(this.points)), mod1(n, length(alongthis.points)))
        )
    end
    return result
end
end # module
