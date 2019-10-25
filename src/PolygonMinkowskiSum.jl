module PolygonMinkowskiSum
export minkowski_sum

using GeometryTypes: Point
using LinearAlgebra: normalize, dot

struct SimplePolygon{VP}
    points::VP
end

using RecipesBase: @recipe

@recipe function f(p::SimplePolygon)
    cycle = p.points[[1:end...,1]]
    seriestype --> :shape
    (x, y) = ([p[1] for p in cycle], [p[2] for p in cycle])
    x, y
end

# TODO doesn't work for `Array`. What is a generic constructor?
rot_quarter(v::P) where {P} = P(-v[1], v[0]) # TODO constrain P

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
function minkowski_sum(this::SimplePolygon, alongthis::SimplePolygon)
    P = typeof(this.points[1] + alongthis.points[1]) # TODO infer eltype
    result = P[]
    n = 1

    d = normalize(alongthis.points[n + 1] - alongthis.points[n])
    # TODO don't create an array for these differences.
    this_diff = [normalize(this.points[mod1(n + 1, end)] - this.points[n]) for n=1:length(this.points)]
    #@show collect(cross(this_diff[i], this_diff[mod1(i + 1, end)])>0 for i = 1:length(this_diff))
    @assert all(cross(this_diff[i], this_diff[mod1(i + 1, end)])>0 for i = 1:length(this_diff))
    @show this_diff

    @show d
    m = mod1(arg_next_direction(this_diff, d) + 1, length(this_diff))
    #@show this_diff[m]
    #m = argmax(map(p->dot(p, d), this.points))
    push!(result, this.points[m] + alongthis.points[n])
    #n += 1
    while length(result) < length(this.points) + length(alongthis.points)
        @show length(result)
        @show n, m
        # TODO avoid recomputing
        # TODO check for repeated points (causes div-by-zero in `normalize`)
        d_path = normalize(alongthis.points[mod1(n+1, end)] - alongthis.points[n])
        d_poly = normalize(this.points[mod1(m+1, end)] - this.points[mod1(m, end)])
        @show d_path
        @show d_poly
        @show d
        @assert cross(d, d_path) >= 0
        @assert cross(d, d_poly) >= 0
        if cross(d, d_path) < cross(d, d_poly)
            n += 1
            d = d_path
        else
            m += 1
            d = d_poly
        end
        push!(result, this.points[mod1(m, end)] + alongthis.points[mod1(n, end)])
    end
    return SimplePolygon(result)
end
end # module
