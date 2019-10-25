using PolygonMinkowskiSum
using Test

import Unitful: °
using GeometryTypes: Point
using PolygonMinkowskiSum: SimplePolygon

function make_regular_ngon(n; offset=0°)
    points = [Point(cos(θ+offset), sin(θ+offset))
        for θ = range(0°, 360°, length=(n+1))][1:end-1]
    SimplePolygon(points)
end

using Plots

for offset = 0°:5°:360°
    @show offset
    diamond = make_regular_ngon(4; offset=0°)
    square = make_regular_ngon(4; offset=offset+45°)


    octogon = minkowski_sum(square, diamond)
    p = plot(aspect_ratio=1.0)
    plot!(p ,diamond, fillalpha=0.2)
    plot!(p, square, fillalpha=0.2)
    plot!(p, octogon, fillalpha=0.2)
end
