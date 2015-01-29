using FactCheck

push!(LOAD_PATH, "/Users/dsanders/Dropbox/papers/billiards/lorentz_gas/")
#println(LOAD_PATH)

using Lorentz
using Vector2d

facts("Particle tests") do

    p = Particle([1, 2], [3, 4])
    @fact p => Particle(Vector2D{Int64}(1,2), Vector2D{Int64}(3,4))

end


facts("Disc tests") do
    d = Disc([1, 2], 3)
    @fact d => Disc(Vector2D{Int64}(1,2),3.0)

    d = Disc([0, 0], 1)
    p = Particle([-2, 0], [1, 0])

    @fact collision_time(p, d)

end


