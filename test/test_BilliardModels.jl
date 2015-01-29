using FactCheck

#push!(LOAD_PATH, "/Users/dsanders/Dropbox/papers/billiards/")

using BilliardModels
#using Vector2d

facts("Particle tests") do

    p = Particle([1, 2.], [3., 4.5])
    @fact p.x => Vector2D{Float64}(1,2)
    @fact p.v => Vector2D{Float64}(3.0,4.5)

end


facts("Disc tests") do
    d = Disc([1, 2], 3)
    @fact d.radius => 3.0

    d = Disc([0, 0], 1)
    p = Particle([-2, 0], [1, 0])

    @fact collision_time(p, d) => 1

end


