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

facts("Plane tests") do
    o = Plane([0,0], [1, 1])
    p = Particle([-1, 0], [1, 0])

    @fact collision_time(p, o) => 1
end

facts("Billiard table tests") do
    table = Sinai_billiard(0.1)
    p = Particle([0.2, 0.], [1., 0.])

    next_collision = calculate_next_collision(p, table, nothing)
    @fact next_collision[2] => Vector2D(-1., 0.)  # new velocity after reflecting on vertical boundary
    @fact next_collision[3] => 0.3  # col]\


    p = Particle([0.2, 0.], [-1., 0.])
    next_collision = calculate_next_collision(p, table, nothing)
    @fact isa(next_collision[4], Disc) => true

end



