using FactCheck

push!(LOAD_PATH, "/Users/dsanders/Dropbox/papers/billiards/lorentz_gas/")
#println(LOAD_PATH)

using Vector2d

facts("Vector2D tests") do
    @fact Vector2D(1, 2) => Vector2D{Int}(1, 2)
    @fact Vector2D([1., 2.]) => Vector2D{Float64}(1., 2.)

    a = Vector2D(1, 2)
    b = Vector2D(3, 4)

    @fact a + b => Vector2D(4, 6)
    @fact a ⋅ b => 11
end
