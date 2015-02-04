using BilliardModels

function run_sinai(num_particles)
    billiard_table = Sinai_billiard(0.354, true, true)  # periodic in x and y

    x = Vector2D(0.3, 0.)
    v = Vector2D(0.99, 0.01)
    l = Vector2D(0, 0)
    p = ParticleOnLattice(x, v, l)

    xs, ls, free_paths = billiard_dynamics_on_lattice(p, billiard_table, num_particles);
end


N = 0
try
    N = int(ARGS[1])
catch
    N = 1_000_000
end

println("# Using N=$N")

run_sinai(1)
@time run_sinai(N)
