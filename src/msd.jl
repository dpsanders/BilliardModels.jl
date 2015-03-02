using BilliardModels

include("StatisticsAccumulator.jl")

function run(N)

    billiard_table = Sinai_billiard(0.354, true, true)  # periodic in x and y

    num_collisions = 1000
    delta_t = 10

    stats = StatisticsObject()

    for i in 1:N
        x, v = initial_condition(billiard_table, -.5, .5, -.5, .5)
        l = Vector2D(0, 0)
        p = ParticleOnLattice(x, v, l)

        xs, ls, free_paths = billiard_dynamics_on_lattice(p, billiard_table, num_collisions)
        positions, times = continuous_time(xs, ls, free_paths, delta_t)

        for (t, xx) in zip(times, positions)
            r = norm(xx)
            add_data!(stats, t, r)
        end
    end

    get_statistics(stats)
end


