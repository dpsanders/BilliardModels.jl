using BilliardModels

include("StatisticsAccumulator.jl")

function run_lorentz(num_particles, num_collisions, delta_t=5.0)

    billiard_table = Sinai_billiard(0.354, true, true)  # periodic in x and y

    stats_usual = StatisticsObject()  # usual statistics
    stats_all_origins = StatisticsObject()  # use all times as origin

    skip = 5

    for i in 1:num_particles
        x0, v0 = initial_condition(billiard_table, -.5, .5, -.5, .5)
        l = Vector2D(0, 0)
        p = ParticleOnLattice(x0, v0, l)

        xs, ls, free_paths = billiard_dynamics_on_lattice(p, billiard_table, num_collisions)
        positions, times = continuous_time(xs, ls, free_paths, delta_t)

        # "usual" statistics: use each time and position only once
        for (t, xx) in zip(times, positions)
            #r = norm(xx - x0)
            r = xx[1] - x0[1]
            add_data!(stats_usual, t, r)
        end

        # "all" statistics: use each (or nearly each) time and position as a new origin
        for i in 1:skip:length(times)
            t00 = times[i]
            x00 = positions[i][1]

            for j in i+1:length(times) #skip:length(times)

                time_displacement = times[j] - t00
                displacement = positions[j][1] - x00

                #r = norm(displacement)
                r = displacement
                add_data!(stats_all_origins, time_displacement, r)
            end
        end

    end

    get_statistics(stats_usual), get_statistics(stats_all_origins)
end


