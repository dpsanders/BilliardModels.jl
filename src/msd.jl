using BilliardModels

include("StatisticsAccumulator.jl")

function run_lorentz(num_particles, num_collisions, delta_t=5.0)

    billiard_table = Sinai_billiard(0.354, true, true)  # periodic in x and y

    stats = StatisticsObject[]

    num_stats_objects = 10
    for i in 1:num_stats_objects+1
        push!(stats, StatisticsObject())  # use all times as origin
    end

    # Divide statistics up according to first free path


    skip = 5

    for i in 1:num_particles
        x0, v0 = initial_condition(billiard_table, -.5, .5, -.5, .5)
        l = Vector2D(0, 0)
        p = ParticleOnLattice(x0, v0, l)

        xs, ls, free_paths = billiard_dynamics_on_lattice(p, billiard_table, num_collisions)
        positions, times = continuous_time(xs, ls, free_paths, delta_t)

        first_free_path = free_paths[1]

        if first_free_path < num_stats_objects
            stats_object = stats[ceil(first_free_path)]
        else
            stats_object = stats[end]
        end


        # "all" statistics: use each (or nearly each) time and position as a new origin
        for i in 1:skip:length(times)
            t00 = times[i]
            x00 = positions[i][1]  # x component

            for j in i+1:length(times) #skip:length(times)

                time_displacement = times[j] - t00
                displacement = positions[j][1] - x00

                #r = norm(displacement)
                #r = displacement
                add_data!(stats_object, time_displacement, displacement)
            end
        end

    end

    map(get_statistics, stats)

end


