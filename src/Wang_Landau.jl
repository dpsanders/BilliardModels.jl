# Try to use Wang-Landau to estimate "density of states" of "energy"=path length

using BilliardModels

function find_free_path(billiard_table, x, v)
    l = Vector2D(0, 0)
    p = ParticleOnLattice(x, v, l)
    num_collisions = 1
    xs, ls, free_paths = billiard_dynamics_on_lattice(p, billiard_table, num_collisions)

    free_paths[1]
end

unif_rand(delta) = delta*(rand()-0.5)

function explore_phase_space(radius=0.1, delta=0.1, num_steps=1000)

    billiard_table = Sinai_billiard(radius, true, true)  # periodic in x and y

    xx, vv = initial_condition(billiard_table, -.5, .5, -.5, .5)

    theta = atan2(vv[2], vv[1])  # this is how I generated the velocity to start with (from an angle) -- a bit ridiculous

    xs = [xx[1]]
    ys = [xx[2]]
    thetas = [theta]
    free_paths = [find_free_path(billiard_table, xx, vv)]

    for i in 1:num_steps
        x, y = xx[1], xx[2]

        x += unif_rand(delta)
        y += unif_rand(delta)

        # periodize:
        (x > 0.5) && (x -= 1)
        (x < -0.5) && (x += 1)
        (y > 0.5) && (y -= 1)
        (y < -0.5) && (y += 1)

        # check outside disc:
        if x*x + y*y < radius^2
            continue
        end

        xx_new = Vector2D(x, y)

        theta += unif_rand(delta)

        vv = Vector2D(cos(theta), sin(theta))


        free_path = find_free_path(billiard_table, xx_new, vv)

        push!(xs, x)
        push!(ys, y)
        push!(thetas, theta)
        push!(free_paths, free_path)

        xx = xx_new

    end

    xs, ys, thetas, free_paths
end



