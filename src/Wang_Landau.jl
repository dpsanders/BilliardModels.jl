# Try to use Wang-Landau to estimate "density of states" of "energy"=path length

using BilliardModels

function find_free_path(billiard_table, x, v)
    l = Vector2D(0, 0)
    p = ParticleOnLattice(x, v, l)
    num_collisions = 1
    xs, ls, free_paths = billiard_dynamics_on_lattice(p, billiard_table, num_collisions)

    free_paths[1]
end


function explore_phase_space(radius=0.1, delta=0.1, num_steps=1000)

    billiard_table = Sinai_billiard(radius, true, true)  # periodic in x and y

    xx, vv = initial_condition(billiard_table, -.5, .5, -.5, .5)

    #theta = atan2(vv[2], vv[1])  # this is how I generated the velocity to start with (from an angle) -- a bit ridiculous

    xx_old, vv_old = xx, vv

    xs = [xx[1]]
    ys = [xx[2]]

    for i in 1:num_steps
        x, y = xx[1], xx[2]

        x += delta*(rand()-0.5)
        y += delta*(rand()-0.5)

        # periodize:
        (x > 0.5) && (x -= 1)
        (x < -0.5) && (x += 1)
        (y > 0.5) && (y -= 1)
        (y < -0.5) && (y += 1)

        # check outside disc:
        if x*x + y*y < radius^2
            continue
        end

        push!(xs, x)
        push!(ys, y)

        xx = [x, y]

    end

    xs, ys
end


#         theta = mod2pi(theta + delta)
#         vv = (cos(theta), sin(theta))
#         find_free_path(billiard_table, x, v)
