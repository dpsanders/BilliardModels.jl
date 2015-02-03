using BilliardModels

billiard_table = Sinai_billiard(0.354,true,true)

x = Vector2D(0.3, 0.)
v = Vector2D(0.99, 0.01)
l = Vector2D(0, 0)
p = ParticleOnLattice(x, v, l)

N = 1_000_000
xs, ls, free_paths = billiard_dynamics_on_lattice(p, billiard_table, 1);
@time xs, ls, free_paths = billiard_dynamics_on_lattice(p, billiard_table, N);
