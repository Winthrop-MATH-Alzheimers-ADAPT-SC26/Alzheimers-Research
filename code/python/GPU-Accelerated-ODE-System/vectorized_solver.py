import jax
import jax.numpy as jnp
import diffrax
from diffrax import ODETerm, Tsit5, SaveAt, PIDController



# alzheimers ode system (no treatment terms)
def alzsys(t, y, args):
    Ab, Ca, tau, N, C = y
    a1, a2, b1, b2, c1, c2, c3, d1, d2, e1, e2, e3, k1, k2, k3, k4, k5, sig1, sig2 = args

    dAb_dt = a1 + a2 * (Ca / (Ca + sig1)) - k1 * Ab
    dCa_dt = b1 + b2 * Ab - k2 * Ca
    dtau_dt = c1 + c2 * Ab + c3 * (Ca / (Ca + sig2)) - k3 * tau
    dN_dt = d1 + d2 * tau - k4 * N
    dC_dt = e1 + e2 * N + e3 * tau - k5 * C

    return jnp.stack([dAb_dt, dCa_dt, dtau_dt, dN_dt, dC_dt])


# solver setup
term = ODETerm(alzsys)
solver = Tsit5() # non stiff solver. stiff is Kvaerno5
t0 = 0
tf = 50
saveat = SaveAt(ts = jnp.linspace(t0, tf, 1000))
stepsize_controller = diffrax.PIDController(rtol = 1e-5, atol = 1e-5) # need for stiff system


# vectorized ODE solving
# takes in initial conditions and a batch of parameters, and solves the ODE system for each set of parameters
@jax.jit
def solve_batch(y0, args_batch):
    fn = jax.vmap(lambda args: diffrax.diffeqsolve(
        terms = term, 
        solver = solver,   
        t0 = t0,
        t1 = tf,
        dt0 = None, # None -> picks automatic dt0
        y0 = y0,
        args = args,
        saveat = saveat,
        stepsize_controller = stepsize_controller,
        max_steps = 100000 # arbitrary large number
    ))

    return fn(args_batch)

# note: args_batch must be a jnp array