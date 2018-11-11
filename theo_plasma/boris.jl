using LinearAlgebra


function E(r :: Vector{Float64}, t :: Float64) :: Vector{Float64}
    [0., 0., .01]
end

function B(r :: Vector{Float64}, t :: Float64) :: Vector{Float64}
    [0., 0., 1.]
end

function boris_step!(r :: Vector{Float64}, v :: Vector{Float64},
                     qm :: Float64, dt :: Float64, t :: Float64) :: Tuple{Vector, Vector}
    dtqm = dt * qm
    r += .5 * dt .* v
    B_rt = B(r, t)
    E_rt = E(r, t)
    p = .5 * dtqm .* B_rt
    a_sq = .25 * dtqm^2 * dot(B_rt, B_rt)
    v += .5 * dtqm .* E_rt
    v_prime = v .+ cross(v, p)
    v += 2. .* cross(v_prime, p) / (1 + a_sq)
    v += .5 * dtqm .* E_rt
    r += .5 * dt .* v
    return r, v
end

qm = 1.
dt = .25
steps = 100
r = zeros(steps, 3)
v = zeros(steps, 3)
r[1,:] = [0. 0. 0.]
v[1,:] = [.5 0. 0.]

function test!()
    t = 0.
    for i=2:steps
        t += .5 * dt
        next = boris_step!(r[i-1,:], v[i-1,:], qm, dt, t)
        r[i,:] = next[1]
        v[i,:] = next[2]
        t += .5 * dt
    end
    return r, v
end
