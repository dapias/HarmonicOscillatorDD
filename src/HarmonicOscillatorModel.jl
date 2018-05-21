module HarmonicOscillatorModel

export initialize, computeenergyandforces!, measurekineticenergy, Atom, dim

const dim = 1

type Atom{T}
  r::Array{T,1}
  p::Array{T,1}
  f::Array{T,1}
end

Atom(r) = Atom(r,zeros(dim), zeros(dim))


function initialize(N::Int64, T::Float64)
  L = 1.0
  atoms = Array{Atom,1}(N)

  oscillator = Atom([rand()], [rand()], [0.0])

  atoms[1] = oscillator

  K = 0 #Kinetic energy

  for i in 1:N
    for j in 1:dim
      K += atoms[i].p[j]^2
    end
  end

  U = computeenergyandforces!(atoms)

  #Instantaneous kinetic temperature and energy
  T = K/(dim*(N))
  K = K/2

  return atoms, T, K, U, L
end


@doc """ Determine the interaction force for each pair of particles (i, j)"""->
function computeenergyandforces!(atoms::Array{Atom,1})

  U = 0.0
  atoms[1].f[1] = 0.
  U += (atoms[1].r[1])^2./2.
  atoms[1].f[1] = -atoms[1].r[1]

  return U
end


function measurekineticenergy(atoms::Array{Atom,1})
  K = 0.
  N = length(atoms)
  for i in 1:N
    for j in 1:dim
      K +=  atoms[i].p[j]^2.
    end
  end

  K = K/2.
  return K
end

end
