#!/usr/bin/env julia

# File description -------------------------------------------------------------
#
# Variational parameter updates for latent dirichlet allocation
#
# author: sankaran.kris@gmail.com
# date:

# libraries
using Distributions

###############################################################################
#                       simulate toy data for debugging                       #
###############################################################################

"""
# Examples

```julia-repl
theta, beta = simulate_params(100, 10, 3)
n = simulate_counts(beta, theta, 100)
```
"""
function simulate_counts(beta::Matrix, theta::Matrix, N::Int)
  V, K = size(beta)
  D, _ = size(theta)

  n = zeros(Int, D, V)
  for i = 1:D
    n[i, :] = rand(Multinomial(N, beta * theta[i, :]))
  end

  return n
end

function simulate_params(D::Int, V::Int, K::Int)
  theta = zeros(D, K)
  beta = zeros(V, K)

  for i = 1:D
    theta[i, :] = rand(Dirichlet(ones(K)))
  end

  for k = 1:K
    beta[:, k] = rand(Dirichlet(ones(V)))
  end

  return theta, beta
end

###############################################################################
#                         Overall variational updates                         #
###############################################################################

function variational_theta(alpha::Vector, n::Matrix, vc::Array)
  D, V = size(n)
  K = length(alpha)

  vtheta = zeros(D, K)
  for i = 1:D
    for k = 1:K
      vtheta[i, k] = alpha[k] + dot(n[i, :], vc[i, :, k])
    end
  end

  return vtheta
end

function variational_c(b::Matrix, vtheta::Matrix)
  V, K = size(b)
  D, _ = size(vtheta)

  c = zeros(D, V, K)
  for i = 1:D
    for v = 1:V
      for k = 1:K
        c[i, v, k] = b[v, k] *  exp(digamma(vtheta[i, k]) - digamma(sum(vtheta[i, :])))
      end
      c[i, v, :] = c[i, v, :] / sum(c[i, v, :])
    end
  end

  return c
end

function update_b(gamma::Vector, n::Matrix, c::Array)
  D, V, K = size(c)

  b = zeros(V, K)
  for k = 1:K
    for v = 1:V
      b[v, k] = gamma[v] + dot(n[:, v], c[:, v, k])
    end
    b[:, k] = b[:, k] / sum(b[:, k])
  end

  return b
end

function initialize_c(K::Int, n::Matrix)
  D, V = size(n)
  c = zeros(D, V, K)
  for i = 1:D
    for v = 1:V
      c[i, v, :] = rand(Multinomial(n[i, v], ones(K) / K))
    end
  end
  return c
end

"""
# Examples
```julia-repl
alpha = ones(3)
gamma = ones(10)
theta, beta = simulate_params(1000, 10, 3)
n = simulate_counts(beta, theta, 1000)
vb_fit = vb(n, alpha, gamma)
cor(vb_fit[3], beta)
cor(vb_fit[1], theta)
```
"""
function vb(n::Matrix, alpha::Vector, gamma::Vector, n_iter::Int = 10)
  D, V = size(n)
  K = length(alpha)

  vtheta = rand(Dirichlet(alpha), D)'
  vc = initialize_c(K, n)
  b = rand(Dirichlet(gamma), K)

  for iter = 1:n_iter
    vtheta = variational_theta(alpha, n, vc)
    vc = variational_c(b, vtheta)
    b = update_b(gamma, n, vc)
  end

  return vtheta, vc, b
end
