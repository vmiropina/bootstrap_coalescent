from simulations_functions import *

#### This code is minimal code to generate figures simular to the ones in the article


##### Choice of the parameters
alpha = 1.2  # characteristic parameter of the beta-coalescent and must be between 1 and 2 to satisfy msprime's assumptions
m = 10 # number of individuals in the sample
n = 100 # number of independent realizations
L = 20 # lenght of the genome



# Generate the matrices A and B
# We want to solve AX = B where B = (c_i) and the ith column of A  = E(Mi(c0)) - E(Mi(bj)) otherwise
# Equation 13 in Proposition 3.1 in the paper 

A, B = estimate_AB(n, m)

##### This can be uncommented to save the matrices A and B and to load them 
# np.savetxt('A_alpha_' + str(alpha) + '_m_' + str(m) + '_n_' + str(n) +  '_L_' + str(L) + '.txt',A)
# np.savetxt('B_alpha_' + str(alpha) + '_m_' + str(m) + '_n_' + str(n) +  '_L_' + str(L) + '.txt',B)
# A = np.loadtxt('A_alpha_' + str(alpha) + '_m_' + str(m) + '_n_' + str(n) +  '_L_' + str(L) + '.txt')
# B = np.loadtxt('B_alpha_' + str(alpha) + '_m_' + str(m) + '_n_' + str(n) +  '_L_' + str(L) + '.txt')


# Solve the AX=B where X contains the estimated coalescent rates (rmj's)
bnds = [(0, 100) for _ in B]
R = solve_constraints(A, B, bnds)


# Compute the true coalescent rates (rmj's)
trueR = coal_rates(m, np.linspace(2, m, m-1), alpha)

# Plot the real nu measure, the one reconstructed from the true coalescent rates (using the reconstruction techniques of the first section)
# and the nu measure reconstructed from the estimated coalescent rates
h = 0.1
plot(h, alpha, m, R, trueR, 1000, 1.5)

