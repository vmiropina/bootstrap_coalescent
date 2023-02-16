import msprime
import numpy as np
import scipy.linalg as ln
import scipy.special as sp
import scipy.stats as st
import scipy.optimize as optimize
import matplotlib.pyplot as plt



#### Generate SNP matrix from Lambda coalescent (Beta-coalescent) using msprime
def SNP_matrix(m, L, alpha):
    tree = msprime.sim_ancestry(samples=m, ploidy =1,  sequence_length=L,model = msprime.BetaCoalescent(alpha = alpha),population_size = 1000000)
    mts = msprime.sim_mutations(tree, rate=1, discrete_genome=False)
    return mts.genotype_matrix()

# weights matrix M, with weights bj
def weighted_sfs(K, m, j): 
    J = K.copy()
    alll = np.linspace(0, m-1,m, dtype = 'int')
    to_remove = np.random.choice(alll, j, replace = False) 
    to_multiply = np.random.choice(to_remove, 1)
    J[:, to_remove] = np.repeat(J[:, to_multiply][:, 0], j).reshape(len(J), j)
    return compute_sfs(J)

### Compute SFS from (weighted) SNP matrix
def compute_sfs(K):
    occ = np.sum(K, axis = 1)
    sfs, bins = np.histogram(occ, bins = np.linspace(1, m, m))  #with these bins, zeros are atomatically removed
    return sfs

#We want to solve AX = B where B = (c_i) and the ith column of A  = E(Mi(c0)) - E(Mi(bj)) otherwise
def matrix_A(M,m):
    sfs = compute_sfs(M)
    A = np.repeat(sfs, m-1).reshape(m-1, m-1)
    for i in range(m-2):
        www =  weighted_sfs(M, m,i+2)
        A[ :, i] = A[ :, i] - www
    return A

def matrix_B(m): 
    B = np.zeros(m-1)
    B[0] = m
    return B

# Average across independent samples to compute the matrices A and B
def estimate_ris( n, m):
    AAA = np.zeros((m-1, m-1))
    for i in range(n):
        M = SNP_matrix(m, L, alpha)
        AAA += matrix_A(M, m)
    AAA = AAA/n
    BBB = matrix_B(m)
    return AAA, BBB


# Functions to compute the "real" nu measure
def nu_1(x, alpha):
    return (x**(-alpha-1))*((1-x)**(alpha-1))/(sp.gamma(2-alpha)*sp.gamma(alpha))

nu = np.vectorize(nu_1)

def nu_estim_1(x1, m, h, R):
    js = np.linspace(2, m, m-1)/m
    S = st.norm.pdf((x1-js)/h) * R
    return np.sum(S/h)

def nu_estim(X, m ,h, R):
    Y = np.zeros(len(X))
    for i in range(len(X)):
        Y[i] = nu_estim_1(X[i], m, h, R)
    return 


# Compute the real coalescent rates from beta-coalescents (r_mj's)
def coal_rates(b, k, alpha):
   return sp.binom(b, k)*sp.beta(k-alpha, b-k + alpha)/sp.beta(alpha, 2-alpha)


# Optimization with constraints
def solve_constraints(A, B, bnds):
    fun = lambda x: np.linalg.norm(np.dot(A,x)-B)   
    X0 = np.zeros(len(B))
    res = optimize.minimize(fun, X0, method = 'L-BFGS-B', bounds = bnds)
    return res.x   


#PLOT
# When using only the true values of R there is no need to normalize
# normalization (by the integral of Lambda) is needed when the mutation rate is unknown
def plot(h, alpha, m, R, trueR, ndots, tau):
    R = R/sum(R)
    trueR = trueR/sum(trueR)
    x = np.linspace(tau*h+1/m, 1, ndots+1)
    y = nu_estim(x,  m, h, R)
    z = nu_estim(x,  m, h, trueR)
    v =  nu(x, alpha)
    intv = np.sum(v*x*x)/ndots
    v = v/intv
    inty = np.sum(y*x*x)/ndots
    y = y/inty
    intz = np.sum(z*x*x)/ndots
    z = z/intz
    plt.plot(x,v, label = 'True')
    plt.plot(x, y, label = 'Estimated')
    plt.plot(x, z, label = 'Estimated true R')
    plt.legend()
    plt.show()
    #plt.savefig('estimated_R_alpha_' + str(alpha) + '_m_' + str(m) + '_n_' + str(n) + '_L_' + str(L) + '.png')
    #plt.close()
