import numpy
import scipy
import matplotlib.pyplot as plt
from collections import defaultdict

#Simulation Function
def wf_sim(n, p_init, n_sim, t):
    simulations = []
    times_to_fixation = []
    for s in range(n_sim):
        p = p_init
        total_n_A = []
        fixated = False
        for i in range(t):
            curr_n_A = numpy.random.binomial(n, p)
            total_n_A.append(curr_n_A)
            p = curr_n_A/n
            if (p == 1 or p == 0) & (fixated == False): 
                times_to_fixation.append(i)
                fixated = True
        simulations.append(total_n_A)
    return(simulations, times_to_fixation)

#Entropy function (for allele A and allele B)
def S(n_A, N):
    n_B = N - n_A
    p_A = (n_A)/N
    p_B = (n_B)/N
    return -(p_A*numpy.log(p_A) + p_B*numpy.log(p_B))

#Time to fixation function (for 2 alleles)
def T(n_pop, n_A):
    return 2*n_pop*S(n_A, n_pop)

#Calculate expectations and variances
def wf_sim_expvar(curr_n, list_p, n_sim, t):
    expectations = []
    variances = []
    for curr_p in list_p:
        simulations, times_to_fixation = wf_sim(curr_n, curr_p, n_sim, t)
        expectations.append(numpy.mean(times_to_fixation))
        variances.append(numpy.var(times_to_fixation))
    return expectations, variances
