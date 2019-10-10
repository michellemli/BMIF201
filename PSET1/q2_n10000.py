import sim_fn

#Simulation Parameters
p = numpy.arange(0.01, 0.99, 0.04)
n_sim = 2000
t = 1000

#Run N = 10000 simulation
expectations_N10000, variances_N10000 = wf_sim_expvar(curr_n = 10000, list_p = p, n_sim = n_sim, t = t)

print(expectations_N10000)

print(variances_N10000)