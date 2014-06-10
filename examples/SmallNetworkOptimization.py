from pymodelica import compile_fmu, compile_fmux
from pyfmi import load_fmu 
from pyjmi import CasadiModel
from pyjmi.common.io import ResultDymolaTextual
import numpy as N
from os import path

'''
	In this example a small network (with 2 wells, 1 subsea production manifold, 1 pipeline, 1 separator, and 1 compressor) 
	is optimized using a infeasible setpoint strategy to maximize the oil production.
	First, the it is run a simulation to obtain the initial conditions using w_gl = 1 kg/s.
	Next, the a first guess using simulation model is generated to initialize the optimization solver (using w_gl = 2 kg/s).
	At last, the optimization problem is solved using the first guess.
'''

extra_lib_dirs = [path.abspath('..\\models'),path.abspath('..\\examples')]

new_model = load_fmu(compile_fmu("SmallNetworkSimulationModel", compiler_options = {'inline_functions':'all','extra_lib_dirs':extra_lib_dirs}))
time = 20000

t = N.linspace(0.,20000.,100) 
u1 = N.ones(100) 
u2 = N.ones(100) 
u_traj = N.transpose(N.vstack((t,u1,u2)))

initial_condition_model_result = new_model.simulate(input = (['w_gl[1]','w_gl[2]'],u_traj),final_time = 20000)

new_model.reset()

states = ['pipeline.m_gp','pipeline.m_op','pipeline.m_wp','pipeline.m_gr','pipeline.m_or','pipeline.m_wr',]
for i in xrange(1,3):
	states.extend(['well['+`i`+'].m_ga','well['+`i`+'].m_gt','well['+`i`+'].m_lt'])
for state in states:
	new_model.set(state+'0',initial_condition_model_result[state][-1])

u_traj = N.transpose(N.vstack((t,2*u1,2*u2)))	
new_model_result = new_model.simulate(input = (['w_gl[1]','w_gl[2]'],u_traj),final_time = 20000)

compiled = compile_fmux("SmallNetworkOptimizationModel", '..\\examples\\SmallNetworkOptimizationModel.mop', compiler_options = {'inline_functions':'all','enforce_bounds':True,'generate_html_diagnostics':False,'extra_lib_dirs':extra_lib_dirs})

new_optimization_model = CasadiModel(compiled)

traj = ResultDymolaTextual(new_model_result.result_file)
opts = new_optimization_model.optimize_options()
opts.update({'init_traj':traj,
			 'nominal_traj': traj, 
			 'n_e': 30,
			 'n_cp': 3,
			 'IPOPT_options': {'tol':1e-4,
							   'dual_inf_tol':1000, 
							   'constr_viol_tol':1000,
							   'compl_inf_tol':100,
							   'max_iter': 1000,
							   'mumps_mem_percent':50
							   }})
result = new_optimization_model.optimize(options = opts)
