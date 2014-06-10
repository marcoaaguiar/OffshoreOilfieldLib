from pymodelica import compile_fmu, compile_fmux
from pyfmi import load_fmu 
from pyjmi import CasadiModel
from pyjmi.common.io import ResultDymolaTextual
import numpy as N
from os import path

extra_lib_dirs = [path.abspath('..\\models'),path.abspath('..\\examples'),path.abspath('..\\test\\oldModels')]
#extra_lib_dirs = ['C:\\Users\\Aguiar\\Documents\\GitHub\\OffshoreOilfieldLib\\models','C:\\Users\\Aguiar\\Documents\\GitHub\\OffshoreOilfieldLib\\examples']
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
new_model_first_guess = new_model.simulate(input = (['w_gl[1]','w_gl[2]'],u_traj),final_time = 20000)

###

compiled = compile_fmux("SmallNetworkOptimizationModel", 'C:\\Users\\Aguiar\\Documents\\GitHub\\OffshoreOilfieldLib\\examples\\SmallNetworkOptimizationModel.mop', compiler_options = {'inline_functions':'all','enforce_bounds':True,'generate_html_diagnostics':False,'extra_lib_dirs':extra_lib_dirs})

new_optimization_model = CasadiModel(compiled)

traj = ResultDymolaTextual(new_model_first_guess.result_file)
opts = new_optimization_model.optimize_options()
opts.update({'init_traj':traj,
			 'nominal_traj': traj, 
			 'n_e': 30,
			 'n_cp': 2,
			 'IPOPT_options': {'tol':1e-4,
							   'dual_inf_tol':1000, 
							   'constr_viol_tol':1e-6,
							   'compl_inf_tol':100,
							   'max_iter': 1000,
							   'mumps_mem_percent':50
							   }})
new_model_result = new_optimization_model.optimize(options = opts)

########################################################33
##########################################################

old_model = load_fmu(compile_fmu("oldSmallNetworkSimulationModel", compiler_options = {'inline_functions':'all','extra_lib_dirs':extra_lib_dirs}))

old_model_first_guess = old_model.simulate(input = (['w_gl[1]','w_gl[2]'],u_traj),final_time = 20000)

####

compiled = compile_fmux("oldSmallNetworkOptimizationModel", 'C:\\Users\\Aguiar\\Documents\\GitHub\\OffshoreOilfieldLib\\test\\oldModels\\oldSmallNetworkOptimizationModel.mop', compiler_options = {'inline_functions':'all','enforce_bounds':True,'generate_html_diagnostics':False,'extra_lib_dirs':extra_lib_dirs})

old_optimization_model = CasadiModel(compiled)

traj = ResultDymolaTextual(old_model_first_guess.result_file)
opts = old_optimization_model.optimize_options()
opts.update({'init_traj':traj,
			 'nominal_traj': traj, 
			 'n_e': 30,
			 'n_cp': 2,
			 'IPOPT_options': {'tol':1e-4,
							   'dual_inf_tol':1000, 
							   'constr_viol_tol':1e-6,
							   'compl_inf_tol':100,
							   'max_iter': 1000,
							   'mumps_mem_percent':50
							   }})
old_model_result = old_optimization_model.optimize(options = opts)

samples = min(len(new_model_result['time']),len(old_model_result['time']))

print "Diff pipeline.w_out[1]: ", sum([abs(new_model_result['pipeline.w_out[1]'][i] - old_model_result['pipeline.w_out[1]'][i]) for i in xrange(0,samples)])
print "Diff pipeline.w_out[2]: ", sum([abs(new_model_result['pipeline.w_out[2]'][i] - old_model_result['pipeline.w_out[2]'][i]) for i in xrange(0,samples)])
print "Diff pipeline.w_out[3]: ", sum([abs(new_model_result['pipeline.w_out[3]'][i] - old_model_result['pipeline.w_out[3]'][i]) for i in xrange(0,samples)])
















