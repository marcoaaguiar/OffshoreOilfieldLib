from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path

'''
	This example simulates the well dynamic response.
'''

extra_lib_dirs = ['..\\models','..\\examples','..\\newModelTest']
new_model = load_fmu(compile_fmu("newWellTest", compiler_options = {'extra_lib_dirs':extra_lib_dirs}))

opts = new_model.simulate_options()
opts['CVode_options']['rtol'] = 1e-6
new_model_result = new_model.simulate(final_time = 20000, options = opts)