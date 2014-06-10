from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path

'''
	This example simulates a network with 2 wells, 1 subsea production manifold, 1 pipeline, 1 separator, and 1 compressor.
'''

extra_lib_dirs = ['..\\models','..\\examples']
new_model = load_fmu(compile_fmu("SmallNetwork", compiler_options = {'extra_lib_dirs':extra_lib_dirs}))

new_model_result = new_model.simulate(final_time = 20000)