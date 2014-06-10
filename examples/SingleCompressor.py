from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path

'''
	This example simulates the compressor static response.
'''

extra_lib_dirs = ['..\\models','..\\examples']
new_model = load_fmu(compile_fmu("SingleCompressor", compiler_options = {'extra_lib_dirs':extra_lib_dirs}))

new_model_result = new_model.simulate(final_time = 100)