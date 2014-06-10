from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path
extra_lib_dirs = ['..\\models','..\\examples']
new_model = load_fmu(compile_fmu("SingleCompressor", compiler_options = {'extra_lib_dirs':extra_lib_dirs}))

new_model_result = new_model.simulate(final_time = 100)


old_model = load_fmu(compile_fmu("oldSingleCompressor", compiler_options = {'extra_lib_dirs':'oldModels'}))

old_model_result = old_model.simulate(final_time = 100)

samples = 1

print "Diff w_out[1]: ", (new_model_result['compressor.outputConnector.massFlow[1]'] - old_model_result['compressor.w_out'])

print "Diff nu_u: ",(new_model_result['compressor.nu_u'] - old_model_result['compressor.nu_u']) 
print "Diff nu_d: ",(new_model_result['compressor.nu_d'] - old_model_result['compressor.nu_d']) 
print "Diff nu_l: ",(new_model_result['compressor.nu_l'] - old_model_result['compressor.nu_l']) 
print "Diff nu_r: ",(new_model_result['compressor.nu_r'] - old_model_result['compressor.nu_r']) 


