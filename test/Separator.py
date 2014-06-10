from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path
extra_lib_dirs = ['..\\models','..\\examples']
new_model = load_fmu(compile_fmu("SingleSeparator", compiler_options = {'extra_lib_dirs':extra_lib_dirs}))

new_model_result = new_model.simulate(final_time = 100)


old_model = load_fmu(compile_fmu("oldSingleSeparator", compiler_options = {'extra_lib_dirs':'oldModels'}))

old_model_result = old_model.simulate(final_time = 100)

samples = 1

print "Diff w_out[1]: ", (new_model_result['separator.outputConnector[1].massFlow[1]'] - old_model_result['separator.w_out[1]'])
print "Diff w_out[2]: ", (new_model_result['separator.outputConnector[2].massFlow[2]'] - old_model_result['separator.w_out[2]'])
print "Diff w_out[3]: ",(new_model_result['separator.outputConnector[3].massFlow[3]'] - old_model_result['separator.w_out[3]']) 
print "Diff x: ",(new_model_result['separator.x'] - old_model_result['separator.x']) 
print "Diff V_s1: ",(new_model_result['separator.V_s1'] - old_model_result['separator.V_s1']) 
print "Diff V_s2: ",(new_model_result['separator.V_s2'] - old_model_result['separator.V_s2']) 
print "Diff theta_1: ",(new_model_result['separator.theta_1'] - old_model_result['separator.theta_1']) 
print "Diff L_1: ",(new_model_result['separator.L_1'] - old_model_result['separator.L_1']) 
print "Diff h_1: ",(new_model_result['separator.h_1'] - old_model_result['separator.h_1']) 
print "Diff v_h: ",(new_model_result['separator.v_h'] - old_model_result['separator.v_h']) 
print "Diff v_v: ",(new_model_result['separator.v_v'] - old_model_result['separator.v_v']) 

