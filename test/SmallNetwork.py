from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path
extra_lib_dirs = ['..\\models','..\\examples']
new_model = load_fmu(compile_fmu("SmallNetwork", compiler_options = {'extra_lib_dirs':extra_lib_dirs}))

new_model_result = new_model.simulate(final_time = 20000)


old_model = load_fmu(compile_fmu("oldSmallNetwork", compiler_options = {'extra_lib_dirs':'oldModels'}))

old_model_result = old_model.simulate(final_time = 20000)

samples = min(len(new_model_result['time']),len(old_model_result['time']))

print "Diff pipeline.w_out[1]: ", sum([abs(new_model_result['pipeline.w_out[1]'][i] - old_model_result['pipeline.w_out[1]'][i]) for i in xrange(0,samples)])
print "Diff pipeline.w_out[2]: ", sum([abs(new_model_result['pipeline.w_out[2]'][i] - old_model_result['pipeline.w_out[2]'][i]) for i in xrange(0,samples)])
print "Diff pipeline.w_out[3]: ", sum([abs(new_model_result['pipeline.w_out[3]'][i] - old_model_result['pipeline.w_out[3]'][i]) for i in xrange(0,samples)])

print "Diff separator.w_out[1]: ", sum([abs(new_model_result['separator.outputConnector[1].massFlow[1]'][i] - old_model_result['separator.w_out[1]'][i]) for i in xrange(0,samples)])
print "Diff separator.w_out[2]: ", sum([abs(new_model_result['separator.outputConnector[2].massFlow[2]'][i] - old_model_result['separator.w_out[2]'][i]) for i in xrange(0,samples)])
print "Diff separator.w_out[3]: ", sum([abs(new_model_result['separator.outputConnector[3].massFlow[3]'][i] - old_model_result['separator.w_out[3]'][i]) for i in xrange(0,samples)])

print "Diff nu_u: ",sum([abs(new_model_result['compressor.nu_u'][i] - old_model_result['compressor.nu_u'][i]) for i in xrange(0,samples)])
print "Diff nu_d: ",sum([abs(new_model_result['compressor.nu_d'][i] - old_model_result['compressor.nu_d'][i]) for i in xrange(0,samples)])
print "Diff nu_l: ",sum([abs(new_model_result['compressor.nu_l'][i] - old_model_result['compressor.nu_l'][i]) for i in xrange(0,samples)])
print "Diff nu_r: ",sum([abs(new_model_result['compressor.nu_r'][i] - old_model_result['compressor.nu_r'][i]) for i in xrange(0,samples)])

#print "Diff m_gt: ", sum([abs(new_model_result['well.m_gt'][i] - old_model_result['well.m_gt'][i]) for i in xrange(0,samples)])
#print "Diff m_lt: ", sum([abs(new_model_result['well.m_lt'][i] - old_model_result['well.m_lt'][i]) for i in xrange(0,samples)])