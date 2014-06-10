from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path
extra_lib_dirs = ['..\\models','..\\examples']
new_model = load_fmu(compile_fmu("FourWellOnePipelineNetwork", compiler_options = {'extra_lib_dirs':extra_lib_dirs}))

new_model_result = new_model.simulate(final_time = 5000)


old_model = load_fmu(compile_fmu("oldFourWellOnePipelineNetwork", compiler_options = {'extra_lib_dirs':'oldModels'}))

old_model_result = old_model.simulate(final_time = 5000)

samples = min(len(new_model_result['time']),len(old_model_result['time']))

print "Diff pipeline.w_out[1]: ", sum([abs(new_model_result['pipeline.w_out[1]'][i] - old_model_result['pipeline.w_out[1]'][i]) for i in xrange(0,samples)])
print "Diff pipeline.w_out[2]: ", sum([abs(new_model_result['pipeline.w_out[2]'][i] - old_model_result['pipeline.w_out[2]'][i]) for i in xrange(0,samples)])
print "Diff pipeline.w_out[3]: ", sum([abs(new_model_result['pipeline.w_out[3]'][i] - old_model_result['pipeline.w_out[3]'][i]) for i in xrange(0,samples)])

#print "Diff m_gt: ", sum([abs(new_model_result['well.m_gt'][i] - old_model_result['well.m_gt'][i]) for i in xrange(0,samples)])
#print "Diff m_lt: ", sum([abs(new_model_result['well.m_lt'][i] - old_model_result['well.m_lt'][i]) for i in xrange(0,samples)])