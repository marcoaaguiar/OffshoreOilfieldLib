from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path
extra_lib_dirs = ['..\\models','..\\examples']
new_model = load_fmu(compile_fmu("SinglePipelineRiser", compiler_options = {'extra_lib_dirs':extra_lib_dirs}))

new_model_result = new_model.simulate(final_time = 100)


old_model = load_fmu(compile_fmu("oldSinglePipelineRiser", compiler_options = {'extra_lib_dirs':'oldModels'}))

old_model_result = old_model.simulate(final_time = 100)

samples = min(len(new_model_result['time']),len(old_model_result['time']))

print "Diff m_gp: ", sum([abs(new_model_result['pipeline.m_gp'][i] - old_model_result['pipeline.m_gp'][i]) for i in xrange(0,samples)])/samples
print "Diff m_op: ", sum([abs(new_model_result['pipeline.m_op'][i] - old_model_result['pipeline.m_op'][i]) for i in xrange(0,samples)])/samples
print "Diff m_wp: ", sum([abs(new_model_result['pipeline.m_wp'][i] - old_model_result['pipeline.m_wp'][i]) for i in xrange(0,samples)])/samples

print "Diff m_gr: ", sum([abs(new_model_result['pipeline.m_gr'][i] - old_model_result['pipeline.m_gr'][i]) for i in xrange(0,samples)])/samples
print "Diff m_or: ", sum([abs(new_model_result['pipeline.m_or'][i] - old_model_result['pipeline.m_or'][i]) for i in xrange(0,samples)])/samples
print "Diff m_wr: ", sum([abs(new_model_result['pipeline.m_wr'][i] - old_model_result['pipeline.m_wr'][i]) for i in xrange(0,samples)])/samples

print "Diff w_out[1]: ", sum([abs(new_model_result['pipeline.w_out[1]'][i] - old_model_result['pipeline.w_out[1]'][i]) for i in xrange(0,samples)])/samples
print "Diff w_out[2]: ", sum([abs(new_model_result['pipeline.w_out[2]'][i] - old_model_result['pipeline.w_out[2]'][i]) for i in xrange(0,samples)])/samples
print "Diff w_out[3]: ", sum([abs(new_model_result['pipeline.w_out[3]'][i] - old_model_result['pipeline.w_out[3]'][i]) for i in xrange(0,samples)])/samples