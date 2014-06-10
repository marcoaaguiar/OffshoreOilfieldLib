from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path
extra_lib_dirs = ['..\\models','..\\examples']
new_model = load_fmu(compile_fmu("SingleWell", compiler_options = {'extra_lib_dirs':extra_lib_dirs}))

new_model_result = new_model.simulate(final_time = 100)


old_model = load_fmu(compile_fmu("oldSingleWellNetwork", compiler_options = {'extra_lib_dirs':'oldModels'}))

old_model_result = old_model.simulate(final_time = 100)

samples = min(len(new_model_result['time']),len(old_model_result['time']))

print "Diff m_ga: ", sum([abs(new_model_result['well.m_ga'][i] - old_model_result['well.m_ga'][i]) for i in xrange(0,samples)])
print "Diff m_gt: ", sum([abs(new_model_result['well.m_gt'][i] - old_model_result['well.m_gt'][i]) for i in xrange(0,samples)])
print "Diff m_lt: ", sum([abs(new_model_result['well.m_lt'][i] - old_model_result['well.m_lt'][i]) for i in xrange(0,samples)])