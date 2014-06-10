from pymodelica import compile_fmu
from pyfmi import load_fmu, FMUModel
from os import path
import matplotlib.pyplot as plt

'''
	This example generate the Gas-lift rate vs. Oil production curve and
	the Gas-lift vs. Bottom Hole Pressure. To achieve this 100 well are
	simulated using gas-lift rates from 0.1 to 20.1 kg/s.
'''

extra_lib_dirs = ['..\\models','..\\examples']
new_model = load_fmu(compile_fmu("multiWellTest", compiler_options = {'extra_lib_dirs':extra_lib_dirs}, compiler_log_level = 'e'))

new_model_result = new_model.simulate(final_time = 10000)

w_gl = []
w_oout = []
p_bh = []
for i in range(1,int(new_model_result['n_wells'])+1):
	w_gl.append(new_model_result['well['+`i`+'].w_in'][-1])
	w_oout.append(new_model_result['well['+`i`+'].w_out[2]'][-1])
	p_bh.append(new_model_result['well['+`i`+'].p_bh'][-1])

plt.figure()
plt.plot(w_gl,w_oout)
plt.xlabel('Gas Lift Mass Rate (kg/s)')
plt.ylabel('Oil Outflow Mass Rate (kg/s)')

plt.figure()
plt.plot(w_gl,p_bh)
plt.xlabel('Gas Lift Mass Rate (kg/s)')
plt.ylabel('Bottom Hole Pressure (Pa)')

plt.show()