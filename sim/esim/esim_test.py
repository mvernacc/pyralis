'''
Engine sim test
'''
import numpy as np
from matplotlib import pyplot as plt
from esim import ESim

esim = ESim()

# load engine params
esim.engine.Isp_opt = 333.1
esim.engine.dry_mass = 3.0
esim.engine.At = 0.00005890
esim.engine.er = 12.0
esim.engine.m_dot_choke = 0.2188
esim.engine.Pc = 4.4e6
esim.engine.Pe = 43379.96

esim.state.m_dot = 0.2188

Pa = np.linspace(101e3, 23e3, 100)
Isp = []
t = 0.1
for p in Pa:
	esim.state.Pa = p
	esim.update(t)
	Isp.append( esim.state.Isp_curr )
	t += 0.1

plt.plot( Pa, Isp )
plt.xlabel('Ambient pressure [Pa]')
plt.ylabel('Isp [sec]')

plt.show()
