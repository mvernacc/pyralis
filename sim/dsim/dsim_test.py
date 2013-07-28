'''
dsim test
'''
import numpy as np
from matplotlib import pyplot as plt
from dsim import DSim

dsim = DSim()
dsim.state.y = 0
dsim.state.dy = 0
dsim.state.ddy =0
dsim.state.t = 0

dsim.vehicle.m = 10
dsim.vehicle.A = 50 * 0.01**2
dsim.vehicle.Cd = 0.76

dsim.force = 900

y = []
Ta = []
Pa = []
a = []
mach = []

time = np.arange(0.01, 50, 0.01)
for t in time:
	if t > 2.5:
		dsim.force = 750
	if t > 5:
		dsim.force = 700
	if t > 5.5:
		dsim.force = 500
	if t > 6:
		dsim.force = 200
	if t > 7:
		dsim.force = 80
	if t > 8:
		dsim.force = 0
	dsim.update(t)
	y.append( dsim.state.y )
	Pa.append( dsim.state.Pa )
	Ta.append( dsim.state.Ta )
	a.append( dsim.state.a )
	mach.append( dsim.state.mach )
time = np.array(time)
y = np.array(y)

plt.subplot(211)
plt.plot( time, y )
plt.xlabel('Time [sec]')
plt.ylabel('Altitude [meter]')

plt.subplot(212)
plt.plot( time, mach, hold=True )
plt.plot( time, np.ones(time.size), 'r' )
plt.xlabel( "Time [sec]")
plt.ylabel('Mach Number [-]')
plt.show()