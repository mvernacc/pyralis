'''
Rocket Flight Simulation
'''
import numpy as np
from matplotlib import pyplot as plt
from dsim.dsim import DSim
from esim.esim import ESim

def get_motor_class(I):
	motor_class = ""
	if I > 1280 and I < 2560:
		motor_class = "K"
	if I > 2560 and I < 5120:
		motor_class = "L"
	if I > 5120 and I < 10240:
		motor_class = "M"
	if I > 10240 and I < 20480:
		motor_class = "N"
	return motor_class

dsim = DSim()
# Vehicle at initial rest
dsim.state.y = 0
dsim.state.dy = 0
dsim.state.ddy =0
dsim.state.t = 0

# Vehicle mass and drag
dsim.vehicle.m_struct = 6
dsim.vehicle.A = 50 * 0.01**2
dsim.vehicle.Cd = 0.76

#Enigne
esim = ESim()
esim.engine.Isp_opt = 333.1
esim.engine.dry_mass = 4.0
esim.engine.At = 0.00005890
esim.engine.er = 12.0
esim.engine.m_dot_choke = 0.2188
esim.engine.Pc = 4.4e6
esim.engine.Pe = 43379.96
# Propellant load
esim.state.prop_mass = 1.7

y = []
Ta = []
Pa = []
a = []
mach = []
thrust = []

# turn on engine
esim.state.flame = True
esim.state.m_dot = esim.engine.m_dot_choke

time = np.arange(0.01, 50, 0.01)
for t in time:
	if t > 1:
		esim.state.m_dot = 0.5*esim.engine.m_dot_choke
	dsim.vehicle.m = dsim.vehicle.m_struct + esim.engine.dry_mass + esim.state.prop_mass
	dsim.update(t)
	esim.state.Pa = dsim.state.Pa
	esim.update(t)
	dsim.force = esim.state.f_thrust
	thrust.append( esim.state.f_thrust )
	y.append( dsim.state.y )
	Pa.append( dsim.state.Pa )
	Ta.append( dsim.state.Ta )
	a.append( dsim.state.a )
	mach.append( dsim.state.mach )
time = np.array(time)
y = np.array(y)

print 'Apogee        = %.1f m above sea level'%(np.max(y))
motor_class = get_motor_class(esim.state.impulse)
print 'Total Impulse = %.1f N sec (class %s)'%(esim.state.impulse, motor_class)
plt.subplot(311)
plt.plot( time, y )
plt.xlabel('Time [sec]')
plt.ylabel('Altitude [meter]')

plt.subplot(312)
plt.plot( time, mach, hold=True )
plt.plot( time, np.ones(time.size), 'r' )
plt.xlabel( "Time [sec]")
plt.ylabel('Mach Number [-]')

plt.subplot(313)
plt.plot( time, thrust, hold=True )
plt.xlabel( "Time [sec]")
plt.ylabel('Thrust [N]')
plt.show()