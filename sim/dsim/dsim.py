'''
dsim.py
one-dimensional, Euler time-step vehicle dynamics simulation.

Matthew Vernacchia
Pyralis
27 July 2013
'''

import numpy as np

class State:
	#### Rocket State ####
	# Altitude. Up is positive. [meter]
	y = 0
	# Vertical velocity [meter sec^-1]
	dy = 0
	# Vertical acceleration [meter sec^-2]
	ddy = 0
	# Mach number [-]
	mach = 0

	#### Simulation State ####
	# Simulation time [sec]
	t = 0

	#### Atmosphere State ####
	# Atmospheric pressure [Pa]
	Pa = 0
	# Atmospheric temperature [K]
	Ta = 0
	# Speed of sound [meter sec^-1]
	a = 0
	def __init__(self):
		pass

class Vehicle:
	# Mass [kg]
	m = 0
	# Coefficient of drag [-]
	Cd = 0
	# Upwards-facing area [meter^2]
	A = 0
	def __init__(self):
		pass

class DSim:
	# simulation state
	state = State()
	# Vehicle parameters
	vehicle = Vehicle()
	# Applied force. Up is positive. [N]
	force = 0
	# Acceleration due to gravity [meter sec^-2]
	g = -9.81
	# Ratio of specific heats for air [-]
	gamma_air = 1.400
	# Gas constant for air [J kg^-1 K^-1]
	R_air = 287.058

	def __int__(self):
		pass

	def atmo(self, y ):
		''' Atmosphere state as a function of altitude.
		@param y Altitude [meter].
		@returns pressure [Pa], temperature [K], speed of sound [meter sec^-1].
		'''
		# Pressure (http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf, eqn 9)
		Pa = 100 * ((44331.514 - y)/(11880.516))**(1/0.1902632)

		# Temperature (1st order lapse rate model, valid in troposphere)
		# lapse rate [K meter^-1]
		lr = -6.5/1000
		# atmo temperature at ground level for Green River, Utah in June (Wolfram|Alpha) [K] 
		Tground = 273.2 + 26
		Ta = Tground + lr*y

		# Speed of sound
		a = ( self.gamma_air * self.R_air * Ta )**0.5

		return Pa, Ta, a

	def update(self, time_new ):
		''' Update the simulation state to a new time.
		@param time_new The new time [sec]. time_new in sequential calls to update must be monotonically increasing.
		For best results the difference between time_new in sequential calls should be small.
		'''
		# update time
		dt = time_new - self.state.t
		assert( dt > 0 )
		self.state.t = time_new

		# Get atmospheric properties
		self.state.Pa, self.state.Ta, self.state.a = self.atmo( self.state.y )

		# Find the drag force on the vehicle
		rho_air = self. state.Pa / ( self.R_air * self. state.Ta ) 
		f_drag = 0.5 * rho_air * (self. state.dy**2) * self.vehicle.A * self.vehicle.Cd * -np.sign(self. state.dy) 

		# Find the net force on the vehicle
		f_net = self.force + f_drag + self.g*self.vehicle.m

		# Find the acceleration
		self. state.ddy = f_net / self.vehicle.m

		# Find the velocity and position
		self. state.dy += self. state.ddy * dt
		self. state.y  += self. state.dy  * dt

		self.state.mach = np.absolute(self.state.dy) / self.state.a
		if self.state.mach > 1:
			print "Warning: supersonic conditions at time %.3f"%(self.state.t)
