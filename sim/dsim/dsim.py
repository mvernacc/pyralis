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
	y
	# Vertical velocity [meter sec^-1]
	dy
	# Vertical acceleration [meter sec^-2]
	ddy
	# Mach number [-]
	mach

	#### Simulation State ####
	# Simulation time [sec]
	t

	#### Atmosphere State ####
	# Atmospheric pressure [Pa]
	Pa
	# Atmospheric temperature [K]
	Ta
	# Speed of sound [meter sec^-1]
	a	

class Vehicle:
	# Mass (non-propulsion system) [kg]
	m
	# Coefficient of drag [-]
	Cd
	# Upwards-facing area [meter^2]
	A

class DSim:
	# simulation state
	State state
	# Vehicle parameters
	Vehicle vehicle
	# Applied force. Up is positive. [N]
	force
	# Acceleration due to gravity [meter sec^-2]
	g = -9.81
	# Ratio of specific heats for air [-]
	gamma_air = 1.400
	# Gas constant for air [J kg^-1 K^-1]
	R_air = 287.058

	def atmo( y ):
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
		a = ( DSim.gamma_air * DSim.R_air * Ta )

		return Pa, Ta, a

	def update( time_new ):
		''' Update the simulation state to a new time.
		@param time_new The new time [sec]. time_new in sequential calls to update must be monotonically increasing.
		For best results the difference between time_new in sequential calls should be small.
		'''
		# update time
		dt = time_new - DSim.state.t
		assert( dt > 0 )
		DSim.state.t = time_new

		# Get atmospheric properties
		DSim.state.Pa, DSim.state.Ta, DSim.state.a = atmo( DSim.state.y )

		# Find the drag force on the vehicle
		rho_air = DSim.state.Pa / ( DSim.R_air * DSim.state.Ta ) 
		f_drag = 0.5 * rho_air * (DSim.state.dy**2) * DSim.vehicle.A * DSim.vehicle.Cd * -np.sign(DSim.state.dy) 

		# Find the net force on the vehicle
		f_net = force + f_drag + DSim.g*DSim.vehicle.m

		# Find the acceleration
		DSim.ddy = f_net / DSim.vehicle.m

		# Find the velocity and position
		DSim.dy += DSim.ddy / dt
		DSim.y  += DSim.dy  / dt
