'''
esim.py
Engine thrust and weight simulation
Matthew Vernacchia
Pylaris
27 July 2013
'''

class State:
	# Simulation time [sec]
	t = 0
	# Mass of propellant in tanks [kg]. 
	prop_mass = 0
	# Ambient Pressure [Pa].
	Pa = 0
	# Mass flow rate [kg sec^-1].
	m_dot = 0
	# Specific impulse (current) [sec].
	Isp_curr = 0
	# Thrust force produced [N].
	f_thrust = 0
	def __init__(self):
		pass

class Engine:
	# Specific impulse (optimal) [sec].
	Isp_opt = 0
	# Dry mass [kg].
	dry_mass = 0
	# Throat area [meter^2].
	At = 0
	# Expansion ratio [-].
	er = 0
	# Max. mass flow rate [kg sec^-1]
	m_dot_choke = 0
	# Chamber pressure [Pa].
	Pc = 0
	# Exit pressue [Pa].
	Pe = 0
	def __init__(self):
		pass

class ESim:
	state = State()
	engine = Engine()
	# Acceleration due to gravity [meter sec^-2]
	g = 9.81

	def update(self, time_new):
		# update time
		dt = time_new - self.state.t
		assert( dt > 0 )
		self.state.t = time_new

		# check for flow choking
		if self.state.m_dot > self.engine.m_dot_choke:
			print 'ESim Warning: choked flow at time %f'%(self.state.t)
			self.state.m_dot = self.engine.m_dot_choke

		# Find c* (source: RPE eqn 2-18)
		c_star = self.engine.Pc * self.engine.At / self.state.m_dot

		# Find the altitude-compensated specific impulse (source: RPE eqn 3-28)
		#TODO May not be valid for aerospike
		self.state.Isp_curr = self.engine.Isp_opt + \
			c_star * self.engine.er / (self.g * self.engine.Pc) * (self.engine.Pe - self.state.Pa)

		# Find the thrust
		self.state.f_thrust = self.state.Isp_curr * self.g * self.state.m_dot

		# Account for propellant mass change
		self.state.prop_mass -= self.state.m_dot * dt