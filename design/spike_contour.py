from math import atan, sin, tan, pi, asin
import numpy as np
from matplotlib import pyplot as plt

def get_Mach( er, y ):
    '''
    Explicit Inversion of Stodola's Area-Mach Equation
    Source: J. Majdalani and B. A. Maickie
    http://maji.utsi.edu/publications/pdf/HT02_11.pdf
    '''
    # order of the aproximation
    n = 5
    X = np.zeros((n,))
    M = np.zeros((n,))

    e = 1/er
    B = (y+1)/(y-1)
    k = np.sqrt( 0.5*(y-1) )
    u = e**(1/B) / np.sqrt( 1+k**2 )
    X[0] = (u*k)**(B/(1-B))
    M[0] = X[0]

    for i in xrange(1,n):
        lamb = 1/( 2*M[i-1]**(2/B)*(B-2) + M[i-1]**2 *B**2*k**2*u**2 )
        X[i] = lamb*M[i-1]*B*( M[i-1]**(2/B) - M[i-1]**2*B*k**2*u**2 \
            + ( M[i-1]**(2+2/B)*k**2*u**2*(B**2-4*B+4) \
            - M[i-1]**2*B**2*k**2*u**4 + M[i-1]**(4/B)*(2*B-3) \
            + 2*M[i-1]**(2/B)*u**2*(2-B) )**0.5 )
        M[i] = M[i-1] + X[i]
    if abs( np.imag( M[n-1] ) ) > 1e-5:
        print 'Warning: Exit Mach Number has nonzero imaginary part'
    Me = np.real( M[n-1] )
    return Me

#### Combustion Chamber Parameters ####
# Combustion chamber output stagnation pressure [Pa]
Pc = 4.4e6
# Combustion chamber output stagnation temperature [K]
# 3200 K = estimate of stoich ethane/N2O
# 1727 K = CEARUN ethane:N20 = 3:1 by mass.
Tc = 1727.0

#### Operating Conditions ####
# Ambient atmospheric pressure [Pa]
# atmo pressure at 11000 m elevation
Pa = 23.0e3

#### Algorithm Inputs ####
# Desired Expansion Ratio []
er = 12.0
# Number of contour points
N = 100
# Gas Temperature at shroud lip [K] - approx equal to comb chamber Tstag
Tl = Tc
# Acceleration due to gravity at surface [m sec^-2]
g = 9.81
# Gas Ratio of specific heats (gamma)
y = 1.2
# Gas constant [J kg^-1 K^-1]
# avg combustion product molar mass
molar_m = (2*(12.01+32) + 3*(18.02) + 7*(14.01*2))/(2+3+7)
R_mol = 8.314 # [J mol^-1 K^-1]
R =  R_mol / molar_m * 1000

#### Compute the Exit Mach Number ####
Me = get_Mach( er, y )
# use this to find the exit (tip-plane) pressure  (assumung isetropic expansion)
Pe = Pc * (1 + (y-1)/2*Me**2)**(-y/(y-1))

#### Find the total flow turning angle [rad] ####
nu_e = ((y+1)/(y-1))**0.5 * atan( ((y-1)/(y+1)*(Me**2-1))**0.5 ) \
    - atan( (Me**2-1)**0.5 )
# use this to find the angle of the shroud edge with the normal to
# the axial direction [rad]
delta = pi/2 - nu_e

#### Find the throat gap width / outer radius ratio [] ####
htRe = (er - (er*(er-sin(delta)))**0.5) / (er*sin(delta))

#### Find the velocity and temperature at the throat (assuming M=1)
# throat thermo temperature [K]
Tt = Tc / (1 + (y-1)/2)
# throat pressure [Pa]
Pt = Pc * (1 + (y-1)/2) ** (-y/(y-1))
# throat fluid velocity [meter sec^-1]
vt = np.sqrt( y*R*Tt )

#### Compute the optimum thrust coefficient [] ####
CF_opt = y*Me*(2/(y+1))**((y+1)/(2*y-2)) * (1+(y-1)/2*Me)**(-0.5)

#### Examine a range of Mach numbers to determine the shape of the nozzle
# the mach numbers to examine
M = np.linspace(1, Me, N)
# the radius of the plug at some point x, normalized by the outer radius
RxRe = np.zeros((N,))
# The axial distance from the lip to point x, normalized by the outer
# radius
X = np.zeros((N,))
# The flow velocity to mach wave angle [rad]
mu = np.zeros((N,))
# the turning angle [rad]
nu = np.zeros((N,))
# the pressure at point x [Pa]
P = np.zeros((N,))
# the temperature at point x [K]
T = np.zeros((N,))
# the cumulative Isp up to point x [sec]
Isp = np.zeros((N,))
# the Isp due to momentum flux and pressure at the throat [sec]
Isp[0] = vt*sin(delta)*( 1 + 1/y*(1-((y+1)/2)**(y/(y-1))*Pa/Pt) ) / g
for x in xrange(N):
    mu[x] = asin( 1/M[x] ) # See Lozano's Fall2012 Lec17 Notes
    # use the Prandtl-Meyer equation to find nu[x]
    nu[x] = ((y+1)/(y-1))**0.5 * atan( ((y-1)/(y+1)*(M[x]**2-1))**0.5 ) \
    - atan( (M[x]**2-1)**0.5 )
    # use CC Lee Eqn (26) to find Rx/Re for the point
    RxRe[x] = ( 1 - ( 2/(y+1)*(1+(y-1)/2*M[x]**2))**((y+1)/(2*y-2)) \
        * sin( nu_e - nu[x] + mu[x] ) / er ) ** 0.5
    # find the X (axial) coordinate of the point. CC Lee Eqn (19)
    X[x] = (1 - RxRe[x]) / tan( nu_e - nu[x] + mu[x] )
    # find the pressure
    P[x] = Pc * (1 + (y-1)/2*M[x]**2)**(-y/(y-1))
    if x > 0:
        Isp[x] = Isp[x-1] + ( vt/y*((y+1)/2)**(y/(y-1)) * \
            er/2 * ( (P[x-1]-Pa)/Pt + (P[x]-Pa)/Pt ) * \
            (RxRe[x-1]**2 - RxRe[x]**2) ) / g
    #find the temperature
    T[x] = Tc / (1+ (y-1)/2 * M[x])

plt.suptitle('Simulated Nozzle Conditions vs Axial Distance from Throat Normalized by Exit Radius')
plt.subplot(2,3,1)
plt.plot( X, RxRe )
plt.ylabel('Rx / Re')
plt.grid(True)

plt.subplot(2,3,4)
plt.plot( X, M )
plt.ylabel( 'Mach Number')
plt.xlabel('Xx / Re')
plt.grid(True)

plt.subplot(2,3,2)
plt.plot( X, P/1.0e6, hold=True )
plt.plot( X, np.ones((N,))*Pa/1.0e6, 'r' )
plt.ylabel('Pressure [MPa]')
plt.grid(True)

plt.subplot(2,3,5)
plt.plot( X, T )
plt.ylabel('Temperature [K]')
plt.xlabel('Xx / Re')
plt.grid(True)

plt.subplot(2,3,3)
plt.plot( X, Isp )
plt.ylabel(r'Cumulative $I_{sp}$')
plt.xlabel('Xx / Re')
plt.grid(True)

plt.subplot(2,3,6)
plt.text(0.05,0.8, 'Chamber Conditions:')
plt.text(0.1,0.7, r'$P_c =$ %.2f MPa'%(Pc/1e6))
plt.text(0.1,0.6, r'$T_c =$ %.0f K'%(Tc))
plt.text(0.1,0.5, r'$M_{molar} =$ %.1f g mol^-1'%(molar_m))
#plt.axes([0,10,0,10])
#### Introduce the dimensional parameter Re to find the thrust
# The radius of the shroud lip [meter]
Re = 0.015
# throat width [meter]
ht = htRe * Re
# throat area [meter^2]
At = pi*ht*(2*Re - ht*sin(delta))
# fluid density at the throat [kg m^-3]
rho_t = Pt / (R*Tt)
# throat mass flow [kg sec^-1]
m_dot = rho_t*At*vt
# Find the thrust [N]
F = Isp[N-1] * m_dot * g

#### Print output ####
print 'Engine Geometry:'
print '\tShroud lip radius,   Re = %.1f mm'%(Re*1000)
print '\tThroat width,        ht = %.2f mm'%(ht*1000)
print '\tThroat area,         At = %.8f m^2'%(At)
print '\tExpansion ratio,     er = %.2f'%(er)
print 'Chamber Conditions:'
print '\tChmaber pressure,    Pc = %.3f MPa'%(Pc/1.0e6)
print '\tChamber temperature, Tc = %.0f K'%(Tc)
print '\tExhaust Avg Molar Mass  = %.1f g mol^-1'%(molar_m)
print 'Engine Performance:'
print '\tMass flow rate,   m_dot = %.4f kg sec^-1'%(m_dot)
print '\tThrust force,         F = %.1f N'%(F)
print '\tSpecific impulse,   Isp = %.1f sec'%(Isp[N-1])
print '\tExit pressure,       Pe = %.2f Pa'%(Pe)
print '\tExit Mach number,    Me = %.2f'%(Me)

#### write the contour to a file
f = open('nozzle_plug_curve.txt', 'w')
for x in xrange(N):
    f.write('%.4f,%.4f,%.4f\r\n'%(0, X[x]*Re*-1000, RxRe[x]*Re*1000))
f.close()
f = open('nozzle_shroud_curve.txt', 'w')
for r in np.arange(0, 0.003, 0.0001):
    f.write('%.4f,%.4f,%.4f\r\n'%( 0, r*tan(delta)*1000, (Re+r)*1000))
f.close()

plt.show()