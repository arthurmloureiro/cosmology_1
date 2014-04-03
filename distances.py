"""
	PROGRAMA PARA CALCULAR DISTANCIAS EM COSMOLOGIA
		Arthur Eduardo da Mote Loureiro
			03/04/2014
"""
import numpy as np
import pylab as pl
import scipy.integrate as inte

h=0.72
#	Fiducial values		#
omf =0.25					# Matter
odef = 0.75					# Dark Energy
orf = 8.2E-5					# Radiation
wf = -1.					# Eq. of state
okf = 1. - omf - odef - orf			# Curvature
H0 = 1.023E-1*h					# Hubble in Gyrs

def H(z,om,ode,orr,ok,w):
	aom = om*np.power((1+z), 3)
	aok = ok*np.power((1+z), 2)
	aor = orr*np.power((1+z), 4)
	aode = ode*np.power((1+z),(3+3*w))
	return 1./(H0*np.sqrt(aom+aok+aor+aode))

def Dx(z,om,ode,orr,ok,w):
	return inte.quad(H,0,z,args=(om,ode,orr,ok,w))[0]

def Da(z,om,ode,orr,ok,w):
	if ok == 0.0:
		return Dx(z,om,ode,orr,ok,w)
	elif ok < 0.0:
		K = H0*H0*np.abs(ok)
		return np.power((K), -1./2)*np.sin(np.sqrt(K)*Dx(z,om,ode,orr,ok,w))
	else:
		K = H0*H0*np.abs(ok)
		return np.power((K), -1./2)*np.sinh(np.sqrt(K)*Dx(z,om,ode,orr,ok,w))

def Dl(z,om,ode,orr,ok,w):
	return Da(z,om,ode,orr,ok,w)*(1+z)
z1 = np.linspace(0.0,10.,1000)
pl.figure()
pl.plot(integral(z1,omf,odef,orf,okf,wf))
pl.show()
















