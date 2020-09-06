TITLE Low threshold calcium current
:
:   Ca++ current responsible for low threshold spikes (LTS)
:   RETICULAR THALAMUS
:   Differential equations
:
:   Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.
:   
:   Written by Alain Destexhe, Salk Institute, Sept 18, 1992
:   
:    - Biophysical properties of the T current were from recordings of
:    - human recombinant Cav3.2 T-channel in HEK-293 cells
:    - see Vitko et al., J. Neurosci 25(19) :4844-4855, 2005
:    - Q10 and shift parameters are fixed 
:   
:
:   

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX it2
	USEION ca READ cai, cao WRITE ica
	RANGE gcabar, ica
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)
	celsius	= 23	(degC)
:	eca	= 120	(mV)
	gcabar	= 1e-10	(mho/cm2)
	shift	= 28 	(mV)
	cai	= 1e-4 (mM)		: adjusted for eca=120 mV
	cao	= 2.5	(mM)
	a01=29	:52 original
	a02=7	:7.4
	a03=50	:80
	a04=5
	a05=23		:27
	a06=5.67	:10
	a07=100		:102
	a08=15.33	:15
	a09=13		:48
	a10=4.32	:4
	a11=395.33	:407
	a12=53.33	:50
	a13=1.44	:1
	a14=26.67	:28.3
}

STATE {
	m h
}

ASSIGNED {
	ica	(mA/cm2)
	carev	(mV)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar * m*m*h * (v-carev)
}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {
:
:   Activation functions and kinetics were obtained from
:   Vitko et al., 2005 at 23-25 deg.
:
:   Activation functions and kinetics were obtained from
:   Huguenard & Prince, and were at 23-25 deg.

	evaluate_fct(v)

	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 
:
:   Time constants were obtained from Vitko, Arias, Perez-Reyes
:
LOCAL q10
	q10 = 3^((celsius - 23)/10)
	m_inf = 1.0 / ( 1 + exp(-(v+shift+a01)/a02) )
	h_inf = 1.0 / ( 1 + exp((v+shift+a03)/a04) )

	tau_m = a13 + (1.0 /(q10* ( exp((v+shift+a05)/a06) + exp(-(v+shift+a07)/a08))))
	tau_h = a14 + (1.0 /(q10*( exp((v+shift+a09)/a10) + exp(-(v+shift+a11)/a12))))
}
UNITSON
