TITLE passive membrane with steady voltage-independent synaptic inputs

UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
}

INDEPENDENT { v FROM -100 TO 50 WITH 50 (mV) }

NEURON {
 SUFFIX PasSyn
 NONSPECIFIC_CURRENT ise, isi
 RANGE gse, ese,ise
 RANGE gsi, esi
}

PARAMETER {
	gse  =   0.0	(mho/cm2)	: steady, voltage-independent 
	ese  =   0.0	(mV)		: excitatory (AMPA) synapse

	gsi  =   0.0	(mho/cm2)  : steady, voltage-independent
	esi  = -80		(mV)       : inhibitory (GABA_A) synapse
}

ASSIGNED { i (mA/cm2)
	ise  (mA/cm2)
	isi  (mA/cm2)
}

BREAKPOINT {
	ise = gse * (v - ese )
	isi = gsi * (v - esi )
}


