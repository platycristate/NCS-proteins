NEURON {
	SUFFIX hpca
	USEION ca READ cao, ica, cai  WRITE cai, ica
	USEION k READ ek WRITE ik
	RANGE HPCA_m_z, HPCA_mut_m_z, ica_pmp, gbar, ik
	GLOBAL Volume,Ra, Rb, caix,q10, temp, tadj, vmin, vmax
	THREADSAFE
}

DEFINE N 4

UNITS {
	(molar) = (1/liter)
	(mol) = (1)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	FARADAY = (faraday) (10000 coulomb)
	PI = (pi) (1)
	(pS) = (picosiemens)
}

PARAMETER {
	DCa = .22 (um2/ms)
	Dhpca = .05(um2/ms)
	Dcahpca = .05 (um2/ms)
	Dca2hpca = .05 (um2/ms)
	Dhpca_m = .001 (um2/ms)
	k1HPCA = 40 (/mM-ms)
	k2HPCA = .01 (/ms)
	k3HPCA = 40 (/mM-ms)
	k4HPCA = .01 (/ms)
	k7HPCA = .01 (/mM-ms)
	k8HPCA = .002 (/ms)
	TotalHPCA = .03821 (mM)
	
	Dhpca_mut = .05(um2/ms)
	Dcahpca_mut = .05 (um2/ms)
	Dca2hpca_mut = .05 (um2/ms)
	Dhpca_mut_m = .001(um2/ms)
	k1HPCA_mut = 12 (/mM-ms) 					
	k2HPCA_mut = .01 (/ms)
	k3HPCA_mut = 9 (/mM-ms)
	k4HPCA_mut = .01 (/ms)
	k7HPCA_mut = .07 (/mM-ms)
	k8HPCA_mut = .003 (/ms)
	
	D_Bufer = .05(um2/ms)
	Bufer0 = 20 (mM)
	k1bufer = 10 (/mM-ms)
	k2bufer = 1 (/ms)
	
	k1Pump = 1 (/mM-ms)
	k2Pump = .0003 (/ms)
	TotalPump = 1e-11 (mol/cm2)
	cai0 = 110e-6 (mM)
	delta = 0.1 (um)
	
	gbar = 600   	(pS/um2)
	caix = 1
	Ra   = .01	(/ms)
	Rb   = .02	(/ms)
	temp = 36	(degC)
	q10  = 2.3
	vmin = -120	(mV)
	vmax = 100	(mV)
}

ASSIGNED {
	diam	(um)
	ica		(mA/cm2)
	cai		(mM)
	Volume[N] (um2)
	cao (mM)
	ica_pmp (mA/cm2)
	parea (um)
	a		(/ms)
	b		(/ms)
	ik 		(mA/cm2)
	gk		(pS/um2)
	ek		(mV)
	ninf
	ntau 		(ms)	
	tadj
}

CONSTANT { volo = 1e10 (um2) }

STATE {
	ca [N]			(mM)	<1e-10>
	HPCA [N] 		(mM)
	CaHPCA [N] 		(mM)
	Ca2HPCA[N] 		(mM)
	HPCA_m[N] 		(mM)
	
	HPCA_mut [N] 	(mM)
	CaHPCA_mut [N] 	(mM)
	Ca2HPCA_mut [N] (mM)
	HPCA_mut_m[N] 	(mM)
	
	pump (mol/cm2)
	pumpca (mol/cm2)
	n
	
	HPCA_z (mM)
	CaHPCA_z (mM)
	Ca2HPCA_z (mM)
	HPCA_m_z (mM)
	HPCA_tot_z (mM)
	
	HPCA_mut_z (mM)
	CaHPCA_mut_z (mM)
	Ca2HPCA_mut_z (mM)
	HPCA_mut_m_z (mM)
	HPCA_mut_tot_z (mM)
	
	Bufer[N] (mM)
	CaBufer[N] (mM)
}

BREAKPOINT { 
	SOLVE state METHOD sparse
	ica = ica_pmp
	SOLVE states METHOD cnexp
	gk = tadj*gbar*n
	ik = (1e-4) * gk * (v - ek)
}

LOCAL factors_done, nexp

DERIVATIVE states {
	rates(cai)
	n' =  (ninf-n)/ntau
}

PROCEDURE rates(cai(mM)) {  
	a = Ra * (HPCA_m[0]+HPCA_mut_m[0])^caix
	b = Rb
	tadj = q10^((celsius - temp)/10)
	ntau = 1/tadj/(a+b)
	ninf = a/(a+b)
}

INITIAL {
	MUTEXLOCK
	if (factors_done == 0) {
		factors_done = 1
		factors()
	}
	MUTEXUNLOCK
	cai = cai0
	FROM i=0 TO N-1 {
		ca[i] = cai
		HPCA[i] = TotalHPCA
		HPCA_mut[i] = TotalHPCA
		Bufer[i] = Bufer0
	}
	parea = PI*diam
	pump = TotalPump/(1 + (cai*k1Pump/k2Pump))
	pumpca = TotalPump - pump
	rates(cai)
	n = ninf
}

LOCAL frat[N]

PROCEDURE factors() {
	LOCAL r, dr2
	r = 1/2
	dr2 = r/(N-1)/2
	Volume[0] = 0
	frat[0] = 2*r
	FROM i=0 TO N-2 {
		Volume[i] = Volume[i] + PI*(r-dr2/2)*2*dr2
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)
		r = r - dr2
		Volume[i+1] = PI*(r+dr2/2)*2*dr2
	}
}

LOCAL dsq, dsqvol

KINETIC state {
	COMPARTMENT (1e10)*parea {pump pumpca}
	COMPARTMENT volo {cao}
	COMPARTMENT i, Volume [i] {ca HPCA CaHPCA Ca2HPCA HPCA_m HPCA_mut CaHPCA_mut Ca2HPCA_mut HPCA_mut_m Bufer}
	
	LONGITUDINAL_DIFFUSION i, DCa*Volume [i] {ca}
	LONGITUDINAL_DIFFUSION i, Dhpca*Volume [i] {HPCA}
	LONGITUDINAL_DIFFUSION i, Dcahpca*Volume [i] {CaHPCA}
	LONGITUDINAL_DIFFUSION i, Dca2hpca*Volume [i] {Ca2HPCA}
	LONGITUDINAL_DIFFUSION i, Dhpca_m*Volume [i] {HPCA_m}
	
	LONGITUDINAL_DIFFUSION i, Dhpca_mut*Volume [i] {HPCA_mut}
	LONGITUDINAL_DIFFUSION i, Dcahpca_mut*Volume [i] {CaHPCA_mut}
	LONGITUDINAL_DIFFUSION i, Dca2hpca_mut*Volume [i] {Ca2HPCA_mut}
	LONGITUDINAL_DIFFUSION i, Dhpca_mut_m*Volume [i] {HPCA_mut_m}
	
	LONGITUDINAL_DIFFUSION i, D_Bufer*Volume [i] {Bufer}
	
	~ ca[0] + pump <-> pumpca (k1Pump*parea*(1e10), k2Pump*parea*(1e10))
	~ pumpca <-> pump + cao (k1Pump*parea*(1e10), k2Pump*parea*(1e10))
	CONSERVE pump + pumpca = TotalPump * parea * (1e10)
	ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea
	~ ca[0] << (-(ica - ica_pmp)/(2*FARADAY*delta))
	
	FROM i=0 TO N-2 {
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
		~ HPCA[i] <-> HPCA[i+1] (Dhpca*frat[i+1], Dhpca*frat[i+1])
		~ CaHPCA[i] <-> CaHPCA[i+1] (Dcahpca*frat[i+1], Dcahpca*frat[i+1])
		~ Ca2HPCA[i] <-> Ca2HPCA[i+1] (Dca2hpca*frat[i+1], Dca2hpca*frat[i+1])
		
		~ Bufer[i] <-> Bufer[i+1] (D_Bufer*frat[i+1], D_Bufer*frat[i+1])
		
		~ HPCA_mut[i] <-> HPCA_mut[i+1] (Dhpca_mut*frat[i+1], Dhpca_mut*frat[i+1])
		~ CaHPCA_mut[i] <-> CaHPCA_mut[i+1] (Dcahpca_mut*frat[i+1], Dcahpca_mut*frat[i+1])
		~ Ca2HPCA_mut[i] <-> Ca2HPCA_mut[i+1] (Dca2hpca_mut*frat[i+1], Dca2hpca_mut*frat[i+1])
	}
	
	FROM i=0 TO N-1 {
		~ ca[i] + HPCA[i] <-> CaHPCA[i] (k1HPCA*Volume[i], k2HPCA*Volume[i])
		~ ca[i] + CaHPCA[i] <-> Ca2HPCA[i] (k3HPCA*Volume[i], k4HPCA*Volume[i])
		
		~ ca[i] + Bufer[i] <-> CaBufer[i] (k1bufer*Volume[i], k2bufer*Volume[i])
		
		~ ca[i] + HPCA_mut[i] <-> CaHPCA_mut[i] (k1HPCA_mut*Volume[i], k2HPCA_mut*Volume[i])
		~ ca[i] + CaHPCA_mut[i] <-> Ca2HPCA_mut[i] (k3HPCA_mut*Volume[i], k4HPCA_mut*Volume[i])
	}
	
	~ Ca2HPCA[0] <-> HPCA_m[0] (k7HPCA*Volume[0], k8HPCA*Volume[0])
	~ Ca2HPCA_mut[0] <-> HPCA_mut_m[0] (k7HPCA_mut*Volume[0], k8HPCA_mut*Volume[0])
	
	cai = ca[0]+ca[1]+ca[2]+ca[3]
	
	HPCA_z = HPCA[0]*Volume[0]+HPCA[1]*Volume[1]+HPCA[2]*Volume[2]+HPCA[3]*Volume[3]
	CaHPCA_z = CaHPCA[0]*Volume[0]+CaHPCA[1]*Volume[1]+CaHPCA[2]*Volume[2]+CaHPCA[3]*Volume[3]
	Ca2HPCA_z = Ca2HPCA[0]*Volume[0]+Ca2HPCA[1]*Volume[1]+Ca2HPCA[2]*Volume[2]+Ca2HPCA[3]*Volume[3]
	HPCA_m_z = HPCA_m[0]*Volume[0]
	HPCA_tot_z = HPCA_z + CaHPCA_z + Ca2HPCA_z + HPCA_m_z
	
	HPCA_mut_z = HPCA_mut[0]*Volume[0]+HPCA_mut[1]*Volume[1]+HPCA_mut[2]*Volume[2]+HPCA_mut[3]*Volume[3]
	CaHPCA_mut_z = CaHPCA_mut[0]*Volume[0]+CaHPCA_mut[1]*Volume[1]+CaHPCA_mut[2]*Volume[2]+CaHPCA_mut[3]*Volume[3]
	Ca2HPCA_mut_z = Ca2HPCA_mut[0]*Volume[0]+Ca2HPCA_mut[1]*Volume[1]+Ca2HPCA_mut[2]*Volume[2]+Ca2HPCA_mut[3]*Volume[3]
	HPCA_mut_m_z = HPCA_mut_m[0]*Volume[0]
	HPCA_mut_tot_z = HPCA_mut_z + CaHPCA_mut_z + Ca2HPCA_mut_z + HPCA_mut_m_z
}
