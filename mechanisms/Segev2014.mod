:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA
	USEION ca READ eca WRITE ica
	RANGE gCa_HVAbar, ica, shift
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gCa_HVAbar = 0.000002 (S/cm2)
	a01=0.055
	a02=0.94
	a03=0.000457
	a04=0.0065
	a05=2.02	:15
	a06=2.26	:5
	shift=5	(mV)
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa_HVA	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
	hInf
	hTau
	hAlpha
	hBeta
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gCa_HVA = gCa_HVAbar*m*m*h
	ica = gCa_HVA*(v-eca)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
	LOCAL q10
	UNITSOFF
		q10 = 3^((celsius - 23)/10)
        if((v == -27) ){
            v = v+0.0001
        }
		mAlpha =  (a01*(-27-(v+shift)))/(exp((-27-(v+shift))/3.8) - 1)					:0.055
		mBeta  =  (a02*exp((-75-(v+shift))/17))									:0.94
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = a05/(q10*(mAlpha + mBeta))		:ORIGINAL 15/(mAlpha + mBeta)
		hAlpha =  (a03*exp((-13-(v+shift))/50))									:0.000457
		hBeta  =  (a04/(exp((-(v+shift)-15)/28)+1))								:0.0065
		hInf = hAlpha/(hAlpha + hBeta)
		hTau = a06/(q10*(hAlpha + hBeta))		:ORIGINAL 5/(hAlpha + hBeta)
	UNITSON
}
