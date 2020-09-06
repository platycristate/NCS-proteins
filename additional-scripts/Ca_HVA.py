import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
style.use("ggplot")


# Plotting activation parameter m_inf of Ca_HVA current
a01 = 0.055
shift = 5
a02 = 0.94
q10 = 1 # for t = 23 C
a05 = 2.02
a03 = 0.000457
a04 = 0.0065
a06 = 2.26
def gating_vars(v):
    mAlpha = a01*(-27-(v+shift))/np.exp((-27-(v+shift)/3.8) - 1)
    mBeta = a02*np.exp((-75-(v+shift))/17)
    mInf = mAlpha / (mAlpha + mBeta)
    hAlpha = a03*np.exp((-13-(v+shift))/50)
    hBeta = a04/(((np.exp((-v+shift)-15)/28)+1))
    hInf = hAlpha / (hAlpha + hBeta)
    return [mInf, hInf]

v = np.linspace(-150, 100, 300)
plt.plot(v, gating_vars(v)[0], label="m", lw=2)
plt.plot(v, gating_vars(v)[1], label="h", lw=2)
plt.xlabel("v (mv)")
plt.title("Ca_HVA")
plt.legend()
plt.show()

def taus_vars(v):
    mAlpha = a01*(-27-(v+shift))/np.exp((-27-(v+shift)/3.8) - 1)
    mBeta = a02*np.exp((-75-(v+shift))/17)
    hAlpha = a03*np.exp((-13-(v+shift))/50)
    hBeta = a04/(((np.exp((-v+shift)-15)/28)+1))
    mTau = a05/q10*(mAlpha + mBeta)
    hTau = a06/q10*(hAlpha + hBeta)
    return [mTau, hTau]

v = np.linspace(-150, 50, 300)
plt.plot(v, taus_vars(v)[0], label="m", lw=2)
#plt.plot(v, taus_vars(v)[1], label="h", lw=2)
plt.xlabel("v (mv)")
plt.ylabel("t (ms)")
plt.title("Ca_HVA taus")
plt.legend()
plt.show()

