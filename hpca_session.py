from neuron import h, gui
from neuron.units import ms, mV
import matplotlib.pyplot as plt
from matplotlib import style
from matplotlib import style, rc
import numpy as np
import pickle

style.use("seaborn")

font = {'family' : 'monospace',
                'weight' : 'normal',
                        'size'   : 12}

rc('font', **font);
plt.rc('axes', titlesize=15);     # fontsize of the axes title
plt.rc('axes', labelsize=14);

h.DCa_hpca = .22
h.Dhpca_hpca = .05
h.Dcahpca_hpca = .05
h.Dca2hpca_hpca = .05
h.Dhpca_m_hpca = .001
h.k1HPCA_hpca = 40
h.k2HPCA_hpca = .01
h.k3HPCA_hpca = 40
h.k4HPCA_hpca = .01
h.k7HPCA_hpca = .01
h.k8HPCA_hpca = .002
h.TotalHPCA_hpca = .03821

for i in range(10, 21):
    sec = h.dend11[i]
    sec.insert("hpca")
    sec.insert("Ca_HVA")
    sec.insert("cal")
    sec.uninsert("cad")
    print('%s: %s' % (sec, ', '.join(sec.psection()['density_mechs'].keys())))

h.celsius = 23
