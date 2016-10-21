# -*- coding: utf-8 -*-

# (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
# INRIA Rhone-Alpes
# STDP model : parameters' values

# ####### Units #######
# time = sec
# frequency = Hz
# consentraion = microM (if not specified differently)
# conductance = nS
# current = pA
# voltage = mV
# capacitance = nF
# temperature = K
# permeability = microL/s
# FARADAY constant in kCol/mol
# R in J K−1 mol−1

# Spine/soma assumed to be spherical
# Spine = SOMA
# spine_radius=8.0 # micrometers
# spine_volume = 4.0*pi/3.0*spine_radius**3*1e-9  # (in microL) ,
#   1 microm^3 = (1e-6)**3/1e-3 L = 1e-15 L = 1e-9 microL
# spine_surface_dm2 = 4*pi*pow(3*spine_volume*1e-6/(4*pi), 2.0/3.0)
#   (in decim^2) = 4*np.pi*pow(3*1e-15/(4*np.pi), 2.0/3.0) decim^2
# spine_surface_cm2 = spine_surface_dm2*100

from collections import OrderedDict

paramets = OrderedDict([
    ("caL13", OrderedDict([
        ("on", 1),
        ("pcaLbar", 1.02e-06),
        ("mslope", -6.7),
        ("hshift", 0.0),
        ("vm", -8.124),
        ("mshift", 0.0),
        ("hslope", 11.9),
        ("kpr", 31.4),
        ("c", 39.8),
        ("k", 9.005),
        ("hvhalf", -13.4),
        ("mvhalf", -33),
        ("cpr", 990.0),
        ("htau", 0.0443),
        ("hqfact", 3.0),
        ("qfact", 3.0)
    ])),
    ("TRPV1", OrderedDict([
        ("on", 1),
        ("gTRPV1", 0.0003),
        ("p_ca", 2.23204558561e-49),
        ("C", 23367.0),
        ("D", 1100.0),
        ("DH", 205000.0),
        ("KD", 0.5),
        ("J0", 0.0169),
        ("L", 0.00042),
        ("P", 750.0),
        ("z", 0.6),
        ("DS", 615.0),
        ("K", 0.00182634305618)
    ])),
    ("I_to_Ca_flux", OrderedDict([
        ("NMDA", 70.0),
        ("VDCC", 84.0),
        ("TRPV1", 310)
    ])),
    ("CaBuff", OrderedDict([
        ("BT", 4.5),
        ("Cab", 0.1),
        ("tauCab", 0.007),
        ("KdB", 0.5)
    ])),
    ("CICR", OrderedDict([
        ("a2", 0.5),
        ("ver", 8.0),
        ("d3", 0.9434),
        ("rc", 4.0),
        ("rl", 0.1),
        ("rhoER", 0.3),
        ("d5", 0.12),
        ("d2", 3.049),
        ("ker", 0.05),
        ("d1", 0.13)
    ])),
    ("IP3", OrderedDict([
        ("kappad", 1.5),
        ("kdelta", 0.1),
        ("r5p", 0.2),
        ("v3k", 0.001),
        ("kd", 1.5),
        ("vdelta", 0.02),
        ("k3", 1.0),
        ("kr", 1.3),
        ("kp", 10.0),
        ("vbeta", 0.8),
        ("kpi", 0.6),
        ("n3", 1.0),
    ])),
    ("DGLandDAG", OrderedDict([
        ("KDGL", 30.0),
        ("kMAGL", 0.5),
        ("rDGL", 20000.0),
        ("kDAGK", 2)
    ])),
    ("KandP_on_DAGLP", OrderedDict([
        ("nK", 6),
        ("rP", 380),
        ("rK", 50)
    ])),
    ("ECb", OrderedDict([
        ("on", 1),
        ("CB1R_on", 1),
        ("kCB1R", 3000.0),
        ("alphaAEACB1", 0.1),
        ("P1", 1e-09),
        ("P2", 1e-05),
        ("P3", 7),
        ("P4", 2.0),
        ("vATAEA", 0.2),
        ("LTDstart", 0.027),
        ("LTDstop", 0.047),
        ("LTDMax", 0.65),
        ("LTPstart", 0.086),
        ("LTPMax", 13.5425),
        ("KFAAH", 1.0),
        ("vFAAH", 4.0),
    ])),
    ("ECb_smooth", OrderedDict([
        ("on", 0),
        ("K", 0.0007),
        ("n", 2),
        ("kw", 10),
        ("tau", 0.0001),
        ("kadd", 0),
        ("kn", 1)
    ])),
    ("post_CaMKII_plast", OrderedDict([
        ("kpkaI1", 4.67),
        ("PP10", 0.2),
        ("KM", 0.4),
        ("K5", 0.1),
        ("k12", 6000.0),
        ("k11", 500.0),
        ("kcan0I1", 0.05),
        ("km11", 0.1),
        ("CaMKT", 16.6),
        ("k7", 6.0),
        ("k6", 6.0),
        ("ncanI1", 3.0),
        ("I10", 1.0),
        ("CaMT", 0.07052),
        ("KdpkaI1", 0.159),
        ("kcanI1", 20.5),
        ("KdcanI1", 0.053),
        ("kpka0I1", 0.0025),
        ("Ka3", 0.32),
        ("Ka2", 0.025),
        ("Ka1", 0.1),
        ("Ka4", 0.4),
        ("npkaI1", 3.0)
    ])),
    ("DA", OrderedDict([
        ("gamma1DA", 0.7),
        ("gamma2DA", 0.07),
        ("DA", 0.01)
    ])),
    ("AMPA", OrderedDict([
        ("on", 1),
        ("gAMPA", 5.1),
        ("Epsilon", 0.0),
        ("Beta", 190.0),
        ("Alpha", 1.02),
        ("Gamma", 0.0)
    ])),
    ("NMDA", OrderedDict([
        ("on", 1),
        ("gNMDA", 1.53),
        ("p_ca", 2.08324254657e-46),
        ("Mg", 1.0),
        ("Alpha", 0.072),
        ("Beta", 100.0)
    ])),
    ("mem", OrderedDict([
        ("EL", -70.0),
        ("gL", 10.0),
        ("Cm", 0.1)
    ])),
    ("CB1R", OrderedDict([
        ("on", 1.0),
        ("Alpha", 0.240194904182),
        ("Beta", 11.0718971839),
        ("Gamma", 416.378884767),
        ("Epsilon", 0.0477956844649)
    ])),
    ("common", OrderedDict([
        ("R", 8.3144621),
        ("RT", 2553.78703401),
        ("T", 307.15),
        ("F", 96.5),
        ("zS", 2.0),
        ("Ca_out", 5000.0)
    ])),
    ("integration", OrderedDict([
        ("t_step", 0.001),
        ("t_start", 0.0),
        ("t_end", 250.0),
        ("ATOL", 1e-07),
        ("RTOL", 1e-07),
        ("MXSTEP", 1000),
        ("HMAX", 50)
    ])),
    ("action", OrderedDict([
        ("APdur", 0.03),
        ("action_as_VDCC", False),
        ("DPmax", 495.0),
        ("APmax", 7020.0),
        ("tausbAP", 0.001)
    ])),
    ("Glu_release", OrderedDict([
        ("steadyrise_on", 0),
        ("Glumax", 2000.0),
        ("BaseLevel", 0.0),
        ("tauGlu", 0.005)
    ])),
    ("stimulation", OrderedDict([
        ("regular_on", 1),
        ("tsdt", 0.015),
        ("tpost", 0.5),
        ("post_on", 1),
        ("dt_stim", 0.02),
        ("num_stim", 100),
        ("pre_on", 1),
        ("Freq", 1.0),
        ("tables_step", 5e-05)
    ]))
])
