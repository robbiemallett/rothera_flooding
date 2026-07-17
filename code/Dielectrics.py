from math import pi

from cmath import exp

 

# %% Constants

PERM_FREE = 8.85419e-12     # Permittivity of free space after Stogryn and Desargant (1985)

 

#freq_dict = {"HydraProbe": 5.0E+07}

#              "C-band RCM": 5.405e+09,

#              "C-band Scatterometer": 5.5e+09,

#              "L-band ALOS-2 PALSAR-2": 1.27e+09,

#              "X-band TerraSAR-X": 9.65e+09,

#              "Ku-band QuickScat Ku": 13.4e+09,

#              "Ku-band Scatterometer": 13.575e+09,

#              "Other": 35.1e+09,

#              "HydraProbe": 5.0E+07}

 

# Begin INPUT >>>>>>>>>>>>>>>

 

# Frequency wanted

freq = 13.575e+09

#freq = 35.1e+09

 

# Frequency at which dielectrics were measured

freq_measured = 5.0E+07

 

# Measured permittivity and loss

Er = 2.3

Eim = 0.2

 

# Calibrations for UofC Steven's meter

# Er = Er - 0.2791

# Eim = Eim - Eim * 0.077937

 

# Assumed or measured snow properties

Ps = 0.320

Ps_dried = 0.320  # initial estimate

Ts = -5

 

# End INPUT <<<<<<<<<<<<<

 

# Salinity of brine

if Ts >= -2.06:

    Sbrine = (1 / (1 - 54.11 / Ts)) * 1000

if Ts < -2.06 and Ts >= -8.2:

    Sbrine = 1.725 - 18.756 * Ts - 0.3964 * Ts**2  # Ulaby pg. 2045

if Ts < -8.2 and Ts > -22.9:

    Sbrine = 57.041 - 9.929 * Ts - 0.16204 * Ts**2 - 0.002396 * Ts**3     # Ulaby pg. 2045

if Ts >= -36.8 and Ts < -22.9:

    Sbrine = 242.94 + 1.5299 * Ts + 0.0429 * Ts**2      # Ulaby page 2045

if Ts >= -43.2 and Ts < -36.8:

    Sbrine = 508.18 + 14.535 * Ts + 0.2018 * Ts**2  #(REF?)

if Ts < -43.2:

    Sbrine = 508.18 + 14.535 * Ts + 0.2018 * Ts**2  # assume same as Ts >= -43.2 and Ts < -36.8

 

# Density of ice, water, brine

Pice = 0.917 - 0.0001403 * Ts   # Pounder 1965 in Backstrom and Eicken 2006

Pwater = 0.9999

Pbrine = 1. + 0.0008 * Sbrine            # Cox and Weeks 1975

   

# Calculate dry snow dielectrics

e_drysnow_p = 1. + 2.55 * Ps_dried

e_drysnow_l = 0.001   # assumed loss for dry snow - need referenced value

e_drysnow = complex(e_drysnow_p, e_drysnow_l)

 

# Calculate brine volume for measured permittivity

Vbsnow = (Er - e_drysnow_p) / 78.65

 

# Recalculate snow properties

Ps_dried = Ps - Vbsnow * Pbrine

porosity = 1. - (Ps - Pbrine * Vbsnow) / Pice            # Denoth (1980)

saturation = Vbsnow / porosity

 

# Using the estimated brine volume from above calculate the permittivity and

# loss at the desired frequency

 

# Depolarization factor (Geldsetzer et al 2009)

sfb = 1.33

 

# Brine dielectrics

e_stat = (939.66 - 19.068 * Ts) / (10.737 - Ts) # Stogryn and Desargant (1985) Eq.10

e_inf = (82.79 + 8.19 * Ts**2) / (15.68 + Ts**2)  # Stogryn and Desargant (1985) Eq.11   ********* compare to water at 4.9 in Tiuri

pit = (0.1099 + 0.13603e-02 * Ts + 0.20894e-03 * Ts**2 +

       0.28167e-05 * Ts**3) * 1e-09   # Stogryn and Desargant (1985) Eq.12

relax_time = pit / (2 * pi)

e_relax_freq = 1 / (2 * pi * relax_time)

 

# Ionic conductivity of brine (Stogryn and Desargant, 1985 Eq.7)

if Ts >= -22.9:

    ion_cond = -Ts * exp(0.5193 + 0.08755 * Ts)

else:

    ion_cond = -Ts * exp(1.0334 + 0.1100 * Ts)

 

# Brine permittivity

e_brine_p = e_inf + ((e_stat - e_inf) / (1. + (pit * freq)**2))   # Ulaby et al (1986)

# Brine loss

e_brine_l = ( ((e_stat - e_inf) * (freq / e_relax_freq))

             / (1+((freq / e_relax_freq)**2)) ) # Tiuri et al 1984, see also Hoekstra and Cappillino 1971 Eq.3

# Brine conductivity term

e_brine_cond = (ion_cond / (2. * pi * freq * PERM_FREE))

 

e_brine = complex(e_brine_p, e_brine_l)

 

# Calculate dielectrics at desired frequency

e_mix = e_drysnow + sfb * Vbsnow * e_brine

 

# Add in conductivity term for loss

cf = 1.89 - 0.14 / (saturation * 100)      # Geldsetzer et al 2009 CRST

if (cf < 1.0):

    cf = 1.0     # outside of data range at very low saturations, so set limit (has very little effect other than avoiding NAN)

 

e_mix_l = (e_mix).imag + (e_brine_cond * ((Vbsnow - 0.0)**cf))   # see Lux 1993 for conductivity

e_mix = complex((e_mix).real, e_mix_l)

 

print(e_mix)