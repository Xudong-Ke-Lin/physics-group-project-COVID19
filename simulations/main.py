from walkers_positions import walkers_positions

positions = walkers_positions


"""
Created on Sunday 28th February 2021

@author: joel pendleton & xudongke
"""



# imports
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.ndimage import gaussian_filter

# start counter
counter = 0

# Number of persons e.g. sweeps over 6400 10000
Nps = [6400, 10000]
# fraction of infected
sicks = [0.005, 0.01, 0.02]
# coughing or not e.g. none once/hour twice/hour
Pcoughs = [0, 1 / 3600, 2 / 3600]
# possibility to investigate superemitters
superswitch = [0]
# walking speed
Upaves = [0.01, 0.1, 0.5, 1]
# ventilation timescale in seconds
taus = [100]
# diffusivity including aerosol mixing due to turbulence
Ds = [0.05]
# mask effectivity
mask_effs_out = [0, 0.6, 0.9]
mask_effs_in = [0, 0.3, 0.5]
# far-UVC susceptibility constant [cm2 mJ−1]
Zs = [0, 4.1]
# fraction of asymptomatic
asymptomatics = [0.2, 0.4]
# fever fraction in symptomatic
fevers = [0.34, 0.45, 0.88]
# use of Non-contact Infrared thermometers (NCIT)
NCITs = [0, 1]
# effectivity of Non-contact Infrared thermometers (NCIT)
NCIT_effs = [0.82]

###################################################
# System dimensions
###################################################
Nx = 20
Ny = 22
#  meters 20m*10m
Lx = 20
Ly = 22
# spacing
dx = Lx / Nx
dy = Ly / Ny


# reference for indices
radius = Nx/(2*Lx)
x = np.arange(0,Nx+1)
y = np.arange(0,Ny+1)
y2d, x2d=np.meshgrid(y,x)
# each element contains list of neighbour indices
indices_array = np.zeros((Ny+1,Nx+1),dtype=object)
# area of neighbours
area_array = np.zeros((Ny+1,Nx+1))
# create neighbour indices for each position


for x0 in x:
    for y0 in y:
        array = np.zeros_like(x2d)
        array[np.where(np.sqrt((x2d-x0)**2+(y2d-y0)**2)<radius)]=1
        # boundary conditions
        if x0<radius:
            array[:,0]=0
        if x0>Lx-radius-1:
            array[:,-1]=0
        if y0<radius:
            array[0,:]=0
        if y0>Ly-radius-1:
            array[-1,:]=0
        indices_array[y0,x0]=np.where(array==1)
        area_array[y0,x0]=len(indices_array[y0,x0])

###################################################
# UV field intensity
###################################################
# grid for intensity field
x_field = np.arange(0,Nx+1,dx)
y_field = np.arange(0,Ny+1,dy)
x2d_field,y2d_field=np.meshgrid(x_field,y_field)
# position of UV lamp
x0_field = x_field[int(len(x_field)/2)]
y0_field = y_field[int(len(y_field)/2)]

# far-UVC maximum intensity [mJ cm-2 s-1]
Ep_0 = 0.0014


denominator = ((x2d_field-x0_field)**2+(y2d_field-y0_field)**2)
denominator = np.where(denominator <= 0, 1, denominator)

# far-UVC intensity field
Ep = Ep_0/denominator
# apply Gaussian filter with sigma=4
Ep = gaussian_filter(Ep, 4)



###################################################
# Parameters selection
###################################################
aaa = 0
# 2% of infected people
bbb = 2
# 2 coughs per hour
ccc = 2
ddd = 0
eee = 2
fff = 0
ggg = 0
# mask eff out 60%
hhh_out = 1
# mask eff in 30%
hhh_in = 1
# with UV light
iii = 1
# 40% asymptomatic
jjj = 1
# 45% fever in symptomatic
kkk = 1
# do not use NCIT
lll = 0
mmm = 0

# total number of people
Np = Nps[aaa]
# fraction of sick persons
sick_fraction = sicks[bbb]
# symptomatic fraction
symptomatic = 1 - asymptomatics[jjj]
# do not use NCIT
if NCITs[lll] == 0:
    # total number of infected people
    total_sick = round(Np * sick_fraction)
# use NCIT
else:
    total_sick = round(Np * sick_fraction * (1 - symptomatic * fevers[kkk] * NCIT_effs[mmm]))

# number of people in the room, assume 100 rooms
N = int(Np / 100)
# Hypergeometric distribution
rv = stats.hypergeom(Np, total_sick, N)
# sick of people in room
N_sicks = np.arange(0, N)
pmf_sick = rv.pmf(N_sicks)
# sick of people in room with nonzero probability
N_sicks = N_sicks[np.round(pmf_sick * 100) != 0]
print(f'possible number of infected people in room: {N_sicks}')
# choose the highest number of sick people in room
N_sick = N_sicks[-1]
# you can pick the fraction how many emit more than others
superemitter = 0.1
super_sick = superswitch[ddd] * max(1, round(superemitter * N_sick));
# sp = 0 -> healthy sp=1 -> sick
sp = np.zeros(N)
sp[0:N_sick] = 1
# super emitter
se = np.zeros(N)
se[np.where(sp == 1)][0:super_sick] = 1
# the healthy persons gain a dose
# i.e. number of aerosols
dose = np.zeros(N)

# velocity of each individual
Upave = Upaves[eee]
# timescale [s] of removing aerosols upwards due to ventilation
tau = taus[fff]
# m^2/s  diffusivity in x,y plane direction due to turbulence
D = Ds[ggg]
# mask effectivity
mask_eff_out = mask_effs_out[hhh_out]
mask_eff_in = mask_effs_in[hhh_in]
# UV light: 0 if not, 1 if yes
Z = Zs[iii]

# concentration, assume zerogradient bc's
C = np.zeros((Ny+1, Nx+1))

# source term in units [# of particles / 0.01 m^3 ]
S = np.zeros((Ny+1, Nx+1))

simuTime = 2000  # seconds
dt = 1  # seconds
simuSteps = round(simuTime / dt)

###################################################
# parameters of sick individuals
###################################################
# lambda aerosols [1/s] e.g. 5 means 5 per second into control volume
lambdanormal = 5 * (1 - mask_eff_out)

# how many particles come from a single cough, assume this is
# immediately spread to the control volume during dt
lambdacough = 40000 * (1 - mask_eff_out)

# how many m^3 you inhale in a second, here 20 l per 60 sec
Vinh = 0.001 * 20 / 60
# volume in which the droplets locate
Vvol = area_array * dx * dy * 1
# Here. e.g. 1 cough per 3600 sec; then during dt probability dt*Pcough
Pcough = Pcoughs[ccc]

# PDF that accumulates the dose
# PDF = np.zeros((1,1+lambdacough*100))
aveinf = []
maxinf = []
aveC = []
Np = 20

# run single simulation with fixed parameters
for k in range(0, simuSteps):

    # zero source term
    S = 0 * S


    for s in range(0, Np):

        xp = positions[s, k, 0]
        yp = positions[s, k, 1]

        # make person coordinate an integer
        xpint = int(xp)
        ypint = int(yp)
        # sick persons generate lambda*dt particles per volume during dt
        if sp[s] == 1:

            S[ypint, xpint] += (
                        lambdanormal * dt + se[s] * lambdanormal * 9 * dt) / area_array[ypint, xpint]
            # coughing occurs at probablity Pcough*dt
            if np.random.random() < Pcough * dt:
                S[xpint, ypint] += lambdacough / area_array[ypint, xpint]
        # healthy individuals accumulate a dose
        else:

            dose[s] += np.sum(C[ypint, xpint]) *(Vinh/Vvol[ypint,xpint])*dt*(1-mask_eff_in)


    # diffusion equation with source and sink terms for aerosol concentration
    # following the paper using expl.Euler and CD2 method (finite difference)
    C = C + dt * D * (np.roll(C, 1, axis=0) + np.roll(C, -1, axis=0) + np.roll(C, 1, axis=1) + np.roll(C, -1,
                                            axis=1) - 4 * C) / dx ** 2
    # hard walls boundary conditions
    C[0, :] = 0
    C[-1, :] = 0
    C[:, 0] = 0
    C[:, -1] = 0
    # [C] number of aerosols in the control volume
    C = C + S - dt * C / tau - Z * Ep * C


    if k % 100 == 0:
        plt.clf()
        plt.plot(xp, yp, '.')
        plt.title('Walkers in a public place')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.xlim(0, Lx)
        plt.ylim(0, Ly)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.pause(1)

        plt.clf()
        plt.imshow(C / 1000, extent=[0, Lx, 0, Ly])
        plt.title('Number of aerosols per liter')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.colorbar()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.pause(1)

        # append dose and concentration
        aveinf.append(np.mean(dose[sp == 0]))
        maxinf.append(np.max(dose))
        aveC.append(np.mean(C))

np.savetxt(f'data_{aaa}_{bbb}_{ccc}_{ddd}_{eee}_{fff}_{ggg}_{hhh_out}_{hhh_in}_{iii}_{jjj}_{kkk}_{lll}_{mmm}.txt'
           , np.column_stack((aveinf, maxinf, aveC)))
