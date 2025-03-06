import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint

#############################
# DETERMINISTIC DYNAMICS CODE
#############################

#Parameter values
fm=1.0
fp=1.0
nequ=2
alpham=100.0*fm
deltam=1.0*fm
alphap=10.0*fp
deltap=0.1*fp

#initial conditions
y0 = [] 
y0.append(0.0)
y0.append(0.0)

#Derivatives
def dy_dt(y,t):    
    dy = []
    dy.append(alpham-deltam*y[0])
    dy.append(alphap*y[0]-deltap*y[1])
   
    return dy   

#Integration of ODEs
t = np.linspace(start=0, stop=50,num=1000)
prot = np.linspace(start=0, stop=15000,num=1000)
sol = odeint(func=dy_dt, y0=y0, t=t)
sol = pd.DataFrame(sol, columns=["mRNA","protein"])

sol.index = t
#sol.index = sol.index/60

sol["mnull"] = alpham/deltam
sol["pnull"] = deltap*prot/alphap

#Plotting
fig, ax = plt.subplots(3,figsize=(12,8))
#sol.plot(subplots=True,layout=(3,1),figsize=(10,6),ax=ax, xlabel="Time (min)", ylabel = "Concentration (nM)")
ax[0].plot(sol.index,sol.mRNA, label="mRNA (nM)")
ax[0].set_xlabel("time (min)")
ax[0].set_ylabel("concentration (nM)")
ax[0].legend()

ax[1].plot(sol.index,sol.protein, label="protein (nM)")
ax[1].set_xlabel("time (min)")
ax[1].set_ylabel("concentration (nM)")
ax[1].legend()


box = ax[2].get_position()
ax[2].set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])

ax[2].plot(sol.protein,sol.mRNA, c="violet", label="protein(nM)")
ax[2].set_xlabel("protein (nM)")
ax[2].set_ylabel("mRNA (nM)")

ax[2].plot(prot,sol.mnull, c="green", label="mnull(x)")
ax[2].plot(prot,sol.pnull, c="blue", label="pnull(x)")

x,y = np.meshgrid(np.linspace(start=0, stop=15000,num=10),np.linspace(start=0, stop=200,num=10))

u = (alphap*y-deltap*x)*fp/fm
v = (alpham-deltam*y)*fp/fm

# vector: (alphap*y-deltap*x):(alpham-deltam*y)

ax[2].quiver(x,y,u,v,color='grey',  angles='xy', scale_units='xy', scale=2, width=0.002, label="Vector field")


ax[2].legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), fancybox=True, shadow=True, ncol=5)

fig.subplots_adjust(top=0.8)
fig.align_ylabels()
plt.tight_layout()
plt.show()



##########################
# STOCHASTIC DYNAMICS CODE
##########################

# Set the parameters for the simulation
alpha_m = 1000  # nM/min
delta_m = 1    # min^-1
alpha_p = 1   # min^-1
delta_p = 0.1  # min^-1
volume = 1     # cubic micrometer
tmax = 100     # min
dt = 0.01      # min
#k = 1000 #negative feedback

# Initialize the simulation
t = 0
m = 0
p = 0
m_list = []
t_list = []
p_list = []

# Run the simulation
while t < tmax:
    # Calculate the feedback function
    #f = k**2 / (k**2 + p**2)
    #rm = 10100 * f
    
    # Calculate the reaction rates
    rm = alpha_m * volume
    rp = alpha_p * m
    dm = delta_m * m
    dp = delta_p * p
    
    # Calculate the total rate of all reactions
    rtotal = rm + rp + dm + dp
    
    # Choose the time until the next reaction event
    tau = np.random.exponential(1/rtotal)
    
    # Choose which reaction event occurs
    r = np.random.uniform(0, rtotal)
    if r < rm:
        m += 1
    elif r < rm + rp:
        p += 1
    elif r < rm + rp + dm and m > 0:
        m -= 1
    elif r < rm + rp + dm + dp and p > 0:
        p -= 1
    
    # Update the time and store the current concentrations
    t += tau
    m_list.append(m / volume)
    p_list.append(p / volume)
    t_list.append(t)
    
# Plot the time evolution of the mRNA concentration
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(t_list, m_list,"C0",label="Stochastic dynamics")
ax.plot(sol.index,sol.mRNA,"k",label="Deterministic dynamics",linewidth=2)
plt.xlim([-1, 50])
ax.set_xlabel('Time (min)')
ax.set_ylabel('Concentration (nM)')
ax.set_title('Time evolution of the mRNA concentration (2)')
ax.legend(loc='lower right')
plt.show()

# Plot the time evolution of the protein concentration
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(t_list, p_list,"C1",label="Stochastic dynamics")
ax.plot(sol.index,sol.protein,"k",label="Deterministic dynamics",linewidth=2)
plt.xlim([-1, 50])
ax.set_xlabel('Time (min)')
ax.set_ylabel('Concentration (nM)')
ax.set_title('Time evolution of the protein concentration (2)')
ax.legend(loc='lower right')
plt.show()

# HISTOGRAM PLOTS

# Calculate the mean and standard deviation of the mRNA concentration data
mean_m = np.mean(m_list)
mean_p = np.mean(p_list)
lan_m = np.mean(m)
lan_p = np.mean(p)
std_m = np.sqrt(lan_m)/volume
std_p = np.sqrt(lan_p)/volume
print("Mean concentration of mRNA: ",mean_m,"nM, with a standard deviation of: ",std_m)
print("Mean concentration of protein: ",mean_p,"nM, with a standard deviation of: ",std_p)

# Define the Gaussian distribution
x_m = np.linspace(0, 2000, 1000)
gaussian_m = 1 / (std_m * np.sqrt(2 * np.pi)) * np.exp(-((x_m - mean_m)**2) / (2 * std_m**2))

x_p = np.linspace(0, 20000, 1000)
gaussian_p = 1 / (std_p * np.sqrt(2 * np.pi)) * np.exp(-((x_p - mean_p)**2) / (2 * std_p**2))

# Plot the histogram of the mRNA concentration data
fig, ax = plt.subplots(figsize=(9,5))
ax.hist(m_list, density=True, bins=50, rwidth=0.9, alpha=0.5, color="C0", label='Simulation')
fig.tight_layout()
# Plot the Gaussian distribution
ax.plot(x_m, gaussian_m, 'k--', label='Gaussian limit of Poisson')

plt.xlim([600, 1400])
ax.set_xlabel('mRNA Concentration (nM)')
ax.set_ylabel('Probability Density')
ax.set_title('Probability density and gaussian distribution for mRNA (Constitutive Expression)')
ax.legend()
plt.show()

# Plot the histogram of the mRNA concentration data
fig, ax = plt.subplots(figsize=(9,5))
ax.hist(p_list, density=True, bins=30, rwidth=0.9, alpha=0.5, color="C1", label='Simulation')
fig.tight_layout()
# Plot the Gaussian distribution
ax.plot(x_p, gaussian_p, 'k--', label='Gaussian limit of Poisson')

plt.xlim([9000, 11000])
ax.set_xlabel('Protein Concentration (nM)')
ax.set_ylabel('Probability Density')
ax.set_title('Probability density and gaussian distribution for protein (Constitutive Expression)')
ax.legend()
plt.show()