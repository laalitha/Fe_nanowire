# read data from file
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def murnaghan(V, E0, B0, BP, V0):
    'From PRB 28,5480 (1983'

    E = E0 + B0 * V / BP * (((V0 / V)**BP) / (BP - 1) + 1) - V0 * B0 / (BP - 1)
    return E
V,E = np.loadtxt('EV_data.txt', unpack=True)
E1 = list(E)
E0 = min(E)
V0 = V[list(E).index(min(E))]

print(E1.index(min(E)))

pars, cov = curve_fit(f=murnaghan, xdata=V, ydata=E, p0=[ E0,1 ,4, V0], bounds=(-np.inf, np.inf))

E0 = pars[0]
B0 = pars[1]
V0 = pars[3]

print("Minimum Energy : ", E0)
print("Minimum Volume : ", V0)
print("Bulk Modulus   : ", B0*160.2176621)
#bulk modulus in Gpa
#minimum energy eV/atom
#minimum volume Å^3/atom
plt.xlabel("Atomic Volume(Å^3)")
plt.ylabel("Energy(eV/atom)")
plt.title("Mendelv")
plt.plot(V,E,'bo')
plt.plot(V,E)
plt.show()
#plt.plot(V,vinet(V,*pars))
