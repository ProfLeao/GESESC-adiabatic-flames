temps = np.linspace(-161.6,537,1000) 

d_mol_real  = PropsSI('Dmolar', 'T', (temps + 273.15), 'P', 101325, "CH4")
d_mol_ideal = dens_molarSI(101325, (temps + 273.15))

plt.plot(temps, d_mol_real, label="101325 Pa - aproximação real")
plt.plot(temps, d_mol_ideal, label="101325 Pa - aproximação ideal")
plt.xlabel("Temperaturas em °C")
plt.ylabel("Densidade molar (mol/m³)")
plt.yscale("log")
plt.title("Relação Densidade Molar do CH4 x Temperatura")
plt.legend()
plt.grid()
plt.show()