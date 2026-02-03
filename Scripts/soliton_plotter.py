import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("soliton.dat", comments="#", delimiter=",")

r = data[:, 0]      
f = data[:, 1]      
phi = data[:, 2]   

plt.figure(figsize=(8, 5))
plt.plot(r, f, label="f(r)")
plt.plot(r, phi, label="phi(r)")
plt.xlabel("r")
plt.ylabel("Value")
plt.title("Soliton Profile")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()