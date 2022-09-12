#alternate code which can be used to simulate individual modes
import numpy as np
import matplotlib.pyplot as plt

n = 0;                  # Beam order
m = 0;                  # Beam order
w0 = 2.0;               # Beam waist
k = 2*np.pi/532.0e-9;   # Wavenumber of light

zR = k*w0**2.0/2;       # Calculate the Rayleigh range

# Setup the cartesian grid for the plot at plane z
z = 0.0;
[xx, yy] = np.meshgrid(np.linspace(-5, 5), np.linspace(-5, 5));

def hermiteH(n, x):
    if n == 0:
        return 1;
    elif n == 1:
        return 2*x;
    else:
        return 2*x*hermiteH(n-1, x) - 2*(n - 1)*hermiteH(n - 2, x);

U00 = 1.0/(1 + 1j*z/zR)*np.exp(-(xx**2 + yy**2)/w0**2/(1 + 1j*z/zR));
Hn = hermiteH(n, xx);
Hm = hermiteH(m, yy);

U = U00*Hn*Hm*np.exp(-1j*(n + m)*np.arctan(z/zR));

plt.figure()
plt.title('Intensity')
plt.pcolor(abs(U)**2);
plt.axis('equal')

plt.figure()
plt.title('Phase')
plt.pcolor(np.angle(U)**2);
plt.axis('equal')

plt.show()
