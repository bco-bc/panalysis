import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
Y = np.arange(0, 8.33332, 1.04166625)
X = np.arange(0, 10.00000, 1.0)

print(f'Shape X: {X.shape}')
print(f'Shape Y: {Y.shape}')

X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

print(f'Shape X: {X.shape}')
print(f'Shape Y: {Y.shape}')
print(f'Shape R: {R.shape}')
print(f'Shape Z: {Z.shape}')

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
