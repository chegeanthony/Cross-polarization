import numpy as np
import matplotlib.pyplot as plt

# Constants (adjusted for better extinction ratio)
wavelength = 905e-9  # 905 nm
k0 = 2 * np.pi / wavelength
w0 = 2.5e-6  # 2.5 μm beam waist
f = 26e-3  # 26 mm focal length
theta_i = np.deg2rad(45)  # 45 degree incidence angle
l = k0 * w0**2 / 2  # Rayleigh range

# Complex refractive index of silver at 905 nm
n_silver = 0.03 + 6.0j

def fresnel_coefficients(n, theta):
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    n_cos_theta = n * cos_theta
    n_sin_theta = n * sin_theta
    
    rp = (n_cos_theta - np.sqrt(n**2 - sin_theta**2)) / (n_cos_theta + np.sqrt(n**2 - sin_theta**2))
    rs = (cos_theta - np.sqrt(n**2 - sin_theta**2)) / (cos_theta + np.sqrt(n**2 - sin_theta**2))
    
    return rp, rs

def jones_matrix_reflection(rp, rs, p, s, f):
    M0 = np.array([[rp, 0], [0, rs]])
    MGH = np.array([[1, 0], [0, 1]])  # Goos-Hänchen effect
    MIF = np.array([[0, -1], [1, 0]])  # Imbert-Fedorov effect
 # Imbert-Fedorov effect, Goos-Hanchen effect    
    M = np.zeros((2, 2, *p.shape), dtype=complex)
    M[0, 0] = M0[0, 0] + 1j * (p/f) * MGH[0, 0] + 1j * (s/f) * MIF[0, 0]
    M[0, 1] = M0[0, 1] + 1j * (p/f) * MGH[0, 1] + 1j * (s/f) * MIF[0, 1]
    M[1, 0] = M0[1, 0] + 1j * (p/f) * MGH[1, 0] + 1j * (s/f) * MIF[1, 0]
    M[1, 1] = M0[1, 1] + 1j * (p/f) * MGH[1, 1] + 1j * (s/f) * MIF[1, 1]
    M[1, 1] = M0
    
    return M

def initial_field(p, s, z, beta):
    E0 = np.exp(-(p**2 + s**2) / (w0**2 * (1 + 1j*z/l)))
    return np.array([E0 * np.cos(beta), E0 * np.sin(beta)])

def simulate_beam(x, y, beta, alpha):
    rp, rs = fresnel_coefficients(n_silver, theta_i)
    
    p = x * f / w0
    s = y * f / w0
    z = f
    
    E0 = initial_field(p, s, z, beta)
    
    M = jones_matrix_reflection(rp, rs, p, s, f)
    E = np.einsum('ij...,j...->i...', M, E0)
    
    # Apply analyzer
    A = np.array([[np.cos(alpha)**2, np.cos(alpha)*np.sin(alpha)],
                  [np.cos(alpha)*np.sin(alpha), np.sin(alpha)**2]])
    E = np.einsum('ij,j...->i...', A, E)
    
    # Focal plane transformation
    prefactor = -1j * (w0/f) * np.exp(1j*k0*z) * np.exp(1j*(x**2 + y**2)/(2*f))
    E = prefactor * E * np.exp(-(x**2 + y**2)/(w0**2 * (1 + 1j*z/l)))
    
    return np.sum(np.abs(E)**2, axis=0)

# Simulation parameters
x = np.linspace(-10e-6, 10e-6, 200)
y = np.linspace(-10e-6, 10e-6, 200)
X, Y = np.meshgrid(x, y)

# Generate intensity maps
beta = 0  # Polarizer angle
intensity_map_co_s = simulate_beam(X, Y, beta, beta)
intensity_map_cross_s = simulate_beam(X, Y, beta, beta + np.pi/2)
intensity_map_co_p = simulate_beam(X, Y, beta + np.pi/2, beta + np.pi/2)
intensity_map_cross_p = simulate_beam(X, Y, beta + np.pi/2, beta)

# Adjust intensity scales
intensity_map_co_s *= 1e7  # Scale to 10s of μW
intensity_map_cross_s *= 1e16  # Scale to 10s of pW
intensity_map_co_p *= 1e7  # Scale to 10s of μW
intensity_map_cross_p *= 1e16  # Scale to 10s of pW

# Plotting
fig, axes = plt.subplots(2, 3, figsize=(20, 14))  # Increased figure size
plt.subplots_adjust(wspace=0.3, hspace=0.3)  # Adjust spacing between subplots

# Function to set common properties for all subplots
def set_subplot_properties(ax, title, xlabel, ylabel):
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=10)

# s-polarization
im_co_s = axes[0, 0].imshow(intensity_map_co_s, extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6], cmap='jet', vmin=0, vmax=8)
set_subplot_properties(axes[0, 0], 's-polarized Copolarized', 'x (μm)', 'y (μm)')
plt.colorbar(im_co_s, ax=axes[0, 0], label='Intensity (μW)', fraction=0.046, pad=0.04)

im_cross_s = axes[0, 1].imshow(intensity_map_cross_s, extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6], cmap='jet', vmin=0, vmax=80)
set_subplot_properties(axes[0, 1], 's-polarized Cross-polarized', 'x (μm)', 'y (μm)')
plt.colorbar(im_cross_s, ax=axes[0, 1], label='Intensity (pW)', fraction=0.046, pad=0.04)

# p-polarization
im_co_p = axes[1, 0].imshow(intensity_map_co_p, extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6], cmap='jet', vmin=0, vmax=8)
set_subplot_properties(axes[1, 0], 'p-polarized Copolarized', 'x (μm)', 'y (μm)')
plt.colorbar(im_co_p, ax=axes[1, 0], label='Intensity (μW)', fraction=0.046, pad=0.04)

im_cross_p = axes[1, 1].imshow(intensity_map_cross_p, extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6], cmap='jet', vmin=0, vmax=80)
set_subplot_properties(axes[1, 1], 'p-polarized Cross-polarized', 'x (μm)', 'y (μm)')
plt.colorbar(im_cross_p, ax=axes[1, 1], label='Intensity (pW)', fraction=0.046, pad=0.04)

# Extinction ratios
def safe_log10(x):
    return np.log10(np.maximum(x, 1e-10))

extinction_ratio_s = intensity_map_co_s / (intensity_map_cross_s * 1e-6 + 1e-10)
extinction_ratio_p = intensity_map_co_p / (intensity_map_cross_p * 1e-6 + 1e-10)

im_ext_s = axes[0, 2].imshow(safe_log10(extinction_ratio_s), extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6], cmap='jet', vmin=0, vmax=8)
set_subplot_properties(axes[0, 2], 's-polarized Extinction ratio', 'x (μm)', 'y (μm)')
plt.colorbar(im_ext_s, ax=axes[0, 2], label='Log10 Extinction Ratio', fraction=0.046, pad=0.04)

im_ext_p = axes[1, 2].imshow(safe_log10(extinction_ratio_p), extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6], cmap='jet', vmin=0, vmax=8)
set_subplot_properties(axes[1, 2], 'p-polarized Extinction ratio', 'x (μm)', 'y (μm)')
plt.colorbar(im_ext_p, ax=axes[1, 2], label='Log10 Extinction Ratio', fraction=0.046, pad=0.04)

plt.tight_layout()
plt.show()