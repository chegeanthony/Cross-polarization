# Beam Simulation Project

## Overview

This project simulates the behavior of a laser beam reflected off a silver surface, taking into account various optical effects such as Fresnel reflection, Goos-Hänchen effect, and Imbert-Fedorov effect. It provides visualizations for both s-polarized and p-polarized light, showing copolarized and cross-polarized intensity distributions as well as extinction ratios.

## Features

- Simulates laser beam reflection with a wavelength of 905 nm
- Considers Fresnel reflection coefficients for s and p polarizations
- Incorporates Goos-Hänchen and Imbert-Fedorov effects
- Generates intensity maps for copolarized and cross-polarized components
- Calculates and visualizes extinction ratios
- Supports both s-polarized and p-polarized incident light

## Requirements

- Python 3.x
- NumPy
- Matplotlib

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/chegeanthony/Cross-polarization.git
   cd Cross-polarization
   ```

2. Install the required packages:
   ```
   pip install numpy matplotlib
   ```

## Usage

Run the main script to generate the simulation and visualizations:

```
python beam_simulation.py
```

This will produce a figure with six subplots:
1. s-polarized Copolarized intensity
2. s-polarized Cross-polarized intensity
3. s-polarized Extinction ratio
4. p-polarized Copolarized intensity
5. p-polarized Cross-polarized intensity
6. p-polarized Extinction ratio

## Customization

You can modify the following parameters in the script to adjust the simulation:

- `wavelength`: Laser wavelength (default: 905 nm)
- `w0`: Beam waist (default: 2.5 μm)
- `f`: Focal length (default: 26 mm)
- `theta_i`: Incidence angle (default: 45 degrees)
- `n_silver`: Complex refractive index of silver at the given wavelength

## Contributing

Contributions to improve the simulation or extend its capabilities are welcome. Please feel free to submit pull requests or open issues for discussion.

## License

[MIT License](LICENSE)

## Contact

For any questions or feedback, please open an issue in the GitHub repository.
