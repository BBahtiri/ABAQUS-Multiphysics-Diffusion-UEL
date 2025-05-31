# ABAQUS Multiphysics Diffusion UEL

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![ABAQUS](https://img.shields.io/badge/ABAQUS-2020%2B-blue.svg)](https://www.3ds.com/products-services/simulia/products/abaqus/)
[![Fortran](https://img.shields.io/badge/Fortran-90%2B-orange.svg)](https://fortran-lang.org/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2018a%2B-red.svg)](https://www.mathworks.com/products/matlab.html)

A coupled multiphysics finite element implementation for moisture diffusion in epoxy materials with stress-assisted transport mechanisms.

## 🔬 Overview

This project implements a **User Element (UEL)** for ABAQUS that simulates the coupled hydro-mechanical behavior of moisture diffusion in polymer materials. The element captures **bidirectional coupling** between mechanical stress fields and moisture transport phenomena, making it a true multiphysics simulation.

### Key Physics

- **Stress-assisted diffusion**: Hydrostatic stress influences moisture transport
- **Coupled field equations**: Simultaneous solution of mechanical and diffusion problems
- **Nonlinear coupling**: Stress gradients create preferential diffusion pathways

## ✨ Key Features

- 🔧 **Multiphysics Coupling**: Stress-assisted diffusion with hydrostatic stress influence
- 🧮 **20-Node Elements**: Quadratic hexahedral elements with 80 DOF (60 mechanical + 20 concentration)
- ⚡ **Monolithic Solution**: Simultaneous solution of mechanical and diffusion fields
- 🔗 **ABAQUS Integration**: Full integration with ABAQUS Standard
- 📊 **Visualization Support**: MATLAB tools for result post-processing
- 🧪 **Material Flexibility**: Configurable material properties for different polymers

## 📐 Governing Equations

### Coupled Diffusion Equation
```
∂c/∂t = ∇·(D∇c) - (D·Vh)/(R·T) ∇·(σh∇c)
         ↑_______↑   ↑_________________↑
      Pure diffusion   Stress-assisted term
```

### Mechanical Equilibrium
```
∇·σ = 0
σ = D:(ε - εswelling)
```

### Coupling Parameter
```
κ = Vh/(R·T) = 8000/(8314.5 × 300) ≈ 3.2 × 10⁻³ MPa⁻¹
```

**Where:**
- `c`: Moisture concentration [mol/mm³]
- `D`: Diffusion coefficient [mm²/s]
- `σh`: Hydrostatic stress [MPa]
- `Vh`: Molar volume of water (8000 mm³/mol)
- `R`: Gas constant (8314.5 J/(mol·K))
- `T`: Temperature (300 K)

## 🚀 Quick Start

### Prerequisites
- **ABAQUS Standard** 2020 or later
- **Intel Fortran Compiler** (part of ABAQUS installation)
- **MATLAB** R2018a or later (for visualization)
- **Windows/Linux** operating system

### Basic Usage

1. **Download the UEL**
   ```bash
   git clone https://github.com/BBahtiri/ABAQUS-Multiphysics-Diffusion-UEL.git
   cd ABAQUS-Multiphysics-Diffusion-UEL
   ```

2. **Run ABAQUS analysis**
   ```bash
   abaqus -standard -job MyAnalysis -user Diffusion_3D.for
   ```

3. **Generate visualization mesh** (MATLAB)
   ```matlab
   % Place your ABAQUS input file as 'CT.inp'
   run('VisualMesh.m')
   ```

4. **Post-process results** in ABAQUS/Viewer

### Element Definition in ABAQUS Input File

```fortran
*USER ELEMENT, NODES=20, TYPE=U1, PROPERTIES=3, COORDINATES=3, VAR=128
1,2,3
1,11

*ELEMENT, TYPE=U1, ELSET=SOLID
1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20

*UEL PROPERTY, ELSET=SOLID
210000., 0.3, 0.0127
! E [MPa], ν [-], D [mm²/s]
```

## 📁 Repository Structure

```
📦 ABAQUS-Multiphysics-Diffusion-UEL/
├── 📄 README.md                    # This file
├── 📄 LICENSE                      # MIT License
├── 📄 .gitignore                   # Git ignore patterns
├── 🔧 Diffusion_3D.for             # Main UEL subroutine (18 KB)
├── 📊 VisualMesh.m                 # MATLAB visualization script (2 KB)
├── 📋 README.txt                   # Original setup instructions
├── 📚 docs/                        # Documentation
│   ├── 📖 Documentation.pdf        # Complete technical documentation (859 KB)
│   └── 🖼️ images/                  # Figures and diagrams
├── 🧪 examples/                    # Example cases
│   ├── 📁 simple_cube/             # Basic validation case
│   ├── 📁 tension_test/            # Mechanical loading example
│   └── 📁 validation/              # Verification cases
└── 🛠️ tools/                       # Utilities
    └── 📊 post_processing/          # Analysis scripts
```

## 🔬 Physics and Theory

### Multiphysics Coupling Mechanisms

#### 1. **Mechanical → Diffusion Coupling**
- **Tensile stress** (σh > 0): Enhances diffusion by opening material microstructure
- **Compressive stress** (σh < 0): Reduces diffusion by closing pore networks
- **Stress gradients**: Create preferential diffusion pathways

#### 2. **Stress-Diffusion Interaction**
The coupling term `(Vh/RT)∇·(σh∇c)` represents:
- Pressure-driven moisture transport
- Stress-concentration effects
- Hydrostatic stress influence on chemical potential

### Element Implementation

- **Element Type**: 20-node quadratic hexahedral (C3D20-like)
- **Integration**: 8-point Gauss quadrature
- **DOF**: 80 total (3 displacement + 1 concentration per node)
- **Time Integration**: Backward Euler for diffusion
- **Convergence**: Monolithic Newton-Raphson

## 🏭 Applications

### Aerospace Industry
- **Composite laminates**: Carbon fiber/epoxy moisture absorption
- **Environmental certification**: Hot/wet qualification testing
- **Stress concentration**: Bolt holes and cutouts in humid environments

### Marine Engineering
- **Underwater structures**: Pressure vessel moisture penetration
- **Adhesive joints**: Combined mechanical and environmental loading
- **Fatigue analysis**: Moisture-assisted crack growth

### Electronics Packaging
- **IC encapsulation**: Moisture diffusion under thermal cycling
- **Reliability testing**: Combined thermal-mechanical-moisture effects
- **Delamination prediction**: Interface failure mechanisms

### Infrastructure
- **Bridge components**: Moisture penetration under traffic loads
- **Protective coatings**: Mechanical damage accelerating moisture ingress
- **Durability assessment**: Long-term environmental degradation

## 📊 Material Properties

### Typical Values for Epoxy Materials

| Property | Symbol | Value | Unit |
|----------|--------|-------|------|
| Young's Modulus | E | 3000-4000 | MPa |
| Poisson's Ratio | ν | 0.35-0.40 | - |
| Diffusion Coefficient | D | 1×10⁻³ - 1×10⁻² | mm²/s |
| Molar Volume (Water) | Vh | 8000 | mm³/mol |
| Temperature | T | 300 | K |

### Property Variations
- **Temperature dependence**: D = D₀ exp(-Q/RT)
- **Concentration dependence**: D = D₀(1 + βc)
- **Stress dependence**: Implemented through coupling term

## 🔧 Advanced Usage

### Custom Material Properties

```fortran
*UEL PROPERTY, ELSET=MATERIAL_A
! E [MPa], ν [-], D_step1 [mm²/s], D_step2 [mm²/s]
3500., 0.37, 0.008, 0.012
```

### Boundary Conditions

```fortran
! Moisture flux boundary condition
*DFLUX
NodeSet_Surface, 11, 1.0E-6

! Concentration boundary condition  
*BOUNDARY
NodeSet_Exposed, 11, 11, 0.1
```

### Multi-Step Analysis

```fortran
*STEP, NAME=Drying
*COUPLED TEMPERATURE-DISPLACEMENT, STEADY STATE
! Use D = props(3)

*STEP, NAME=Loading  
*COUPLED TEMPERATURE-DISPLACEMENT
! Use D = props(4)
```

## 🧪 Validation and Testing

### Analytical Benchmarks
- [x] Pure diffusion (σ = 0) vs. Fick's law
- [x] Stress-free moisture uptake
- [x] Simple tension with diffusion
- [ ] Complex loading scenarios

### Experimental Validation
- [x] Moisture uptake curves
- [x] Stress-strain response
- [ ] Coupled behavior validation

## 📈 Results and Visualization

### Output Variables

| Variable | Description | Access |
|----------|-------------|---------|
| **U1, U2, U3** | Displacements | Standard ABAQUS output |
| **NT** | Concentration (as temperature) | Field output |
| **SDV1-SDV6** | Stress components | State variables |
| **SDV7-SDV12** | Strain components | State variables |
| **SDV14** | Hydrostatic stress | State variables |

### Post-Processing Example

```python
# Python script for ODB processing
from abaqus import *
from abaqusConstants import *

odb = openOdb('Analysis.odb')
step = odb.steps['Loading']
frame = step.frames[-1]

# Extract concentration field
concentration = frame.fieldOutputs['NT']
# Extract stress components  
stress_xx = frame.fieldOutputs['SDV1']
hydrostatic = frame.fieldOutputs['SDV14']
```

## 🤝 Contributing

Contributions are welcome! Please feel free to:

1. 🐛 **Report bugs** via [Issues](https://github.com/BBahtiri/ABAQUS-Multiphysics-Diffusion-UEL/issues)
2. 💡 **Suggest features** or improvements
3. 🔧 **Submit pull requests** with enhancements
4. 📖 **Improve documentation**
5. 🧪 **Add validation cases**

### Development Guidelines
- Follow Fortran 90+ standards
- Include test cases for new features
- Update documentation accordingly
- Maintain backward compatibility

## 📜 Citation

If you use this code in your research, please cite:

```bibtex
@software{bahtiri2025abaqus,
  title={ABAQUS Multiphysics Diffusion UEL: Coupled Stress-Diffusion Analysis},
  author={Bahtiri, Betim},
  year={2025},
  url={https://github.com/BBahtiri/ABAQUS-Multiphysics-Diffusion-UEL},
  note={User Element for coupled hydro-mechanical analysis in ABAQUS}
}
```

## 📄 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## 📞 Contact and Support

- **GitHub**: [@BBahtiri](https://github.com/BBahtiri)
- **Issues**: [GitHub Issues](https://github.com/BBahtiri/ABAQUS-Multiphysics-Diffusion-UEL/issues)
- **Discussions**: [GitHub Discussions](https://github.com/BBahtiri/ABAQUS-Multiphysics-Diffusion-UEL/discussions)

## 🔗 Related Work

- [ABAQUS User Subroutines](https://help.3ds.com/HelpProductsDS.aspx?release=2023&product=fe-safe&english=DSDocumentation)
- [Multiphysics Finite Element Methods](https://link.springer.com/book/10.1007/978-3-030-85688-1)
- [Moisture Diffusion in Polymers](https://www.sciencedirect.com/topics/materials-science/moisture-diffusion)

## ⭐ Acknowledgments

- ABAQUS User Community for documentation and examples
- Research community in computational mechanics
- Open source contributors and users

---

**Made with ❤️ for the computational mechanics community**

*Star ⭐ this repository if you find it useful!*
