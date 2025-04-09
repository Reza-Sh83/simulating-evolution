# Simulating Evolution by Mutation, Selection, and Genetic Drift

This project models the core mechanisms of biological evolution using Python, focusing on allele frequency changes over time in haploid and diploid populations. It simulates **natural selection**, **mutation**, **genetic drift**, and **ploidy effects** through mathematical models and visualizations.

## Overview

Through this simulation, you will explore:
- How natural selection and fitness shape allele frequencies.
- The role of genetic drift in small vs. large populations.
- Differences between haploid and diploid inheritance.
- Scenarios like heterozygote advantage/disadvantage.
- The interaction of mutation and selection.

## Key Features

- **Haploid and Diploid Selection Models**: Simulate allele frequency changes under different fitness regimes and ploidy levels.
- **Genetic Drift Simulation**: Model drift in finite populations and visualize random fluctuations.
- **Drift + Selection**: Combine both processes to observe stochastic vs. deterministic effects.
- **Heterozygote Advantage/Disadvantage**: Demonstrate maintenance or elimination of harmful alleles.
- **Mutation-Selection Balance**: Model how mutations can sustain harmful alleles.
- **Interactive Visualizations**: Matplotlib-based plots show evolutionary dynamics over time.

## Installation

Install required packages via pip:

```bash
pip install numpy matplotlib
```

## Usage

1. **Run the Script:**
   ```bash
   python simulating_evolution.py
   ```

2. **Examine Outputs:**
   - Allele frequency plots across generations.
   - Subplots showing how selection and drift interact.
   - Color-mapped simulations demonstrating the effect of initial conditions.

3. **Adjust Parameters:**
   You can customize:
   - Initial allele frequencies (e.g., `p`, `q`, `f_A1`, `f_A2`)
   - Fitness values (e.g., `Wpp`, `Wpq`, `Wqq`, `Wp`, `Wq`)
   - Population size
   - Number of generations
   - Selection and mutation modes

## Educational Value

This project is ideal for:
- **Students of evolutionary biology and population genetics**.
- **Researchers** wanting to prototype selection/mutation scenarios.
- **Educators** looking to demonstrate genetic drift and ploidy effects.

## Sample Visuals

- **Haploid vs Diploid Selection:**
  ![Diploid vs Haploid](haploid_vs_diploid.png)

- **Genetic Drift in Different Pop Sizes:**
  ![Genetic Drift](genetic_drift.png)

- **Drift + Selection:**
  ![Drift Selection Interaction](drift_selection.png)

- **Heterozygote Advantage/Disadvantage:**
  ![Heterozygote Fitness](heterozygote_dynamics.png)

## License

This project is released under the [MIT License](LICENSE).

## Acknowledgements

- Developed with inspiration from evolutionary dynamics courses and computational biology.
- Educational foundations based on Hardy-Weinberg models and modern synthesis theory.

---

Feel free to fork this repository and modify the simulation models for research or educational purposes!

