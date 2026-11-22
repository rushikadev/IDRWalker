# IDRWalker Study: Monte Carlo Random Walks for Disordered Protein Regions

A study of **IDRWalker**, a random-walk tool for modeling intrinsically disordered regions (IDRs) in proteins, applied to Spe-18, a highly disordered *C. elegans* protein. Completed as a Spring 2025 Random Walks course project by Rushika Devineni, Jonah Oxman, and Sophia Sharma.

**Paper:** Chen, G., & Zhang, Z. (2024). *IDRWalker: A Random Walk Based Tool for Generating Intrinsically Disordered Regions in Large Protein Complexes.* ACS Omega, 9(29), 32059–32065. https://doi.org/10.1021/acsomega.4c04161

## What this project does

1. **Applies IDRWalker to Spe-18.** Spe-18 is over 70% disordered, with AlphaFold confidence below 50% in those regions. We removed low-confidence regions, used IDRWalker to fill the resulting gaps, and refined the output with GROMACS energy minimization (run outside these scripts). The report describes a model with no spatial or connectivity conflicts.
2. **Studies the Monte Carlo random walk behind the tool.** A MATLAB simulation of protein backbone random walks at 7 chain lengths (5–200 residues, 1,000 trials each) examines end-to-end distance distributions, following the paper's analysis.
3. **Assesses gaps in a PDB file.** A MATLAB tool detects missing-residue gaps and checks whether random walks of the matching length can span the distance between the anchor residues.

## Repository contents

| File | Description | Author |
|---|---|---|
| `IDRWalker_github.py` | IDRWalker, the paper authors' tool | Original authors ([link to original repo]) |
| `remove_IDRs.py` | Filters atoms below a confidence threshold from an AlphaFold PDB | This project |
| `Spe18_IDRWalker.py` | Runs IDRWalker on the filtered Spe-18 structure | This project |
| `study_monte_carlo_project.m` | Monte Carlo random-walk study across chain lengths | This project |
| `idrwalker.m` | Detects PDB gaps and tests whether random-walk chains can span them | This project |
| `Random_Walks_Project.pdf` | Full project report | This project |

## Key results (from the report)

- **Distributions:** short chains (10 residues) give skewed end-to-end distance distributions. Longer chains (100 residues) are well approximated by a normal distribution.
- **Scaling:** RMS end-to-end distance against √N was close to linear (y = 6.43x − 2.20), more so than mean distance against N (y = 0.35x + 16.64).
- **Gap analysis:** `idrwalker.m` analyzed 11 gaps in Spe-18, ranging from 1 to 105 residues.

## Requirements and usage

- **Python:** `numpy`, `tqdm`, and IDRWalker's dependencies (`kearsley`, `phi_psi.txt`, amino acid templates)
- **MATLAB:** Statistics Toolbox (`fitdist`) and Bioinformatics Toolbox (`pdbread`)
- **GROMACS:** for energy minimization (not included)

```
python remove_IDRs.py        # edit input_pdb, output_pdb, and threshold at the bottom of the file
python Spe18_IDRWalker.py    # requires spe18.fasta and the filtered PDB
```
In MATLAB, run `study_monte_carlo_project` for the chain-length study, and `idrwalker('<file>.pdb', <num_trials>)` for gap analysis.

## Limitations

- IDRWalker's spatial conflict check is limited, so its output needs energy minimization before further use.
- The simulations use idealized backbone geometry with uniformly sampled dihedral angles and no self-avoidance, so they model an ideal random coil, not a physically realistic IDR.

## Credit

IDRWalker was developed by Chen and Zhang. Cite their paper if you use it. [Add the original repo's license and link here.]
