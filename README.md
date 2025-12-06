# Brane-Structured Dark Sector (BSDS)

A purely geometric alternative to dark matter arising from the electric part of the 
five-dimensional Weyl tensor projected onto a soft-wall brane. This repository contains 
the LaTeX manuscripts, analysis code, rotation-curve data, and reproducible workflows 
corresponding to the BSDS research series.

---

## Repository Structure

| Folder | Description |
|--------|-------------|
| papers/ | LaTeX manuscripts for the BSDS paper series |
| papers/paper1/ | Completed Paper 1 + figures + final PDF |
| code/ | Python/Colab scripts used to generate rotation-curve fits |
| data/ | SPARC rotation curve tables and processed inputs |
| figures/ | Rotation curves, BTFR plots, residual histograms |
| notebooks/ | Reproducible Jupyter/Colab notebooks for fitting and plotting |

---

## Key Result (Paper 1)

- Fits **175 SPARC galaxies** using **one universal constant**: k₀ = 0.069 ± 0.002  
- RMS residual velocity scatter: **6.8 km/s**, reduced χ² ≈ **1.06**  
- Reproduces the **baryonic Tully–Fisher relation** with **zero intrinsic scatter**
- Requires **no particle dark matter** — geometry alone generates the effect  
- Predictive scaling: **k(Mᵦ) ∝ Mᵦ⁻¹ᐟ²**, derived from soft-wall graviton normalization  

---

## How to Reproduce the Fits (when code is uploaded)

1. Clone the repository:2. Open `notebooks/BSDS_SPARK_fit.ipynb` (or Colab link once added)
3. Install dependencies:4. Load SPARC data from `/data/`
5. Run rotation-curve fitting cells
6. Output plots will save to `/figures/`

---

## Paper Series Roadmap

| Paper | Focus | Status |
|-------|-------|--------|
| Paper 1 | Galactic rotation curves & BTFR from geometry | **Complete** |
| Paper 2 | Full derivation of k(Mᵦ) ∝ Mᵦ⁻¹ᐟ² scaling from soft-wall KK modes | Drafting |
| Paper 3 | Cluster lensing + Bullet Cluster reconstruction | Drafting |
| Paper 4 | Cosmological background, growth, and CMB signatures | Planned |
| Paper 5+ | Predictions, falsifiable tests, simulations | Planned |

---

## Citation

If referencing this work:

> Patterson, K. (2025), *A Purely Geometric Dark Sector from Warped Extra Dimensions: The Brane-Structured Dark Sector (BSDS)*  
> https://github.com/<YOUR-USERNAME>/BSDS

---

## License

Recommended license for academic sharing while retaining ownership:

**CC-BY-NC-4.0 (Creative Commons Attribution–NonCommercial 4.0)**  
(You can apply this later via "Add license" in GitHub)

---

## Contact

For collaboration, discussion, or endorsement inquiries:  
**KennePatterson@gmail.com**
