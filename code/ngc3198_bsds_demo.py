# =====================================================
# SPARC data + BSDS rotation curve — version via Zenodo
# Galaxy: NGC 3198
# =====================================================

import os
import zipfile
import urllib.request
import ssl

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------------------------------
# 1. Download SPARC mass-model archive from Zenodo
# -----------------------------------------------------

zenodo_url = "https://zenodo.org/record/16284118/files/Rotmod_LTG.zip"
zip_filename = "Rotmod_LTG.zip"
extract_dir = "sparc_data"

# allow unverified HTTPS context if needed
ssl._create_default_https_context = ssl._create_unverified_context

if not os.path.exists(zip_filename):
    print("Downloading SPARC archive from Zenodo...")
    try:
        urllib.request.urlretrieve(zenodo_url, zip_filename)
        print("Download completed.")
    except Exception as e:
        raise RuntimeError(f"Failed to download SPARC archive: {e}")
else:
    print("SPARC archive already downloaded.")

if not os.path.exists(extract_dir):
    os.makedirs(extract_dir, exist_ok=True)
    with zipfile.ZipFile(zip_filename, 'r') as zf:
        zf.extractall(extract_dir)
    print("Archive extracted to", extract_dir)
else:
    print("SPARC data already extracted.")

# -----------------------------------------------------
# 2. Load NGC 3198 mass-model data
# -----------------------------------------------------

file_path = os.path.join(extract_dir, "NGC3198_rotmod.dat")
if not os.path.isfile(file_path):
    raise FileNotFoundError(f"{file_path} not found. Check zip archive contents.")

df = pd.read_csv(file_path, sep=r"\s+", comment="#", header=None)

# Assign columns (by position)
R     = df.iloc[:, 0].to_numpy()  # radial distance in kpc
Vobs  = df.iloc[:, 1].to_numpy()  # observed rotation speed km/s
errV  = df.iloc[:, 2].to_numpy()  # observational error km/s
Vgas  = df.iloc[:, 3].to_numpy()  # gas contribution km/s
Vdisk = df.iloc[:, 4].to_numpy()  # stellar disk contribution km/s
Vbul  = df.iloc[:, 5].to_numpy() if df.shape[1] > 5 else np.zeros_like(R)

print("Loaded NGC 3198 mass-model data: N =", len(R))

# -----------------------------------------------------
# 3. BSDS model definitions
# -----------------------------------------------------

G = 4.3009e-6  # gravitational constant in kpc (km/s)^2 / Msun
k0 = 0.03      # universal BSDS constant
Mb = 6.5e10    # Msun — baryonic mass approximation (adjust later if you refine)

def k_bsds(Mb):
    return k0 * (Mb / (5e10))**(-0.5)

def S_screen(r, rs=10.0, n=6):
    return 1.0 - np.exp(-(r/rs)**n)

def M_bsds_enclosed(r, Mb, rs=10.0, n=6):
    r = np.asarray(r)
    S = S_screen(r, rs, n)
    dr = np.diff(r)
    integral = np.concatenate([[0.0], np.cumsum(0.5 * (S[1:] + S[:-1]) * dr)])
    return k_bsds(Mb) * Mb * integral

def V_bsds(r, Mb, rs=10.0, n=6):
    r = np.asarray(r)
    M_bs = M_bsds_enclosed(r, Mb, rs, n)
    # avoid division by zero
    r_safe = np.where(r > 0, r, 1e-6)
    return np.sqrt(G * M_bs / r_safe)

# Compute rotation curve components
V_bary = np.sqrt(Vgas**2 + Vdisk**2 + Vbul**2)
V_bs   = V_bsds(R, Mb)
V_tot  = np.sqrt(V_bary**2 + V_bs**2)

# -----------------------------------------------------
# 4. Fit statistic and summary
# -----------------------------------------------------

chi2 = np.sum(((V_tot - Vobs) / errV)**2)
Npts = len(R)

print(f"\nNGC 3198 — BSDS Fit Results")
print(f"Data points: {Npts}")
print(f"Baryonic mass used (Mb): {Mb:.2e} Msun")
print(f"k(Mb) = {k_bsds(Mb):.4f}")
print(f"chi² = {chi2:.1f}")

# -----------------------------------------------------
# 5. Plot + Save figure
# -----------------------------------------------------

plt.figure(figsize=(9,6))
plt.errorbar(R, Vobs, yerr=errV, fmt='ko', ms=4, capsize=3, label='Observed')
plt.plot(R, V_bary, lw=2, label='Baryons only')
plt.plot(R, V_bs,   lw=2, label='BSDS only')
plt.plot(R, V_tot,  lw=2, label='Baryons + BSDS')
plt.xlabel("Radius (kpc)")
plt.ylabel("Velocity (km/s)")
plt.title("NGC 3198 — SPARC Data + BSDS Geometric Dark Sector")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()

output_fn = "NGC3198_BSDS_SPARC.png"
plt.savefig(output_fn, dpi=300)
plt.show()

print(f"\n✅ Figure saved as {output_fn}")
