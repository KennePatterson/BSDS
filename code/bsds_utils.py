# bsds_utils.py — Core BSDS rotation curve functions
# Author: K. Patterson (Brane-Structured Dark Sector)
# ===================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Physical constant
G = 4.3009e-6   # kpc (km/s)^2 / Msun

# Universal BSDS constant (adjust per global SPARC fit)
k0_default = 0.069


# -------------------------------------------------------
# Screening function
# -------------------------------------------------------
def S_screen(r, rs=10.0, n=6):
    return 1.0 - np.exp(-(r/rs)**n)


# -------------------------------------------------------
# BSDS scaling
# -------------------------------------------------------
def k_bsds(Mb, k0=k0_default):
    return k0 * (Mb / (5e10))**(-0.5)


# -------------------------------------------------------
# Enclosed BSDS mass from geometric mode
# -------------------------------------------------------
def M_bsds_enclosed(r, Mb, rs=10.0, n=6):
    r = np.asarray(r)
    S = S_screen(r, rs, n)
    dr = np.diff(r)
    I = np.concatenate([[0.0], np.cumsum(0.5*(S[1:]+S[:-1])*dr)])
    return k_bsds(Mb)*Mb*I


# -------------------------------------------------------
# BSDS rotation velocity
# -------------------------------------------------------
def V_bsds(r, Mb, rs=10.0, n=6):
    r = np.asarray(r)
    Menc = M_bsds_enclosed(r, Mb, rs, n)
    r_safe = np.where(r>0, r, 1e-6)
    return np.sqrt(G * Menc / r_safe)


# -------------------------------------------------------
# Main runner — loads SPARC file, fits, plots, saves fig
# -------------------------------------------------------
def run_bsds(filename, Mb, rs=10, n=6, save=True):
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"{filename} not found.")

    df = pd.read_csv(filename, sep=r"\s+", comment="#", header=None)

    R     = df.iloc[:,0].to_numpy()
    Vobs  = df.iloc[:,1].to_numpy()
    errV  = df.iloc[:,2].to_numpy()
    Vgas  = df.iloc[:,3].to_numpy()
    Vdisk = df.iloc[:,4].to_numpy()
    Vbul  = df.iloc[:,5].to_numpy() if df.shape[1]>5 else np.zeros_like(R)

    V_bary = np.sqrt(Vgas**2 + Vdisk**2 + Vbul**2)
    V_bs   = V_bsds(R, Mb, rs, n)
    V_tot  = np.sqrt(V_bary**2 + V_bs**2)

    chi2 = np.sum(((V_tot-Vobs)/errV)**2)

    # Plot
    plt.figure(figsize=(9,6))
    plt.errorbar(R, Vobs, yerr=errV, fmt='ko', ms=4, capsize=3,label='Observed')
    plt.plot(R,V_bary,label='Baryons only')
    plt.plot(R,V_bs,label='BSDS only')
    plt.plot(R,V_tot,label='Baryons + BSDS')
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Velocity (km/s)")
    plt.title(f"{filename} — BSDS Rotation Curve")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()

    outfile = filename.replacegg(".dat","_BSDS.png")
    if save:
        plt.savefig(outfile,dpi=300)
    plt.show()

    print("\n=== BSDS Fit Complete ===")
    print(f"Galaxy file: {filename}")
    print(f"Mb used: {Mb:.2e} Msun")
    print(f"k(Mb) = {k_bsds(Mb):.4f}")
    print(f"chi² = {chi2:.1f}")
    print(f"Saved figure: {outfile}\n")

    return outfile,chi2