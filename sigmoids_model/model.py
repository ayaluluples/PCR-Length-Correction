# sigmoids_model/model.py

import numpy as np
import gzip
import math
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.optimize import curve_fit, brentq

class SigmoidsModel:
    """
    A model to estimate and correct PCR length bias in cfDNA sequencing data
    using a double sigmoid (band) function.
    """
    
    def __init__(self, input_path, bin_size=5, max_len=600):
        """
        Initialize and fit the model directly from a tagAlign file.
        
        Args:
            input_path (str): Path to the .tagAlign.gz file.
            bin_size (int): Size of the length bins in base pairs.
            max_len (int): Maximum fragment length to consider.
        """
        self.bin_size = bin_size
        self.max_len = max_len
        self.popt = None
        self.is_fitted = False
        
        # Data storage for plotting and analysis
        self.raw_bin_starts = []
        self.raw_lambdas = []
        self.M_per_bin = [] # Total reads observed (M)
        self.Mu_per_bin = [] # Unique molecules observed (Mu)

        # Execute the pipeline
        print(f"Processing file: {input_path}...")
        bin_data = self._create_bins(input_path)
        
        print("Estimating lambdas...")
        self._estimate_all_lambdas(bin_data)
        
        print("Fitting model...")
        self._fit_model()
        print("Done!")

    @staticmethod
    def sigmoid(t):
        """Standard logistic sigmoid function."""
        return 1.0 / (1.0 + np.exp(-t))

    @staticmethod
    def band_sigmoid(x, A, x1, k1, x2, k2, C):
        """Double sigmoid function defining a 'band' of efficiency."""
        return C + A * SigmoidsModel.sigmoid((x - x1) / k1) * SigmoidsModel.sigmoid((x2 - x) / k2)

    def _create_bins(self, input_path):
        """Read tagAlign file and group molecules into length bins."""
        bins = defaultdict(lambda: {'rows': 0, 'unique_molecules': set()})
        with gzip.open(input_path, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                length = end - start
                if 0 <= length < self.max_len:
                    bin_start = (length // self.bin_size) * self.bin_size
                    bins[bin_start]['rows'] += 1
                    bins[bin_start]['unique_molecules'].add((chrom, start, end))
        return bins

    def _estimate_lambda(self, M, Mu):
        """Estimate lambda for a single bin using the Poisson equation."""
        if Mu <= 0: return float("inf")
        if Mu >= M: return 0.0 if Mu == M else float("nan")
        g = lambda lam: M * (1.0 - math.exp(-lam)) / lam - Mu
        lo, hi = 1e-12, 1.0
        while g(hi) > 0:
            hi *= 2.0
            if hi > 1e6: return float("nan")
        return brentq(g, lo, hi)

    def _estimate_all_lambdas(self, bin_data):
        """Calculate initial lambda estimates for all available bins."""
        for b_start in sorted(bin_data.keys()):
            M = bin_data[b_start]['rows']
            Mu = len(bin_data[b_start]['unique_molecules'])
            lam = self._estimate_lambda(M, Mu)
            
            if math.isfinite(lam):
                self.raw_bin_starts.append(b_start)
                self.raw_lambdas.append(lam)
                self.M_per_bin.append(M)
                self.Mu_per_bin.append(Mu)

    def _fit_model(self):
        """Fit the continuous double-sigmoid curve to the raw lambda estimates."""
        x = np.array(self.raw_bin_starts, dtype=float)
        y = np.array(self.raw_lambdas, dtype=float)
        p0 = [max(y) - min(y), 150.0, 20.0, 430.0, 30.0, min(y)]
        bounds = ([0, 0, 1e-3, 0, 1e-3, -np.inf], [np.inf, 1000, 500, 1000, 500, np.inf])
        self.popt, _ = curve_fit(self.band_sigmoid, x, y, p0=p0, bounds=bounds, maxfev=20000)
        self.is_fitted = True

    def predict_lambda(self, length):
        """Predict lambda for any given fragment length using the fitted model."""
        return self.band_sigmoid(length, *self.popt)

    def plot_lambda_vs_length(self, show_fit=True):
        """Plot raw lambda estimates vs. fragment length."""
        plt.figure(figsize=(10, 6))
        plt.scatter(self.raw_bin_starts, self.raw_lambdas, alpha=0.5, label="Raw Estimates (Brentq)", color='C0')
        
        if show_fit and self.is_fitted:
            x_range = np.linspace(0, self.max_len, 500)
            y_fit = self.predict_lambda(x_range)
            plt.plot(x_range, y_fit, color='red', linewidth=2, label="Sigmoid Fit")
        
        plt.xlabel("Fragment Length (bp)")
        plt.ylabel("Estimated Lambda (λ)")
        plt.title(f"λ Estimation vs. Length {'(with Fit)' if show_fit else ''}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

    def plot_m_and_n(self):
        """Plot observed read counts (M) vs. estimated original molecules (N)."""
        M = np.array(self.M_per_bin)
        N_est = [m / self.predict_lambda(b) if self.predict_lambda(b) > 0 else 0 
                 for m, b in zip(self.M_per_bin, self.raw_bin_starts)]

        plt.figure(figsize=(10, 6))
        plt.bar(self.raw_bin_starts, M, width=self.bin_size, alpha=0.4, label="Observed Reads (M)", color='gray')
        plt.step(self.raw_bin_starts, N_est, where='mid', label="Estimated Original Molecules (N)", color='green', linewidth=2)
        
        plt.xlabel("Fragment Length (bp)")
        plt.ylabel("Count")
        plt.title("Read Count (M) vs. Estimated Original Library (N)")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
    
    def estimate_total_unique_molecules(self):
        """Calculate the total estimated size of the original library."""
        total_n = 0
        for m, b_start in zip(self.M_per_bin, self.raw_bin_starts):
            lam = self.predict_lambda(b_start)
            if lam > 0:
                total_n += m / lam
        return total_n

    def plot_normalized_density(self):
        """Plot relative density of observed reads vs. corrected molecules."""
        M = np.array(self.M_per_bin)
        N_values = np.array([m / self.predict_lambda(b) if self.predict_lambda(b) > 0 else 0 
                            for m, b in zip(self.M_per_bin, self.raw_bin_starts)])
        
        m_density = M / np.sum(M)
        n_density = N_values / np.sum(N_values)

        plt.figure(figsize=(10, 6))
        plt.fill_between(self.raw_bin_starts, m_density, step="mid", alpha=0.3, 
                         label="Observed Reads Density (M)", color='gray')
        plt.step(self.raw_bin_starts, n_density, where='mid', 
                 label="Corrected Molecule Density (N)", color='green', linewidth=2.5)
        
        plt.xlabel("Fragment Length (bp)")
        plt.ylabel("Density (Relative Abundance)")
        plt.title("Corrected cfDNA Length Distribution (Density)")
        plt.legend()
        plt.grid(True, alpha=0.2)
        
        total_n = self.estimate_total_unique_molecules()
        plt.text(0.65, 0.5, f"Est. Total Library Size:\n{total_n:,.0f} molecules", 
                 transform=plt.gca().transAxes, fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
        plt.show()

    def plot_validation_experiment(self, gmm_params=None):
        """Plot comparison between Ground Truth, Naive, and Corrected densities."""
        if not self.is_fitted:
            print("Error: Model must be fitted first.")
            return

        bin_starts = np.array(self.raw_bin_starts)
        M = np.array(self.M_per_bin)
        N_est = np.array([m / self.predict_lambda(b) if self.predict_lambda(b) > 0 else 0 
                         for m, b in zip(M, bin_starts)])
        
        m_density = M / (np.sum(M) * self.bin_size)
        n_density = N_est / (np.sum(N_est) * self.bin_size)

        plt.figure(figsize=(12, 7))
        plt.fill_between(bin_starts, m_density, step="mid", alpha=0.2, 
                         label="Naive Estimate (Observed Reads M)", color='gray')
        plt.step(bin_starts, n_density, where='mid', 
                 label="SigmoidsModel Correction (Estimated N)", color='green', linewidth=2.5)

        if gmm_params:
            x_range = np.linspace(min(bin_starts), max(bin_starts), 500)
            true_gmm_y = np.zeros_like(x_range)
            for mu, sigma, w in gmm_params:
                term = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x_range - mu) / sigma)**2)
                true_gmm_y += w * term
            plt.plot(x_range, true_gmm_y, '--', label="Ground Truth (Original GMM)", color='red', linewidth=2)

        plt.xlabel("Fragment Length (bp)")
        plt.ylabel("Probability Density")
        plt.title("Model Validation: Correcting Bias to Recover Ground Truth")
        plt.legend()
        plt.grid(True, alpha=0.2)
        plt.show()