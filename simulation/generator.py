# simulation/generator.py

import numpy as np
import gzip

class SyntheticDataGenerator:
    """
    Generates synthetic sequencing data with controlled fragment length 
    distribution and length-dependent PCR bias.
    """
    def __init__(self, N, depth_c, length_dist_func, lambda_func):
        """
        N: Number of unique fragments in the library.
        depth_c: Desired sequencing depth (M/N).
        length_dist_func: Function that returns an array of N lengths.
        lambda_func: Function that returns lambda for a given length.
        """
        self.N = N
        self.depth_c = depth_c
        self.length_dist_func = length_dist_func
        self.lambda_func = lambda_func
        self.fragments = []

    def generate_library(self):
        """Stage 1: Generate N unique fragments with lengths and normalized lambdas."""
        print(f"Generating library with {self.N} unique fragments...")
        lengths = self.length_dist_func(self.N)
        raw_lambdas = np.array([self.lambda_func(l) for l in lengths])
        
        # Normalize lambdas so their sum equals N
        sum_lambdas = np.sum(raw_lambdas)
        normalized_lambdas = raw_lambdas * (self.N / sum_lambdas)
        
        self.fragments = list(zip(lengths, normalized_lambdas))
    
    def simulate_sequencing(self, output_path):
        """Stage 2: Simulate sequencing using Poisson distribution."""
        if not self.fragments:
            self.generate_library()
            
        print(f"Simulating sequencing with target depth c={self.depth_c}...")
        
        with gzip.open(output_path, 'wt') as f:
            total_reads = 0
            unique_observed = 0
            lines_printed = 0
            
            for i, (length, lam) in enumerate(self.fragments):
                # Ki ~ Poiss(c * lam)
                ki = np.random.poisson(self.depth_c * lam)
                
                if ki > 0:
                    unique_observed += 1
                    total_reads += ki
                    
                    chrom = "chr1"
                    start = i * 1000 
                    end = start + int(length)
                    line = f"{chrom}\t{start}\t{end}\n"
                    
                    for _ in range(ki):
                        f.write(line)
                        if lines_printed < 10:
                            print(line.strip())
                            lines_printed += 1
            
            c_measured = total_reads / self.N
            print("-" * 30)
            print(f"Measured Depth (c_measured): {c_measured:.4f}")
            return total_reads, unique_observed