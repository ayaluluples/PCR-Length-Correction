# SigmoidsModel: PCR Length Bias Correction for cfDNA

This project provides a computational tool to estimate and correct fragment length biases induced by PCR amplification and sequencing processes, specifically designed for cell-free DNA (cfDNA) data.

## 🧬 The Problem
In cfDNA sequencing, PCR amplification and library preparation often exhibit significant bias based on fragment length. Since the length distribution of cfDNA is a critical feature for clinical diagnostics (e.g., cancer detection), these technical biases can distort biological signals and lead to inaccurate quantification of the original molecule population.

## 🛠️ The Solution
`SigmoidsModel` implements a mathematical approach to "un-bias" the data:

1.  **Efficiency Estimation**: It models the sequencing efficiency ($\lambda$) per length bin using a Poisson-based approach by comparing total reads ($M$) and unique molecules ($M_u$).
2.  **Double Sigmoid Fit**: It fits a continuous "band-pass" sigmoid function to these estimates to smooth out noise and capture the rising and falling efficiency gates.
3.  **Data Correction**: It applies a correction factor ($N \approx M/\lambda$) to recover the original molecule counts ($N$) from the observed read counts ($M$).

## 🧪 Model Validation (Example at $c=1.0$)
The model was validated using a **Synthetic Data Generator** simulating a 3-Gaussian Mixture Model (GMM) library. As shown below, the corrected distribution (green) effectively recovers the ground truth (red dashed) from the biased sequencing data.

![Model Validation](Figure_1.png)
*Figure 1: Comparison between Observed Reads (Gray), Corrected Molecules (Green), and the GMM Ground Truth (Red Dashed).*

## 📉 Sensitivity to Sequencing Depth ($c$)
To determine the model's reliability, we analyzed its performance across a range of sequencing depths ($c$). This analysis identifies the **Threshold Depth**—the minimum data volume required for a stable and effective correction.

| Sequencing Depth (c) | Naive MSE (Uncorrected) | Model MSE (Corrected) | Improvement % |
|:--------------------|:-----------------------:|:---------------------:|:-------------:|
| 0.05                | 9.65e-07                | 7.48e-06              | -675.09%      |
| 0.10                | 8.52e-07                | 2.53e-07              | 70.35%        |
| 0.20                | 8.04e-07                | 1.34e-07              | 83.36%        |
| 0.50                | 7.97e-07                | 9.21e-08              | 88.44%        |
| 1.00                | 7.96e-07                | 7.28e-08              | 90.86%        |

**Key Finding:** The correction becomes highly effective at $c \ge 0.1$. While very low depths ($c=0.05$) may lead to overfitting due to sampling noise, depths of $c \ge 0.2$ consistently eliminate over 80% of the length-dependent bias.

## 🚀 Future Roadmap
* **GC-Content Integration**: Extending the model to handle multi-variate bias correction, incorporating both fragment length and GC content.
* **Real-world Benchmarking**: Evaluating the model on clinical cfDNA datasets (cancer patients vs. healthy controls) to demonstrate its impact on fragmentomics-based cancer detection and diagnostics.

## 📂 Project Structure
* `sigmoids_model/`: Core package containing the `SigmoidsModel` class.
* `simulation/`: Tools for generating synthetic cfDNA data for testing and validation.
* `notebooks/`: Demonstration of the validation experiment and usage examples.

## 🚀 Getting Started
1. **Clone the repository**:
   ```bash
   git clone [https://github.com/YourUsername/SigmoidsModel.git](https://github.com/YourUsername/SigmoidsModel.git)

2. **Install dependencies**
   Ensure you have Python installed, then run:
   ```bash
   pip install -r requirements.txt

3. **Run the demo**
   Open the provided Jupyter Notebook to see the model in action:
   jupyter notebook notebooks/validation_demo.ipynb