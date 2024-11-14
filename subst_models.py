import numpy as np
from scipy.optimize import minimize


def extract_base_frequencies(Q):
    """
    Estimate base frequencies from the rate matrix Q.

    Parameters:
    - Q: 4x4 rate matrix (numpy array) with rows/columns in the order [A, C, G, T].

    Returns:
    - freqs: Estimated base frequencies as a 4-element list.
    """
    # Calculate base frequencies by ensuring rows sum to zero
    # Then normalize the diagonal entries as proportions
    row_sums = -np.diag(Q)  # Base frequencies should be proportional to diagonal elements
    return row_sums / row_sums.sum()

def gtr_loss(params, Q, freqs):
    """
    Loss function for fitting GTR rates to a general Q matrix.

    Parameters:
    - params: Array of six GTR rate parameters.
    - Q: Original 4x4 rate matrix.
    - freqs: Estimated equilibrium base frequencies.

    Returns:
    - loss: Sum of squared differences between the elements of Q and the GTR model approximation.
    """
    # Unpack the six GTR rates
    r_AC, r_AG, r_AT, r_CG, r_CT, r_GT = params
    # Construct the GTR model matrix based on the rates and frequencies
    GTR = np.array([
        [-freqs[1] * r_AC - freqs[2] * r_AG - freqs[3] * r_AT, r_AC * freqs[1], r_AG * freqs[2], r_AT * freqs[3]],
        [r_AC * freqs[0], -freqs[0] * r_AC - freqs[2] * r_CG - freqs[3] * r_CT, r_CG * freqs[2], r_CT * freqs[3]],
        [r_AG * freqs[0], r_CG * freqs[1], -freqs[0] * r_AG - freqs[1] * r_CG - freqs[3] * r_GT, r_GT * freqs[3]],
        [r_AT * freqs[0], r_CT * freqs[1], r_GT * freqs[2], -freqs[0] * r_AT - freqs[1] * r_CT - freqs[2] * r_GT]
    ])
    # Calculate sum of squared differences between GTR approximation and Q
    return np.sum((Q - GTR)**2)

def fit_gtr_model(Q):
    """
    Fit a GTR model to a given 4x4 rate matrix Q.

    Parameters:
    - Q: 4x4 rate matrix (numpy array) with rows/columns in the order [A, C, G, T].

    Returns:
    - rates: List of fitted GTR rates [r_AC, r_AG, r_AT, r_CG, r_CT, r_GT].
    - freqs: List of equilibrium base frequencies [pi_A, pi_C, pi_G, pi_T].
    """
    freqs = extract_base_frequencies(Q)
    # Initial guess for the rates (arbitrary starting values)
    initial_rates = np.ones(6)
    # Minimize the loss function to fit GTR rates to Q
    result = minimize(gtr_loss, initial_rates, args=(Q, freqs), method='L-BFGS-B', bounds=[(0, None)]*6)
    return result.x, freqs


def save_gtr_model(rates, freqs, output_file):
    """
    Save the fitted GTR model to a file in a format compatible with IQ-TREE.

    Parameters:
    - rates: List of GTR rates [r_AC, r_AG, r_AT, r_CG, r_CT, r_GT].
    - freqs: List of equilibrium base frequencies [pi_A, pi_C, pi_G, pi_T].
    - output_file: Path to the output .mod file.
    """
    with open(output_file, 'w') as f:
        f.write(f"GTR{{rates={','.join(map(str, rates))}; freqs={','.join(map(str, freqs))}}}\n")

# Example usage
Q = np.array([
    [-1.183773, 0.186562, 0.859435, 0.137776],
    [0.143636, -0.849065, 0.217282, 0.488148],
    [0.661688, 0.217282, -1.032053, 0.153083],
    [0.137775, 0.634030, 0.198833, -0.970637]
])
rates, freqs = fit_gtr_model(Q)
save_gtr_model(rates, freqs, "fitted_gtr.mod")
print("Fitted GTR model saved to fitted_gtr.mod")