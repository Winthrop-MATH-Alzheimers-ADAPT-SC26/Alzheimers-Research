import pandas as pd
import matplotlib.pyplot as plt
# helper function for decimal from sci notation
def decimal_converter(list):
    return [f"{float(x):.4f}" for x in list]

x = [4,6,8,10,12,14,16,18,20]
k5_y = [2.102129e-01,2.660479e-01,2.596922e-01,2.436954e-01,2.429086e-01,2.429657e-01,2.431698e-01,2.432004e-01,2.431817e-01]
k3_y = [2.332754e-01,2.714512e-01,3.011183e-01,2.828300e-01,2.756811e-01,2.780218e-01,2.782584e-01,2.783015e-01,2.782781e-01]
e3_y = [1.079469e-01,1.211876e-01,1.314145e-01,1.253703e-01,1.243911e-01,1.242627e-01,1.243811e-01,1.243939e-01,1.243843e-01]
k1_y = [1.925263e-01,5.580824e-02,5.883532e-02,6.064919e-02,6.083066e-02,6.078846e-02,6.080791e-02,6.020412e-02,6.020095e-02]
c1_y = [6.262175e-02,7.975393e-02,8.832893e-02,8.043565e-02,8.073674e-02,8.050334e-02,8.057575e-02,8.058228e-02,8.047175e-02]
c2_y = [4.688306e-02,6.279251e-02,6.221974e-02,5.912209e-02,5.820809e-02,5.818323e-02,5.824234e-02,5.824374e-02,5.823958e-02]

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

axes = axes.flatten()

datasets = [
    (k3_y, "K3"),
    (k5_y, "K5"),
    (e3_y, "E3"),
    (k1_y, "K1"),
    (c1_y, "C1"),
    (c2_y, "C2")
]

for i, (y_data, title) in enumerate(datasets):
    axes[i].plot(x, y_data, marker="o", markeredgecolor="black" ,linestyle = "--", color="red") # Use raw floats for numerical plotting
    axes[i].set_title(title)
    axes[i].set_xlabel("$log_2(N)$")
    axes[i].set_ylabel("Total Sensitivity")
    axes[i].grid(True, which="both", ls="--")

plt.tight_layout()
plt.show()
