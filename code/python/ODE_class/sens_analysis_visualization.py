import pandas as pd
import matplotlib.pyplot as plt
# helper function for decimal from sci notation
def decimal_converter(list):
    return [f"{float(x):.4f}" for x in list]

x = [4,6,8,10,12,14,16]
k5_y = [5.172885e-01,2.869495e-01,2.707002e-01,2.776512e-01,2.792339e-01,2.789572e-01,2.791388e-01]
k3_y = [3.992311e-01,2.726747e-01,2.695417e-01,2.760486e-01,2.762013e-01,2.754814e-01,2.756130e-01]
e3_y = [1.996194e-01,1.247840e-01,1.198194e-01,1.225195e-01,1.225337e-01,1.223394e-01,1.223913e-01]
k1_y = [1.858128e-01,5.488212e-02,5.771790e-02,5.957540e-02,5.975806e-02,5.970557e-02,5.973277e-02]
c1_y = [9.481663e-02,5.700392e-02,5.479694e-02,5.490583e-02,5.541765e-02,5.502604e-02,5.505464e-02]
c2_y = [8.358298e-02,4.611036e-02,4.693576e-02,4.775756e-02,4.714883e-02,4.744855e-02,4.745079e-02]

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

axes = axes.flatten()

datasets = [
    (k5_y, "K5"),
    (k3_y, "K3"),
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
