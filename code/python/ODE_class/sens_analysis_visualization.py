import pandas as pd
import matplotlib.pyplot as plt
# helper function for decimal from sci notation
def decimal_converter(list):
    return [f"{float(x):.4f}" for x in list]

x = [4,6,8,10,12,14,16]
k5_y = [4.690061e-01,2.572608e-01,2.412119e-01,2.485901e-01,2.499403e-01,2.497006e-01,2.498575e-01]
k3_y = [4.227757e-01,2.847610e-01,2.817681e-01,2.874942e-01,2.876493e-01,2.869227e-01,2.870554e-01]
e3_y = [2.131046e-01,1.309609e-01,1.255694e-01,1.284154e-01,1.284127e-01,1.282216e-01,1.282703e-01]
k1_y = [1.925263e-01,5.580824e-02,5.883532e-02,6.064919e-02,6.083066e-02,6.078846e-02,6.080791e-02]
c1_y = [1.012035e-01,6.000662e-02,5.758015e-02,5.770420e-02,5.824292e-02,5.784382e-02,5.786813e-02]
c2_y = [8.781849e-02,4.806099e-02,4.861209e-02,4.948258e-02,4.884286e-02,4.915714e-02,4.916192e-02]

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
