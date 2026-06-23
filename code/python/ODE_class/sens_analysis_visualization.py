import pandas as pd
import matplotlib.pyplot as plt
 
RESULTS_CSV = "sobol_results.csv"
TOP_N = 6
 
df = pd.read_csv(RESULTS_CSV)
 
max_n = df['N'].max()
top_params = (
    df[df['N'] == max_n]
    .sort_values('Sobol Total (ST)', ascending=False)
    .head(TOP_N)['Params']
    .tolist()
)
 
n_cols = 3
n_rows = -(-len(top_params) // n_cols)  # ceil division, in case TOP_N != 6
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
axes = axes.flatten()
 
for i, param in enumerate(top_params):
    sub = df[df['Params'] == param].sort_values('N')
    axes[i].plot(
        sub['N'], sub['Sobol Total (ST)'],
        marker="o", markeredgecolor="black", linestyle="--", color="red"
    )
    axes[i].set_title(param)
    axes[i].set_xlabel("$log_2(N)$")
    axes[i].set_ylabel("Total Sensitivity")
    axes[i].grid(True, which="both", ls="--")
 
for j in range(len(top_params), len(axes)):
    axes[j].set_visible(False)
 
plt.tight_layout()
plt.savefig("sobol_sensitivity_plot.png", dpi=150)
plt.show()
