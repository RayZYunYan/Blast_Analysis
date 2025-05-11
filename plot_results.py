import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("blast_summary.csv")

# Set style
sns.set(style="whitegrid")

# Plot 1: length vs identity (line plot)
plt.figure(figsize=(10, 6))
sns.lineplot(data=df, x="length", y="identity_percent", hue="algorithm", marker="o", ci="sd")
plt.title("Sequence Length vs Identity %")
plt.xlabel("Sequence Length (bp)")
plt.ylabel("Identity (%)")
plt.savefig("plot_length_vs_identity.png")
plt.close()

# Plot 2: average identity by algorithm (bar plot)
plt.figure(figsize=(8, 5))
sns.barplot(data=df, x="algorithm", y="identity_percent", estimator="mean", ci="sd")
plt.title("Average Identity % by Algorithm")
plt.ylabel("Mean Identity (%)")
plt.savefig("plot_algorithm_identity_avg.png")
plt.close()

# Plot 3: grouped bar plot by length and algorithm
plt.figure(figsize=(10, 6))
sns.barplot(data=df, x="length", y="identity_percent", hue="algorithm", ci="sd")
plt.title("Identity % by Length and Algorithm")
plt.ylabel("Identity (%)")
plt.savefig("plot_grouped_bar_length_algo.png")
plt.close()

print("✅ All plots saved: PNGs are in your folder.")
