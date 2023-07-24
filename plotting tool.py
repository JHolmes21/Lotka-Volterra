import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# Load CSV data into a Pandas DataFrame
df = pd.read_csv('opper_onset_vary_r_redo.csv', skiprows=4)

# Filter data by ROW values 0, 1, and 5
df_filtered = df[df['GAMMA'].isin([1,1.5,2])]

# Set up plot
fig, ax = plt.subplots()

# Plot each GAMMA as a separate line
for GAMMA in df_filtered['GAMMA'].unique():
    data = df_filtered[df_filtered['GAMMA'] == GAMMA]
    ax.plot(data['VARIANCE'], data['Divergence Percentage'], label=f'$\gamma$ = {GAMMA}')

# opper r 
# ax.axvline(x=2, linestyle='--', color='C0')
# ax.axvline(x=2.5558283, linestyle='--', color='orange')
# ax.axvline(x=4.1678944, linestyle='--', color='green')

#MIFTY r
#ax.vlines(x=[3.18, 2.56,1.44], ymin=0, ymax=1, colors=['C0', 'orange', 'green'], linestyles='dashed')

# opper lg 
# ax.axvline(x=2.555828318070956, linestyle='--', color='grey')

#mifty lg
ax.vlines(x=[.88,0.68,.54], ymin=0, ymax=1, colors=['C0', 'orange','green'], linestyles='dashed')

# Add labels and legend
ax.set_xlabel(r'Variance of Interaction, $\sigma^2$', fontsize = 13)
ax.set_ylabel('Divergence Percentage', fontsize = 13)
ax.set_title(r'Numerical Onset of the $M \to\infty$ Transition, Varying $\gamma$', fontsize = 14.5)
ax.legend()

ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
plt.ylim(-0.005,1.005)
plt.xlim(0.2,1.5)
# Display plot
plt.show()