import matplotlib.pyplot as plt

fig, ax = plt.subplots()
output['pt1'].plot1d(ax=ax, overlay='dataset')
ax.set_yscale('log')
ax.set_ylim(1, None)
