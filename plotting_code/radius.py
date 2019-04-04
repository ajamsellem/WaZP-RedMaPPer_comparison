# Initial info from data set comparison
max_rad = max(max(R), max(R_u))
min_rad = min(min(R), min(R_u))

print("Number of Matched Galxies: " + str(len(r_diff)))
print("Max Radius Measurement of both catalogues: " + str(max_rad))
print("Min Radius Measurement of both catalogues: " + str(min_rad))

count_1 = 0
count_2 = 0
count_3 = 0
count_4 = 0

# Count how many cluster radii are off by more than 'x' and print info
for value in r_diff:
    if value >= 0.1:
        count_1 = count_1 + 1
    if value >= 0.2:
        count_2 = count_2 + 1
    if value >= 0.3:
        count_3 = count_3 + 1
    if value >= 0.4:
        count_4 = count_4 + 1

print("Number of Matched Radii off by 0.1 (MPC?): " + str(count_1))
print("Number of Matched Radii off by 0.2 (MPC?): " + str(count_2))
print("Number of Matched Radii off by 0.3 (MPC?): " + str(count_3))
print("Number of Matched Radii off by 0.4 (MPC?): " + str(count_4))
print("Conclusion: These measurements of radius are not measuring the "
      "same thing, or we have pretty large discrepancies "
      "between these two data sets.")
      

# Plot the radius of each cluster versus the redshift difference between data sets distinguishing between
# big (green) and small (purple) differences in redshift
%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
ax1.scatter(r, z_diff, color = 'purple', marker = '.', zorder = 0, label = "Small Z Difference")
ax1.scatter(r_big, z_diff_big, color = "green", marker = '.', zorder = 1, label = "Big Z Difference")
ax2.scatter(r_u, z_diff, color = 'purple', marker = '.', zorder = 0)
ax2.scatter(r_u_big, z_diff_big, color = "green", marker = '.', zorder = 1)

ax1.set_xlabel("Redmapper Radius (Mpc)")
ax2.set_xlabel("WaZP Radius (Mpc)")
ax1.set_ylabel("Redshift Difference")

ax1.set_yticks([5,10,15,20,25,30,35])
ax2.set_yticks([5,10,15,20,25,30,35])
ax1.set_yticklabels(['0.05','0.10','0.15','0.20','0.25','0.30','0.35'])
ax2.set_yticklabels(['0.05','0.10','0.15','0.20','0.25','0.30','0.35'])

title = 'Z Difference \nv. Radius'
font = FontProperties()
font.set_style('italic')
font.set_weight('bold')
ax1.text(0.08, 1.03, title, transform=ax1.transAxes, fontproperties=font, fontsize=12)

f.legend(ncol=2)

plt.savefig('./Redshift_Difference_v_Radius.png')
plt.show()


# Discretized histogram of cluster radii for RM and WaZP
# Arbitrary discretized binning
bins = [.2, .3, .4, .5, .6, .7, .8,  .9, 1.]

# Plot Histogram of Radii Distribution for both data sets
fig, ax = plt.subplots()
colors = ['blue', 'orange']
plt.hist([R, R_u], bins, histtype='bar', color=colors, label= ['Redmapper', 'WaZP'])
plt.legend(fontsize = 14)
plt.xlabel('Radius (Mpc)', fontsize = 14)
plt.ylabel('Number of Clusters', fontsize = 14)

ax.set_xticks([0.2, 0.3, 0.4, .5, 0.6, .7, 0.8, .9, 1.])
ax.set_xticklabels(['0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'])

ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.savefig('./Radii_Distribution.png', dpi = 600)

plt.show()

# Show that there are some radii values that lie beyond the x-axis bounds in this plot
print('WaZP radii below 0.2: ' +  str(R_u[np.argsort(R_u)[0:15]]))
print('\nWaZP radii above 1.0: ' + str(R_u[np.argsort(R_u)[-10:]]))
print('\nRM radii above 1.0: ' + str(R[np.argsort(R)[-50:]]))
