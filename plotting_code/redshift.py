# Plot all redshift differences as a function of ra and dec
%matplotlib inline
import matplotlib.pyplot as plt

plt.scatter(ra_plot, dec_plot, z_diff, c= z_diff, label="Redmapper", marker = 'o')
plt.xlabel("Right Ascension ($\\degree$)")
plt.ylabel("Declination ($\\degree$)")
plt.title('Galaxy Catalogue Difference in Redshift')

cbar = plt.colorbar()
cbar.set_ticks([5,10,15,20,25,30,35])
cbar.set_ticklabels(['0.05','0.10','0.15','0.20','0.25','0.30','0.35'])
plt.savefig('./Redshift_Compare_1.png')
plt.show()


# Plot large redshift differences as a function of ra and dec
z_diff_big = np.asarray(z_diff_big)

plt.scatter(ra_big, dec_big, z_diff_big, c = z_diff_big,label="Redmapper", marker = 'o')
plt.xlabel("Right Ascension ($\\degree$)", fontsize = 13)
plt.ylabel("Declination ($\\degree$)", fontsize = 13, labelpad = 0)
plt.title('Large Redshift Differences Comparison (Sky Location)', fontsize = 13)

cbar = plt.colorbar()
cbar.set_ticks([15,20,25,30,35])
cbar.set_ticklabels(['0.15','0.20','0.25','0.30','0.35'])
cbar.set_label('Redshift Discrepancy', rotation = 270, labelpad = 19, fontsize = 13)

plt.savefig('./Large_Diff_Redshift_Compare_(RA_Dec).png', dpi = 500)
plt.show()


# Plot RM and WaZP redshift ranges
import matplotlib.pyplot as plt

height = 0.5

plt.hlines(3, np.min(z), np.max(z), 'r', label = 'Redmapper')
plt.vlines(np.min(z), 3 - height / 2., 3 + height / 2., 'r')
plt.vlines(np.max(z), 3 - height / 2., 3 + height / 2., 'r')

plt.hlines(1, np.min(z_u), np.max(z_u), 'g', label = 'WaZP')
plt.vlines(np.min(z_u), 1 - height / 2., 1 + height / 2., 'g')
plt.vlines(np.max(z_u), 1 - height / 2., 1 + height / 2., 'g')

plt.legend(fontsize = 12)
plt.xlabel("Redshift", fontsize = 12)

plt.xlim(-0.25,1.25)
plt.ylim(0,5)

plt.yticks([])
plt.text(0.14, 2.41, "Range: (" + str('{:.2e}'.format(np.min(z))) + ", " + str('{:.2e}'.format(np.max(z))) + ")", fontsize=12)
plt.text(0.025, .42, "Range: (" + str('{:.2e}'.format(np.min(z_u))) + ", " + str('{:.2e}'.format(np.max(z_u))) + ")", fontsize=12)

plt.savefig('./Redshift_Distribution.png', dpi = 600)
plt.show


# Graphic showing at what redshifts large differences in redshift occur
z_diff_big = np.asarray(z_diff_big)

line_x = np.linspace(0.33, 0.77, 100)
line_y = line_x

plt.scatter(z_big, z_u_big, z_diff_big, c = z_diff_big,label="Cluster", marker = 'o')
plt.plot(line_x, line_y, 'dimgray', label = 'Line of Data \nSet Agreement')
plt.xlabel("Redshift (Redmapper)", fontsize = 14)
plt.ylabel("Redshift (WaZP)", fontsize = 14)
plt.title('Large Redshift Differences Comparison (By Redshift)', fontsize = 14)
plt.xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

cbar = plt.colorbar()
cbar.set_ticks([15,20,25,30,35])
cbar.set_ticklabels(['0.15','0.20','0.25','0.30','0.35'])
cbar.set_label('Redshift Discrepancy', rotation = 270, labelpad = 19, fontsize = 14)

plt.legend(loc = 'lower right')

plt.savefig('./Large_Diff_Redshift_Compare_(Redshift).png', dpi = 500)
plt.show()
