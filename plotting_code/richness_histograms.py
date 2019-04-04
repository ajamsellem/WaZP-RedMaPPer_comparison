# 2D Histogram with arbitrary binning
bins = [5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 25., 30., 50.,70.]
plt.hist2d(rich_u, rich, [bins, np.append(bins, 100)], cmap='GnBu')
cb = plt.colorbar()
cb.set_label('Counts in Bin', rotation = 270, labelpad = 13, fontsize = 14)
plt.xticks(bins, ['5', '', '', '', '', '10', '', '', '', '', '15', '', '', '', '', '20', '25', '30', '50','70'])
plt.yticks(np.append(bins, 100), ['5', '', '', '', '', '10', '', '', '', '', '15', '', '', '', '', '20', '25', '30', '50','70', '100'])
plt.xlabel("Number of Galaxies (WaZP)", fontsize = 15, labelpad = 10)
plt.ylabel("Richness (Redmapper)", fontsize = 15)

plt.savefig('./Richness_Ngal_compare_(Hist).png', bbox_inches="tight", dpi = 500)
plt.show()


# 2D Histogram with uniform binning and logarithmic color bar
binsx = np.linspace(5, 75, 10)
binsy = np.linspace(20,90,10)
plt.hist2d(rich_u_new, rich_new, [binsx, binsy], cmap='RdYlGn_r', norm=matplotlib.colors.LogNorm())
cb = plt.colorbar()
cb.set_label('Counts in Bin', rotation = 270, labelpad = 13, fontsize = 14)
plt.xlabel("Number of Galaxies (WaZP)", fontsize = 15, labelpad = 10)
plt.ylabel("Richness (Redmapper)", fontsize = 15)

plt.savefig('./Richness_Ngal_compare_(Hist)_Cut.png', bbox_inches="tight", dpi = 500)
plt.show()
