# Construct Fiduciary Data Set similar to Chang et. al.
rich_fid = []
z_fid = []
rich_u_fid = []
z_u_fid = []

for i, richness in enumerate(rich):
    if richness >= 20. and richness <= 100. and z_match[i] >= 0.2 and z_match[i] <= 0.55:
        rich_fid.append(richness)
        z_fid.append(z_match[i])
        
        rich_u_fid.append(rich_u[i])
        z_u_fid.append(z_u_match[i])
        

# Redmapper
binsx = np.linspace(0.25, 0.55, 15)
binsy = np.linspace(20, 95, 16)

rich_fid = np.asarray(rich_fid)
z_fid = np.asarray(z_fid)   

cmap = plt.cm.RdYlGn_r
cmap.set_bad(color='white')

plt.hist2d(z_fid, rich_fid, [binsx, binsy], cmap=cmap)
cb = plt.colorbar()
cb.set_label('Counts in Bin', rotation = 270, labelpad = 13, fontsize = 14)

plt.xlabel("Redshift (Redmapper)", fontsize = 15, labelpad = 10)
plt.ylabel("Richness (Redmapper)", fontsize = 15)

plt.savefig('./Redshift_Richness_compare_(RM).png', bbox_inches="tight", dpi = 500)
plt.show()


# WaZP
plt.hist2d(z_u_fid, rich_u_fid, cmap='RdYlGn_r')
cb = plt.colorbar()
cb.set_label('Counts in Bin', rotation = 270, labelpad = 13, fontsize = 14)
plt.xlabel("Redshift (WaZP)", fontsize = 15, labelpad = 10)
plt.ylabel("Number of Galaxies (WaZP)", fontsize = 15)

plt.savefig('./Redshift_Ngal_compare_(WAZP).png', bbox_inches="tight", dpi = 500)
plt.show()
