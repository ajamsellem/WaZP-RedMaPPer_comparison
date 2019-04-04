from astropy.io import fits
import numpy as np
%matplotlib inline
import matplotlib.pyplot as plt

# Read in Redmapper data
hdul = fits.open('./redmapper.fit')
data = hdul[1].data
cols = hdul[1].columns
print(cols.names)

# Coordinates and Data of RM clusters
ra = data.field('RA')
dec = data.field('DEC')
R = data.field('R_LAMBDA')
z = data.field('Z_LAMBDA')
rich = data.field('LAMBDA_CHISQ')
drich = data.field('LAMBDA_CHISQ_E')

# Read in the WaZP (Everything from WaZP is denoted with a '_u' at the end)
hdu2 = fits.open('./unpublished.fits')
data_u = hdu2[1].data
cols_u = hdu2[1].columns
print(cols_u.names)

# Coordinates and Data of WAZP clusters
ra_u = data_u.field('ra')
dec_u = data_u.field('dec')
R_u = data_u.field('RADIUS_ISO_MPC')
z_u = data_u.field('zp')
rich_u = data_u.field('NGals')
drich_u = data_u.field('ngals_err')

# Create an array of tuples of the form (ra, dec) for each data set
data = list(zip(ra, dec))
data_u = list(zip(ra_u, dec_u))
data = np.array(data)
data_u = np.array(data_u)

# Find where the the two data sets have equal (ra, dec) tuples and create an array of these tuples
aset = set([tuple(x) for x in data])
bset = set([tuple(x) for x in data_u])
ra_dec = np.array([x for x in aset & bset])

# Find indices where the RM tuples are equal to the ra_dec tuples
idx = []
for i,j in enumerate(data):
    if (j in ra_dec) == True:
        idx.append(i)

# Find indices where the WaZP tuples are equal to the ra_dec tuples
idx_u = []
for i,j in enumerate(data_u):
    if (j in ra_dec) == True:
        idx_u.append(i)
        
# Make lists of all the needed data at matched galaxies in the catalogues, using the above index lists
z_list = []
ra_list = []
dec_list = []
r_list = []
rich_list = []
drich_list = []

z_u_list = []
ra_u_list = []
dec_u_list = []
r_u_list = []
rich_u_list = []
drich_u_list = []


for i in idx:
    z_list.append(z[i])
    ra_list.append(ra[i])
    dec_list.append(dec[i])
    r_list.append(R[i])
    rich_list.append(rich[i])
    drich_list.append(drich[i])
    
for i in idx_u:
    z_u_list.append(z_u[i])
    ra_u_list.append(ra_u[i])
    dec_u_list.append(dec_u[i])
    r_u_list.append(R_u[i])
    rich_u_list.append(rich_u[i])
    drich_u_list.append(drich_u[i])
    
# Match all the data so at each index value, the information in each list is from the same cluster
ra_plot = []
dec_plot = []
z_diff = []
r_diff = []
r_diff_plot = []
r = []
r_u = []
rich = []
drich = []
rich_u = []
drich_u = []

for i,j in enumerate(ra_list):
    for k,l in enumerate(ra_u_list):
        if j == l and dec_list[i] == dec_u_list[k]:
            ra_plot.append(j)
            dec_plot.append(dec_list[i])
            z_diff.append((np.absolute(z_list[i] - z_u_list[k]))*100)
            r_diff_plot.append((np.absolute(r_list[i] - r_u_list[k]))*50)
            r_diff.append((np.absolute(r_list[i] - r_u_list[k])))
            r.append(r_list[i])
            r_u.append(r_u_list[k])
            rich.append(rich_list[i])
            drich.append(drich_list[i])
            rich_u.append(rich_u_list[k])
            drich_u.append(drich_u_list[k])

# Make a list of clusters with largest redshift differences between data sets
big_guys = []
for value in z_diff:
    # Take the 100 largest values in z_diff, which happens to be all values over a threshold of .1071 difference
    if value >= 10.71:
        big_guys.append(value)

# Make lists of data corresponding to the large redshift difference cluster
ra_big = []
dec_big =[]
z_diff_big = []
r_big = []
r_u_big = []

for big_val in big_guys:
    for i,diff in enumerate(z_diff):
        if big_val == diff:
            ra_big.append(ra_plot[i])
            dec_big.append(dec_plot[i])
            z_diff_big.append(diff)
            r_big.append(r[i])
            r_u_big.append(r_u[i])
