
import astropy.io.fits as pf
import numpy as np
from astropy import units as u
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import sys
import treecorr
import healpy as hp

# setting up config file and correlator for treecorr
config_file = 'default.params'
config = treecorr.config.read_config(config_file)

zmin = float(sys.argv[1])
zmax = float(sys.argv[2])
Nz = int(sys.argv[3])
Maglim1 = float(sys.argv[4])
Maglim2 = float(sys.argv[5])
lambmin = float(sys.argv[6])
lambmax = float(sys.argv[7])

clusters = pf.open(sys.argv[8])[1].data
randoms = pf.open(sys.argv[9])[1].data

galaxies = pf.open(sys.argv[10])[1].data
galaxy_randoms = pf.open(sys.argv[11])[1].data
jkid = int(sys.argv[12])
tot_area = float(sys.argv[13])
outfile = sys.argv[14]
nR = int(sys.argv[15])

# Measurement parameters
h = 0.7
Rmin = 0.1/h   #Mpc
Rmax = 20.0/h
lnrperp_bins = np.linspace(np.log(Rmin), np.log(Rmax), num = nR+1)
R_edge = np.exp(lnrperp_bins)
R_mid = np.sqrt(R_edge[:-1] * R_edge[1:])
bslop = 0.03
nside = 8

# read and mask JK sub-sample
JK = clusters['JK']
mask = (JK==jkid)*(clusters['Z']>=zmin)*(clusters['Z']<zmax)*(clusters['LAMBDA']>=lambmin)*(clusters['LAMBDA']<lambmax)
RA = clusters['RA'][mask]
DEC = clusters['DEC'][mask]
Z = clusters['Z'][mask]
LAMB = clusters['LAMBDA'][mask]
RA[RA>180] -= 360.

JK_ran = randoms['JK']
mask = (JK_ran==jkid)*(randoms['Z']>=zmin)*(randoms['Z']<zmax)*(randoms['LAMBDA']>=lambmin)*(randoms['LAMBDA']<lambmax)
RA_ran = randoms['RA'][mask]
DEC_ran = randoms['DEC'][mask]
Z_ran = randoms['Z'][mask]
W_ran = randoms['W'][mask]
RA_ran[RA_ran>180] -= 360.

# this is just to count the galaxies
JK_gal = galaxies['JK']
mask = (JK_gal==jkid)
mag = galaxies['MAG_AUTO_I'][mask]

ra_gal_all = galaxies['RA']
dec_gal_all = galaxies['DEC']
mag_gal_all = galaxies['MAG_AUTO_I']
theta = (90.0 - dec_gal_all)*np.pi/180.
phi = ra_gal_all*np.pi/180.
pix = hp.ang2pix(nside,theta,phi)
pixjk = hp.ang2pix(nside,((90.0 - np.median(DEC_ran))*np.pi/180.),np.median(RA_ran)*np.pi/180.)
pixsjk = hp.get_all_neighbours(nside,pixjk)
pixsjk = np.append(pixsjk,pixjk)

ra_gal_all = ra_gal_all[np.in1d(pix,pixsjk)]
dec_gal_all = dec_gal_all[np.in1d(pix,pixsjk)]
mag_gal_all = mag_gal_all[np.in1d(pix,pixsjk)]
ra_gal_all[ra_gal_all>180] -= 360.

# this is just to count the area
JK_gal_ran = galaxy_randoms['JK']
N_allgal = len(JK_gal_ran)
mask = (JK_gal_ran==jkid)
num_gal_ran = galaxy_randoms['RA'][mask]
N_jkgal = len(num_gal_ran)
area_smalljk = tot_area*(N_jkgal*1.0/N_allgal)
print("area:", area_smalljk)

ra_gal_ran_all = galaxy_randoms['RA']
dec_gal_ran_all = galaxy_randoms['DEC']

#ids_ran = np.arange(len(ra_gal_ran_all))
#np.random.shuffle(ids_ran)
print('original gal', len(ra_gal_ran_all))
#ids_ran2 = ids_ran[:70000000]
#print('here')
ra_gal_ran_all2 = ra_gal_ran_all[:70000000]
dec_gal_ran_all2 = dec_gal_ran_all[:70000000]
dec_gal_ran_all = dec_gal_ran_all2.copy()
ra_gal_ran_all = ra_gal_ran_all2.copy()
N_allgal = len(dec_gal_ran_all)
print('new gal', N_allgal)

theta = (90.0 - dec_gal_ran_all)*np.pi/180.
phi = ra_gal_ran_all*np.pi/180.
pix_ran = hp.ang2pix(nside,theta,phi)
ra_gal_ran_all = ra_gal_ran_all[np.in1d(pix_ran, pixsjk)]
dec_gal_ran_all = dec_gal_ran_all[np.in1d(pix_ran, pixsjk)]
ra_gal_ran_all[ra_gal_ran_all>180] -= 360.

area_bigjk = tot_area*(len(ra_gal_ran_all)*1.0/N_allgal)
print("area:", area_bigjk)

# make sure we're not occupying memory
clusters=0
randoms=0
galaxies=0
galaxy_randoms=0

n1 = np.histogram(Z, range=(zmin,zmax), bins=Nz)
zmid = (n1[1][1:]+n1[1][:-1])/2

Area = []
Nclust = []
Nclust_ran = []
Nclust_ran_w = []
Ngal = []
Ngal_ran = []
Ngal_jk = []
DD = []
RR = []
DR = []
RD = []

for i in range(Nz):
    # treecorr    
    print("bin:", i, "mean z:", zmid[i])    

    mask = (Z>=n1[1][i])*(Z<n1[1][i+1])
    mask_ran = (Z_ran>=n1[1][i])*(Z_ran<n1[1][i+1])
    
    z_cen_bin = Z[mask]
    lamb_cen_bin = LAMB[mask]
    ra_cen_bin = RA[mask]
    dec_cen_bin = DEC[mask] 

    z_ran_bin = Z_ran[mask_ran]
    ra_ran_bin = RA_ran[mask_ran]
    dec_ran_bin = DEC_ran[mask_ran]  
    weight_ran_bin = W_ran[mask_ran] 

    M = mag - 5*(np.log10(cosmo.luminosity_distance(zmid[i]).value*1e6) - 1)
    M = M - 5.0*np.log10(0.7)
    mask_gal = (M>Maglim1)*(M<Maglim2)
    ngal_jk = len(mag[mask_gal])

    M = mag_gal_all - 5*(np.log10(cosmo.luminosity_distance(zmid[i]).value*1e6) - 1)
    M = M - 5.0*np.log10(0.7)
    mask_gal_all = (M>Maglim1)*(M<Maglim2)
    ra_gal_bin = ra_gal_all[mask_gal_all]
    dec_gal_bin = dec_gal_all[mask_gal_all]

    # Convert physical to angular distance at this zl
    D_l = cosmo.comoving_distance(zmid[i]).value # comiving distance
    #D_l = cosmo.comoving_distance(zmid[i]) / (1.+ zmid[i]) # physical distance
    thmin = np.arctan(Rmin / D_l) * (180./np.pi) * 60.      #arcmin
    thmax = np.arctan(Rmax / D_l) * (180./np.pi) * 60.	

    n_clust = len(ra_cen_bin)
    n_clust_ran = len(ra_ran_bin)
    n_clust_ran_w = np.sum(weight_ran_bin)
    n_gal = len(ra_gal_bin)
    n_gal_ran = len(ra_gal_ran_all)

    #area_Mpch = area_bigjk*(np.pi/180.)**2*(cosmo.comoving_distance(zmid[i]).value)**2
    area_Mpch = area_smalljk*(np.pi/180.)**2*(cosmo.comoving_distance(zmid[i]).value)**2

    Area.append(area_Mpch)
    Nclust.append(n_clust)
    Nclust_ran.append(n_clust_ran)
    Nclust_ran_w.append(n_clust_ran_w)
    Ngal.append(n_gal)
    Ngal_ran.append(n_gal_ran)
    Ngal_jk.append(ngal_jk)

    # Define catalogs
    print(len(ra_gal_bin), len(ra_gal_ran_all))

    ran_cat = treecorr.Catalog(ra=ra_ran_bin, dec=dec_ran_bin,
                                    ra_units='degrees', dec_units='degrees', w=weight_ran_bin)
    gal_cat = treecorr.Catalog(ra=ra_gal_bin, dec=dec_gal_bin,
                                   ra_units='degrees', dec_units='degrees')
    gal_ran_cat = treecorr.Catalog(ra=ra_gal_ran_all, dec=dec_gal_ran_all,
                                   ra_units='degrees', dec_units='degrees')

    rd = treecorr.NNCorrelation(nbins = nR, min_sep = thmin, max_sep = thmax,
                                        bin_slop = bslop, sep_units = 'arcmin', verbose=2) 

    rr = treecorr.NNCorrelation(nbins = nR, min_sep = thmin, max_sep = thmax,
                                        bin_slop = bslop, sep_units = 'arcmin', verbose=2) 

    rd.process(ran_cat, gal_cat)
    rr.process(ran_cat, gal_ran_cat)

    RR.append(rr.weight)
    RD.append(rd.weight)

    if len(z_cen_bin)>0:
    
        cen_cat = treecorr.Catalog(ra=ra_cen_bin, dec=dec_cen_bin,
                                    ra_units='degrees', dec_units='degrees')
        dd = treecorr.NNCorrelation(nbins = nR, min_sep = thmin, max_sep = thmax,
                                        bin_slop = bslop, sep_units = 'arcmin', verbose=2)
        dr = treecorr.NNCorrelation(nbins = nR, min_sep = thmin, max_sep = thmax,
                                        bin_slop = bslop, sep_units = 'arcmin', verbose=2)
        dd.process(cen_cat, gal_cat)    # For cross-correlation.
        dr.process(cen_cat, gal_ran_cat)
        DD.append(dd.weight)
        DR.append(dr.weight)


    else:
        DD.append(np.zeros(nR))
        DR.append(np.zeros(nR))


Area = np.array(Area)
RR = np.array(RR)
DD = np.array(DD)
DR = np.array(DR)
RD = np.array(RD)
Nclust = np.array(Nclust)
Nclust_ran = np.array(Nclust_ran)
Nclust_ran_w = np.array(Nclust_ran_w)
Ngal = np.array(Ngal)
Ngal_ran = np.array(Ngal_ran)
Ngal_jk = np.array(Ngal_jk)

np.savez(outfile, R=R_mid, area=Area, 
nclust=Nclust, nclust_ran=Nclust_ran, nclust_ran_w=Nclust_ran_w, ngal=Ngal, ngal_ran=Ngal_ran, ngal_jk=Ngal_jk,
dd=DD, rr=RR, dr=DR, rd=RD)


