
import healpy as hp
import numpy as np
import scipy 
import astropy.io.fits as pf
import kmeans_radec


def radec2thetaphi(ra, dec, nside):
    """
    Convert RA DEC in degrees to THETA and PHI in Healpix 
    convention. 
    """

    theta = (90.0 - dec)*np.pi/180.0
    phi = ra*np.pi/180.0
    return theta, phi

def make_mask(ra, dec, nmin=1, nside=4096):
    """
    Take RA, DEC, build a binary mask just by assigning 
    1 to pixels with count>=nmin and 0 otherwise. Mask 
    is in Healpix format with assigned nside. 
    """
    mask = np.zeros(hp.nside2npix(nside))
    theta, phi = radec2thetaphi(ra, dec, nside)
    pix = hp.ang2pix(nside, theta, phi, nest=False)
    for i in range(len(pix)):
        mask[pix[i]] += 1
    mask[mask>=nmin] = 1
    mask[mask!=1] = 0

    return mask

def make_random(mask, ramin, ramax, decmin, decmax, N=1000, nside=4096, seed=100):
    """
    Create N random points within a given mask.
    """

    #print 'Building the random catalog...'

    np.random.seed(seed)
    ra_rand = (np.random.random(N)* (ramax - ramin))+ramin
    v = np.random.random(N)
    #print len(v)
    vmin = np.cos((90.0+decmin)/180.*np.pi)
    vmax = np.cos((90.0+decmax)/180.*np.pi)
    v *= (vmax-vmin)
    v += vmin
    #v *= 2
    #v -= 1
    dec_rand = np.arccos(v)
    np.rad2deg(dec_rand,dec_rand)
    dec_rand -= 90.0
    #dec_rand_se =((dec_rand < decmax)&(dec_rand > decmin))
    #ra_rand, dec_rand = ra_rand[dec_rand_se],dec_rand[dec_rand_se]
    #print len(ra_rand)
    #print 'Masking the random catalog...'

    #Converting degrees into radians
    theta_rand = (90.0 - dec_rand)*np.pi/180.
    phi_rand = ra_rand*np.pi/180.
    pix_rand = hp.ang2pix(nside, theta_rand, phi_rand, nest=False)

    goodm, = np.where(mask[pix_rand]==1)
    ra_rand = ra_rand[goodm]
    dec_rand = dec_rand[goodm]
    #print len(ra_rand)

    return ra_rand, dec_rand


def make_jk(ra_ran, dec_ran, ra, dec, N=100, dilute_factor=1, rand_out=1, large_mem=True, maxiter=500, tol=1e-05, seed=100):
    """
    Given coordinate of random points, generate JK indecies 
    for another catalog of positions. Include the possibility 
    of diluting the random catalog. Return an array of JK 
    indicies the same length of ra and dec.  
    """

    RADEC_ran = np.zeros((len(ra_ran),2))
    RADEC_ran[:,0] = ra_ran
    RADEC_ran[:,1] = dec_ran

    RADEC = np.zeros((len(ra),2))
    RADEC[:,0] = ra
    RADEC[:,1] = dec

    np.random.seed(seed)
    ids = np.arange(len(ra_ran))
    np.random.shuffle(ids)
    RADEC_ran_dilute = np.zeros((len(ra_ran)/dilute_factor,2))
    RADEC_ran_dilute[:,0] = ra_ran[ids[:len(ra_ran)/dilute_factor]]
    RADEC_ran_dilute[:,1] = dec_ran[ids[:len(ra_ran)/dilute_factor]]

    km = kmeans_radec.kmeans_sample(RADEC_ran_dilute, N, maxiter=500, tol=1e-05)
    print(np.unique(km.labels))

    if large_mem == True:
        Ntotal = len(RADEC)
        Ntotal_ran = len(RADEC_ran)

        JK = np.array([])
        JK_ran = np.array([])

        for i in range(99):
            #print i
            JK = np.concatenate((JK, km.find_nearest(RADEC[i*(Ntotal/100):(i+1)*(Ntotal/100)])), axis=0)
            print(np.unique(JK))

            if rand_out==1:
                print(i)
                JK_ran = np.concatenate((JK_ran, km.find_nearest(RADEC_ran[i*(Ntotal_ran/100):(i+1)*(Ntotal_ran/100)])), axis=0)

        JK = np.concatenate((JK, km.find_nearest(RADEC[99*(Ntotal/100):])), axis=0)
        if rand_out==1:
            JK_ran = np.concatenate((JK_ran, km.find_nearest(RADEC_ran[99*(Ntotal_ran/100):])), axis=0)
        print('len of random', len(ra_ran))
        print('len of JK', len(JK_ran))

    else:
        JK = km.find_nearest(RADEC)
        if rand_out==1:
            JK_ran = km.find_nearest(RADEC_ran)
    
    if rand_out==1:    
        return JK_ran, JK
    else:
        return JK, JK




def make_jk_from_random(ra_ran, dec_ran, N=100, dilute_factor=1, maxiter=500, tol=1e-05, seed=100):
    """
    Given coordinate of random points, generate JK indecies.
    """

    RADEC_ran = np.zeros((len(ra_ran),2))
    RADEC_ran[:,0] = ra_ran
    RADEC_ran[:,1] = dec_ran

    np.random.seed(seed)
    ids = np.arange(len(ra_ran))
    np.random.shuffle(ids)
    RADEC_ran_dilute = np.zeros((len(ra_ran)/dilute_factor,2))
    RADEC_ran_dilute[:,0] = ra_ran[ids[:len(ra_ran)/dilute_factor]]
    RADEC_ran_dilute[:,1] = dec_ran[ids[:len(ra_ran)/dilute_factor]]

    km = kmeans_radec.kmeans_sample(RADEC_ran_dilute, N, maxiter=500, tol=1e-05)
    
    return km, RADEC_ran_dilute[:,0], RADEC_ran_dilute[:,1]

