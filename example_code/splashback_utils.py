import numpy as np
import pdb
from scipy import interpolate
import scipy.integrate as integrate
import scipy.signal as sig

def Sigmag(R, z, params, h0, splash):
    
    minr = 0.01
    maxr = 100.
    numr = 500  
    rr = np.exp(np.linspace(np.log(minr), np.log(maxr), num = numr))   

    if splash==1:
        ln_alpha, ln_beta, ln_gamma, ln_r_s, ln_r_t, ln_rho_O, ln_rho_s, se = params
        beta = 10**ln_beta
        gamma = 10**ln_gamma
        r_t = 10**ln_r_t
        f_trans = (1.+(rr/r_t)**beta)**(-1*gamma/beta)

    if splash==0:
        ln_alpha, ln_r_s, ln_rho_O, ln_rho_s, se = params
        f_trans = 1.0

    alpha = 10**ln_alpha
    r_s = 10**ln_r_s
    rho_O = 10**ln_rho_O
    rho_s = 10**ln_rho_s    
    r_o = 1.5/h0

    rho_gi = rho_s*np.exp(-2./alpha*((rr/r_s)**alpha-1))
    rho_go = rho_O*(rr/r_o)**(-1*se)
    rho_g = rho_gi * f_trans + rho_go
    rho_g_func = interpolate.interp1d(rr, rho_g)

    sigmag = []
    for i in range(len(R)):
        func_evals = rho_g_func(np.sqrt(R[i]**2.+z**2.))
        sigmag.append(2*integrate.simps(func_evals, z))
        # it appears to make a difference how this integration is done... 

    return np.array(sigmag)


def priors(params, h0, splash):

    if splash==1:
        ln_alpha, ln_beta, ln_gamma, ln_r_s, ln_r_t, ln_rho_O, ln_rho_s, se = params
        beta = 10**ln_beta
        gamma = 10**ln_gamma
        r_t = 10**ln_r_t

    if splash==0:
        ln_alpha, ln_r_s, ln_rho_O, ln_rho_s, se = params

    alpha = 10**ln_alpha
    r_s = 10**ln_r_s
    rho_O = 10**ln_rho_O
    rho_s = 10**ln_rho_s  

    min_logrs = np.log10(0.1/h0)
    max_logrs = np.log10(5.0/h0)
    min_logrt = np.log10(0.1/h0)
    max_logrt = np.log10(5.0/h0)
    min_se = -10.0
    max_se = 10.0
    
    if splash==1:
        if (ln_r_s > max_logrs) or (ln_r_s < min_logrs) or (ln_r_t > max_logrt) or (ln_r_t < min_logrt)  or (se < min_se) or (se > max_se):
            lnprior = -1.0e10
        
        else:
            mean_logalpha = np.log10(0.2)
            sigma_logalpha = 0.6
            mean_logbeta = np.log10(4.)
            sigma_logbeta = 0.2
            mean_loggamma = np.log10(6.0)
            sigma_loggamma = 0.2
            
            lnprior_alpha = -0.5*np.log(2.*np.pi*sigma_logalpha**2.)-0.5*((ln_alpha - mean_logalpha)**2.)/sigma_logalpha**2.
            lnprior_beta =  -0.5*np.log(2.*np.pi*sigma_logbeta**2.)-0.5*((ln_beta - mean_logbeta)**2.)/sigma_logbeta**2.
            lnprior_gamma =  -0.5*np.log(2.*np.pi*sigma_loggamma**2.)-0.5*((ln_gamma- mean_loggamma)**2.)/sigma_loggamma**2.
            lnprior = lnprior_alpha + lnprior_beta + lnprior_gamma

    if splash==0:

        if ((np.log10(r_s) > max_logrs) or (np.log10(r_s) < min_logrs) or (se < min_se) or (se > max_se)):
            lnprior = -1.0e10
        
        else:
            mean_logalpha = np.log10(0.2)
            sigma_logalpha = 0.6
            lnprior_alpha = -0.5*np.log(2.*np.pi*sigma_logalpha**2.)-0.5*((np.log10(alpha) - mean_logalpha)**2.)/sigma_logalpha**2.
            lnprior = lnprior_alpha 

    if (np.isnan(lnprior)):
        pdb.set_trace()

    return lnprior


def lnlikelihood(params, R, z, data_vec, invcov, h0, splash):

    lnlike_priors = priors(params, h0, splash)
    lnlike_data = 0.0
    
    if (lnlike_priors > -1.0e5):
        
        model = Sigmag(R, z, params, h0, splash)
        diff = data_vec - model
        detinvcov = np.linalg.det(invcov)
        detcov = 1./detinvcov
        lnlike_data = -0.5*(len(data_vec)*np.log(2.*np.pi) + np.log(detcov)) -0.5*np.dot(diff, np.dot(invcov, diff))
    
    lnlike = lnlike_data + lnlike_priors
       
    return lnlike


def derivative_savgol(R, data, N=10000, window_length=5, polyorder=3):
    
    data_sm = sig.savgol_filter(np.log10(data), window_length=window_length, polyorder=polyorder)
    f = interpolate.interp1d(np.log10(R), data_sm, kind='cubic')
    
    # Evaluate spline across a very fine grid or radii
    lnrad_fine = np.linspace(np.log10(np.min(R)), np.log10(np.max(R)), num=N)
    lnsigma_fine = f(lnrad_fine)
    
    # Calculate derivative using finite differencing
    dlnsig_dlnr_fine = (lnsigma_fine[1:] - lnsigma_fine[:-1])/(lnrad_fine[1:] - lnrad_fine[:-1])
       
    return (lnrad_fine[1:]+lnrad_fine[:-1])/2, dlnsig_dlnr_fine

# Alternative to DSigmag, results are basically identical...
def DelSigmag(R, z, params, h0, splash, N=100):
    
    sigma = Sigmag(R, z, params, h0, splash)
    # interpolate onto a finer grid just so the integration is smooth
    
    f = interpolate.interp1d(np.log10(R), np.log10(sigma), kind='cubic')
    rad_fine = np.linspace(np.log10(np.min(R)), np.log10(np.max(R)), num=N)    

    Dsigma = []
    for i in range(len(rad_fine)-1):
        func_evals = f(rad_fine[:i+1])*2.*np.pi*rad_fine[:i+1]
        sigmag_sum = integrate.simps(func_evals, rad_fine[:i+1])
        sigmag_mean = sigmag_sum/(np.pi*rad_fine[i+1]**2)
        Dsigma.append(sigmag_mean - f(rad_fine[i+1]))

    f = interpolate.interp1d(rad_fine[1:], Dsigma, kind='cubic', bounds_error=False, fill_value=0)
    lnsigma_coarse = f(R) 
    
    return lnsigma_coarse

def DSigmag(R, z, params, h0, splash, N=100):
    
    sigma = Sigmag(R, z, params, h0, splash)
    # interpolate onto a finer grid just so the integration is smooth
    f = interpolate.interp1d(np.log10(R), np.log10(sigma), kind='cubic')
    lnrad_fine = np.linspace(np.min(np.log10(R)), np.max(np.log10(R)), num=N)
    lnsigma_fine = f(lnrad_fine)   
    
    R_fine = 10**lnrad_fine
    sigma_fine = 10**lnsigma_fine   
    R_fine_mid = (R_fine[1:]+R_fine[:-1])/2
    dR_fine = R_fine[1:]-R_fine[:-1]
    sigma_fine_mid = (sigma_fine[1:]+sigma_fine[:-1])/2
    
    Dsigma = []
    
    for i in range(len(R_fine_mid)):
        Mean = np.sum(sigma_fine_mid[:i+1]*2*np.pi*R_fine_mid[:i+1]*dR_fine[:i+1])/np.sum(2*np.pi*R_fine_mid[:i+1]*dR_fine[:i+1])
        Dsigma.append(Mean-sigma_fine_mid[i])
             
    # interpolate back to original R grid
    f = interpolate.interp1d(R_fine_mid, Dsigma, kind='cubic', bounds_error=False, fill_value=0)
    lnsigma_coarse = f(R) 
    
    return lnsigma_coarse


def lnlikelihoodD(params, R, z, data_vec, invcov, h0, splash):

    lnlike_priors = priors(params, h0, splash)
    lnlike_data = 0.0
    
    if (lnlike_priors > -1.0e5):
        
        model = DSigmag(R, z, params, h0, splash)
        diff = data_vec - model
        detinvcov = np.linalg.det(invcov)
        detcov = 1./detinvcov
        lnlike_data = -0.5*(len(data_vec)*np.log(2.*np.pi) + np.log(detcov)) -0.5*np.dot(diff, np.dot(invcov, diff))
    
    lnlike = lnlike_data + lnlike_priors
       
    return lnlike
