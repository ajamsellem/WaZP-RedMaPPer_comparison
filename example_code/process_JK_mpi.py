
import numpy as np
import os 
import sys

# this is to combine all the JK files and output the mean and covariance

name = sys.argv[1]
Njk = int(sys.argv[2])
Nz = int(sys.argv[3])
result_dir = sys.argv[4]

R = np.load(result_dir+name+'/Sigmag_0.npz')['R']

sum_RR = np.zeros((Nz, len(R)))
sum_DD = np.zeros((Nz, len(R)))
sum_DR = np.zeros((Nz, len(R)))
sum_RD = np.zeros((Nz, len(R)))
sum_nclust = np.zeros(Nz)
sum_nclust_ran = np.zeros(Nz)
sum_ngal = np.zeros(Nz)
sum_ngal_ran = np.zeros(Nz)
sum_area = np.zeros(Nz)
sum_ngal_jk = np.zeros(Nz)

for i in range(Njk):
    infile = np.load(result_dir+name+'/Sigmag_'+str(i)+'.npz')
    sum_RR += infile['rr']
    sum_DD += infile['dd']
    sum_DR += infile['dr']
    sum_RD += infile['rd']
    sum_area += infile['area']
    sum_nclust += infile['nclust']
    sum_nclust_ran += infile['nclust_ran_w']
    sum_ngal += infile['ngal']
    sum_ngal_ran += infile['ngal_ran']
    sum_ngal_jk += infile['ngal_jk']


Sg = []

for i in range(Njk):
    infile = np.load(result_dir+name+'/Sigmag_'+str(i)+'.npz')

    DD = sum_DD - infile['dd']
    RR = sum_RR - infile['rr']
    DR = sum_DR - infile['dr']
    RD = sum_RD - infile['rd']

    RR = RR*1.0*((sum_nclust - infile['nclust'])*1.0/(sum_nclust_ran - infile['nclust_ran_w'])*(sum_ngal - infile['ngal'])*1.0/(sum_ngal_ran - infile['ngal_ran'])).reshape(Nz,1)
    DR = DR*1.0*((sum_ngal - infile['ngal'])*1.0/(sum_ngal_ran - infile['ngal_ran'])).reshape(Nz,1)
    RD = RD*1.0*((sum_nclust - infile['nclust'])*1.0/(sum_nclust_ran - infile['nclust_ran_w'])).reshape(Nz,1)
    ave_dens = (sum_ngal_jk-infile['ngal_jk'])*1.0/(sum_area-infile['area'])
    nclust = (sum_nclust - infile['nclust'])*1.0

    #print ave_dens
    dens = np.sum(ave_dens*nclust)/np.sum(nclust)
    mean = np.sum(DD-DR-RD+RR, axis=0)/np.sum(RR, axis=0)*dens  

    Sg.append(mean)

Sg = np.array(Sg)
C = np.cov(Sg.T)*(len(Sg)-1)
Sg_mean = np.mean(Sg, axis=0)
Sg_sig = np.sum((Sg-Sg_mean)**2, axis=0)**0.5/len(Sg)**0.5*(len(Sg)-1)**0.5

np.savez('splashback_cov_'+str(name)+'.npz', cov=C, r_data=R, sg_mean=Sg_mean, sg_sig=Sg_sig) 


