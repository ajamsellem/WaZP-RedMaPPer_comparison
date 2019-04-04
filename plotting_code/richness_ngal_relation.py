# Plot Richness (RM) versus Number of Galaxies (WaZP) and find relationship (using errors from RM)
%matplotlib inline
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

# Arrays with richness and richness error from RM and WaZP
rich_u = np.asarray(rich_u)
rich = np.asarray(rich) 
drich = np.asarray(drich)
drich_u = np.asarray(drich_u)

# Fitting and Residual Functions
def linear(p, xvar):
    return p[0] + p[1]*xvar
def linear_residual(p, xvar, yvar, err):
    return (linear(p, xvar) - yvar)/err
def quadratic(p, xvar):
    return p[0] + p[1]*xvar + p[2]*(xvar**2)
def quadratic_residual(p, xvar, yvar, err):
    return (quadratic(p, xvar) - yvar)/err

# Initial guesses for fit parameters
p01 = [1., 1.]
p02 = [1., 1., 1.]

# Perform fit
pf1, cov1, info1, mesg1, success1 = optimize.leastsq(linear_residual, p01,
                                                args=(rich_u, rich, drich), full_output=1)
pf2, cov2, info2, mesg2, success2 = optimize.leastsq(quadratic_residual, p02,
                                                args=(rich_u, rich, drich), full_output=1)

# Linear Fit
# If the fit failed, print the reason
if cov1 is None:
    print('Fit did not converge')
    print('Success code:', success1)
    print(mesg1)
else:
    chisq1 = sum(info1['fvec']*info1['fvec'])
    dof1 = len(rich_u) - len(pf1)
    pferr1 = [np.sqrt(cov1[i,i]) for i in range(len(pf1))]
    print('These are the results of the Linear Fit:')
    print('Converged with chi-squared', chisq1)
    print('Number of degrees of freedom, dof =', dof1)
    print('Reduced chi-squared ', chisq1/dof1)
    print('Inital guess values:')
    print('  p0 =', p01)
    print('Best fit values:')
    print('  pf =', pf1)
    print('Uncertainties in the best fit values:')
    print('  pferr =', pferr1)
      
# Quadratic Fit
# If the fit failed, print the reason
if cov2 is None:
    print('Fit did not converge')
    print('Success code:', success2)
    print(mesg2)
else:
    chisq2 = sum(info2['fvec']*info2['fvec'])
    dof2 = len(rich_u) - len(pf2)
    pferr2 = [np.sqrt(cov2[i,i]) for i in range(len(pf2))]
    print('\nThese are the results of the Quadratic Fit:')
    print('Converged with chi-squared', chisq2)
    print('Number of degrees of freedom, dof =', dof2)
    print('Reduced chi-squared ', chisq2/dof2)
    print('Inital guess values:')
    print('  p0 =', p02)
    print('Best fit values:')
    print('  pf =', pf2)
    print('Uncertainties in the best fit values:')
    print('  pferr =', pferr2)

# Plot Richness (RM) v. Number of Galaxies (WaZP)
X = np.linspace(rich_u.min(), rich_u.max(), 500)

fig1 = plt.figure(figsize = (11,8))
ax1 = fig1.add_subplot(111)
ax1.errorbar(rich_u, rich, yerr = drich, fmt = 'r.', capsize = 2, zorder = 0, label = "Cluster")
ax1.plot(X, linear(pf1, X), 'g-', label = 'Linear Fit', zorder = 1)
ax1.plot(X, quadratic(pf2, X), 'b-', label = 'Quadratic Fit', zorder = 3)

plt.xlabel("Number of Galaxies (WaZP)", fontsize = 15)
plt.ylabel("Richness (Redmapper)", fontsize = 15)

texquad = 'Quadratic Fit: \n $r_{RM} = %.2f \\cdot N_{W}^2$ + $%.2f \\cdot N_{W} + %.2f$ \n $Reduced \: \chi^2 = % .2f$ \n' % (pf2[2], pf2[1], pf2[0], chisq2/dof2)
texline = 'Linear Fit: \n $r_{RM} = %.2f \\cdot N_{W} + %.2f$ \n $Reduced \: \\chi^2 = % .2f$' % (pf1[1], pf1[0], chisq1/dof1)
facts = '$r_{RM} \\equiv$ Redmapper Richness \n' \
        '$N_{W} \\equiv$ WaZP Number of Galaxies \n' \

ax1.text(0.04, .95, texline, transform=ax1.transAxes, fontsize=17, verticalalignment='top')
ax1.text(0.38, .95, texquad, transform=ax1.transAxes, fontsize=17, verticalalignment='top')
ax1.text(0.34, .752, facts, transform=ax1.transAxes, fontsize=17, verticalalignment='top')

txt= "Errors used were from Redmapper data."
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='right', fontsize=17)

plt.legend(loc = 'center left', fontsize =17, bbox_to_anchor=(0., 0.65))

fig1.savefig('./Richness_Ngal_compare_(RM_Err).png', bbox_inches="tight")
plt.show()


# Plot Richness (RM) versus Number of Galaxies (WaZP) and find relationship (using errors from WaZP)
%matplotlib inline
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

rich_u = np.asarray(rich_u)
rich = np.asarray(rich) 
drich = np.asarray(drich)
drich_u = np.asarray(drich_u)

def linear(p, xvar):
    return p[0] + p[1]*xvar
def linear_residual(p, xvar, yvar, err):
    return (linear(p, xvar) - yvar)/err
def quadratic(p, xvar):
    return p[0] + p[1]*xvar + p[2]*(xvar**2)
def quadratic_residual(p, xvar, yvar, err):
    return (quadratic(p, xvar) - yvar)/err

p01 = [1., 1.]
p02 = [1., 1., 1.]

pf1, cov1, info1, mesg1, success1 = optimize.leastsq(linear_residual, p01,
                                                args=(rich, rich_u, drich_u), full_output=1)
pf2, cov2, info2, mesg2, success2 = optimize.leastsq(quadratic_residual, p02,
                                                args=(rich, rich_u, drich_u), full_output=1)

# Linear Fit
# If the fit failed, print the reason
if cov1 is None:
    print('Fit did not converge')
    print('Success code:', success1)
    print(mesg1)
else:
    chisq1 = sum(info1['fvec']*info1['fvec'])
    dof1 = len(rich) - len(pf1)
    pferr1 = [np.sqrt(cov1[i,i]) for i in range(len(pf1))]
    print('These are the results of the Linear Fit:')
    print('Converged with chi-squared', chisq1)
    print('Number of degrees of freedom, dof =', dof1)
    print('Reduced chi-squared ', chisq1/dof1)
    print('Inital guess values:')
    print('  p0 =', p01)
    print('Best fit values:')
    print('  pf =', pf1)
    print('Uncertainties in the best fit values:')
    print('  pferr =', pferr1)
     
# Quadratic Fit
# If the fit failed, print the reason
if cov2 is None:
    print('Fit did not converge')
    print('Success code:', success2)
    print(mesg2)
else:
    chisq2 = sum(info2['fvec']*info2['fvec'])
    dof2 = len(rich_u) - len(pf2)
    pferr2 = [np.sqrt(cov2[i,i]) for i in range(len(pf2))]
    print('\nThese are the results of the Quadratic Fit:')
    print('Converged with chi-squared', chisq2)
    print('Number of degrees of freedom, dof =', dof2)
    print('Reduced chi-squared ', chisq2/dof2)
    print('Inital guess values:')
    print('  p0 =', p02)
    print('Best fit values:')
    print('  pf =', pf2)
    print('Uncertainties in the best fit values:')
    print('  pferr =', pferr2)

# Plot Richness v. Ngal
X = np.linspace(rich.min(), rich.max(), 500)

fig1 = plt.figure(figsize = (11,8))
ax1 = fig1.add_subplot(111)
ax1.errorbar(rich, rich_u, yerr = drich_u, fmt = 'r.', capsize = 2, zorder = 0, label = "Cluster")
ax1.plot(X, linear(pf1, X), 'g-', label = 'Linear Fit', zorder = 1)
ax1.plot(X, quadratic(pf2, X), 'b-', label = 'Quadratic Fit', zorder = 3)
plt.ylabel("Number of Galaxies (WaZP)", fontsize = 15)
plt.xlabel("Richness (Redmapper)", fontsize = 15)

texquad = 'Quadratic Fit: \n $r_{RM} = %.2f \\cdot N_{W}^2$ + $%.2f \\cdot N_{W} + %.2f$ \n $Reduced \: \chi^2 = % .2f$ \n' % (pf2[2], pf2[1], pf2[0], chisq2/dof2)
texline = 'Linear Fit: \n $r_{RM} = %.2f \\cdot N_{W} + %.2f$ \n $Reduced \: \\chi^2 = % .2f$' % (pf1[1], pf1[0], chisq1/dof1)
facts = '$r_{RM} \\equiv$ Redmapper Richness \n' \
        '$N_{W} \\equiv$ WaZP Number of Galaxies \n' \

ax1.text(0.04, .95, texline, transform=ax1.transAxes, fontsize=17, verticalalignment='top')
ax1.text(0.38, .95, texquad, transform=ax1.transAxes, fontsize=17, verticalalignment='top')
ax1.text(0.1, .78, facts, transform=ax1.transAxes, fontsize=17, verticalalignment='top')

txt= "Errors used were from WaZP data."
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='right', fontsize=17)

plt.legend(loc = 'lower right', fontsize =17)

fig1.savefig('./Richness_Ngal_compare_(WaZP_Err).png', bbox_inches="tight")
plt.show()


# Plot Richness (RM) (in range 20. to 100.) versus Number of Galaxies (WaZP) and find relationship (using errors from RM)
rich_new = []
drich_new = []
rich_u_new = []

# Obtain richnesses from RM data between 20. and 100.
for i, value in enumerate(rich):
    if value >= 20. and value <= 100.:
        rich_new.append(value)
        rich_u_new.append(rich_u[i])
        drich_new.append(drich[i])
        
rich_u_new = np.asarray(rich_u_new)
rich_new = np.asarray(rich_new) 
drich_new = np.asarray(drich_new)
        
pf1, cov1, info1, mesg1, success1 = optimize.leastsq(linear_residual, p01,
                                                args=(rich_u_new, rich_new, drich_new), full_output=1)
pf2, cov2, info2, mesg2, success2 = optimize.leastsq(quadratic_residual, p02,
                                                args=(rich_u_new, rich_new, drich_new), full_output=1)

# Linear Fit
# If the fit failed, print the reason
if cov1 is None:
    print('Fit did not converge')
    print('Success code:', success1)
    print(mesg1)
else:
    chisq1 = sum(info1['fvec']*info1['fvec'])
    dof1 = len(rich_u_new) - len(pf1)
    pferr1 = [np.sqrt(cov1[i,i]) for i in range(len(pf1))]
    print('These are the results of the Linear Fit:')
    print('Converged with chi-squared', chisq1)
    print('Number of degrees of freedom, dof =', dof1)
    print('Reduced chi-squared ', chisq1/dof1)
    print('Inital guess values:')
    print('  p0 =', p01)
    print('Best fit values:')
    print('  pf =', pf1)
    print('Uncertainties in the best fit values:')
    print('  pferr =', pferr1)
    
# Quadratic Fit
# If the fit failed, print the reason
if cov2 is None:
    print('Fit did not converge')
    print('Success code:', success2)
    print(mesg2)
else:
    chisq2 = sum(info2['fvec']*info2['fvec'])
    dof2 = len(rich_u_new) - len(pf2)
    pferr2 = [np.sqrt(cov2[i,i]) for i in range(len(pf2))]
    print('\nThese are the results of the Quadratic Fit:')
    print('Converged with chi-squared', chisq2)
    print('Number of degrees of freedom, dof =', dof2)
    print('Reduced chi-squared ', chisq2/dof2)
    print('Inital guess values:')
    print('  p0 =', p02)
    print('Best fit values:')
    print('  pf =', pf2)
    print('Uncertainties in the best fit values:')
    print('  pferr =', pferr2)

# Plot Richness v. Ngal
X = np.linspace(rich_u_new.min(), rich_u_new.max(), 500)

fig1 = plt.figure(figsize = (11,8))
ax1 = fig1.add_subplot(111)
ax1.errorbar(rich_u_new, rich_new, yerr = drich_new, fmt = 'r.', capsize = 2, zorder = 0, label = "Cluster")
ax1.plot(X, linear(pf1, X), 'g-', label = 'Linear Fit', zorder = 1)
ax1.plot(X, quadratic(pf2, X), 'b-', label = 'Quadratic Fit', zorder = 3)

plt.xlabel("Number of Galaxies (WaZP)", fontsize = 15)
plt.ylabel("Richness (Redmapper)", fontsize = 15)

texquad = 'Quadratic Fit: \n $r_{RM} = %.2f \\cdot N_{W}^2$ + $%.2f \\cdot N_{W} + %.2f$ \n $Reduced \: \chi^2 = % .2f$ \n' % (pf2[2], pf2[1], pf2[0], chisq2/dof2)
texline = 'Linear Fit: \n $r_{RM} = %.2f \\cdot N_{W} + %.2f$ \n $Reduced \: \\chi^2 = % .2f$' % (pf1[1], pf1[0], chisq1/dof1)
facts = '$r_{RM} \\equiv$ Redmapper Richness \n' \
        '$N_{W} \\equiv$ WaZP Number of Galaxies \n' \

txt = "Errors used were from Redmapper data."
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='right', fontsize=17)

plt.legend(loc = 'upper left', fontsize =17)

fig1.savefig('./Richness_Ngal_compare_(20-100).png', bbox_inches="tight")
plt.show()
