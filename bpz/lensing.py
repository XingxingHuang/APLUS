#Este fichero modela lensing.
from bpz_tools import *
from cosmology import *
from biggles import *
from LinearAlgebra import *
#This module below also need the cephes module by Travis Oliphant
from quadpack import *

cho=2.99e3       #c/H_0 in Mpc
g=6.672e-8       #cm**3/g/s**2
c=2.99e10        #cm/s
mpc=3e24         #cm
solar_mass=2e33  #g

def d_filled_beam(z_l,z_s,cosmo=(0.3,0.7,.7)):
    #Reference Fukugita et al. 1992,393,3
    #Filled beam approximation
    k=cosmo[0]+cosmo[1]
    o=cosmo[0]
    l=cosmo[1]
    def ez(x,k=cosmo[0]+cosmo[1],o=cosmo[0],l=cosmo[1]): 
	return 1./sqrt(o*(1.+x)**3+(1.-k)*(1.+x)**2+l)    
    ch21=quad(ez,z_l,z_s)[0]
    
    if k>1:
	ch21=sqrt(abs(k-1.))*ch21
	return cho/cosmo[2]/sqrt(k-1.)/(1.+z_s)*sin(ch21)
    elif k==1:
	return cho/cosmo[2]/(1.+z_s)*ch21
    else:
	ch21=sqrt(abs(k-1.))*ch21
	return cho/cosmo[2]/sqrt(1.-k)/(1.+z_s)*sinh(ch21)	    

def d_dyer_roeder(z_l,z_s,cosmo=(0.3,0.7,.7)):
    #Reference Fukugita et al. 1992,393,3
    #Dyer-Roeder aproximation (empty beam)
    k=cosmo[0]+cosmo[1]
    o=cosmo[0]
    l=cosmo[1]
    def ez(x,k=cosmo[0]+cosmo[1],o=cosmo[0],l=cosmo[1]): 
	return 1./sqrt((o*(1.+x)**3+(1.-k)*(1.+x)**2+l))/(1.+x)**2    
    return cho/cosmo[2]*(1.+z_l)*quad(ez,z_l,z_s)[0]

#d_ang=d_dyer_roeder
d_ang=d_filled_beam

def einstein_radius(sigma,z_l,z_s,cosmo=(0.3,0.7,.7)):
    #Assuming Omega+Lambda=1
    theta=4.*pi*(sigma/3e5)**2*d_ang(z_l,z_s,cosmo)/d_ang(0.,z_s,cosmo)
    return theta/pi*180.*3600.

def sigma_einstein(e_radius,z_l,z_s,cosmo=(0.3,0.7,.7)):
    e_radius=e_radius/3600./180.*pi
    sigma=3e5*sqrt(e_radius/4./pi/d_ang(z_l,z_s,cosmo)*d_ang(0.,z_s,cosmo))
    return sigma

def einstein_mass(e_radius,z_l,z_s,cosmo=(0.3,0.7,.7)):
    sigma=sigma_einstein(e_radius,z_l,z_s,cosmo)
    e_radius=e_radius/3600./180.*pi
    mass=pi*(sigma*1e5)**2/g*e_radius*d_ang(0.,z_l,cosmo)*mpc/solar_mass
    return mass

def iso_lensing(x,y,m,z,theta_e=51.,z_e=1.5,z_l=0.183):
    #Transforms x,y and m for an isothermal sphere taking 
    #into accout the mass parameters given above
    masa=einstein_mass(theta_e,z_l,z_e)
    r=sqrt(x**2+y**2)
    cos_theta=x/r
    sin_theta=y/r
    thetae=einstein_radius(masa,z_l,z)
    mu=1.+thetae/r
    m=m-2.5*log10(mu)
    r=r+thetae
    x=r*cos_theta
    y=r*sin_theta
    return x,y,m

def kappafromg(r=array([  60., 90.,120.,150.,180.,230.,310.,380.,510.,640.,820.]),
               g=array([.195,.075,.14,  .07,  .1,  .04, .06, .05, .02,.015,.01])
               ,alpha=1):
    """
    Generates kappa from the values of the reduced shear g
    Equation 5.1 in Schneider & Seitz 1995

    1.-k(x)=(1-g(x))**-1 * exp(-\int^{r_0}_{r} 2 g(x) /x/(1-g(x)) dx -2 g(r_0)/alpha)

    Can be used if g<1

    """
    #Check whether the values of g are not too big
    if max(g)>=1:
        print "Values of g too big"
        sys.exit()


    #Sort the radial data if needed
    if not ascend(r): r,g=multisort(r,(r,g))


    kappa=zeros(len(r))*1.
    for i in range(len(r)):
        expo=exp(-trapz(2.*g[i:]/r[i:]/(1.-g[i:]),r[i:])-2.*g[-1]/alpha)
        print expo
        kappa[i]=1.-1./(1.-g[i])*expo
    return kappa
    
class lensing_iso:
    def __init__(self,x,y,z_s,z_l=.6,v=200.,cosmo=(0.3,0.7,.7)):
	#Given an isothermal sphere with redshift z_l
	#and velocity dispersion v, estimate the 
	#kappa, shear and magnification at a position x,y
	#from its center
	self.r=sqrt(x**2+y**2)
	self.phi=arctan2(y,x)
	self.r_e=einstein_radius(v,z_l,z_s,cosmo)
	#kappa
	self.kappa=0.5*self.r_e/self.r
	#shear
	self.shear1=self.kappa*cos(2.*self.phi)
	self.shear2=self.kappa*sin(2.*self.phi)
	self.shear=sqrt(self.shear1**2+self.shear2**2)
	self.magnification=1./((1.-self.kappa)**2+self.shear**2)
	self.complex_shear=self.kappa*exp(2j*self.phi)

class multiple_iso:
    def __init__(self,z_s,parameters,cosmo=(0.7,0.3,1.)):
        #Given a set of lenses with positions x,y (arcsec)
        #redshifts z and velocity dispersions velocity
        #estimate the shear and magnification in the origin
        #of coordinates
	#Eliminate objects with z smaller than or equal that z_source
	#OJO!: returns the results ordered by redshift!
	good=less(parameters[2],z_s)
	pars=multicompress(good,parameters)
	#Sort them by redshift
	pars=multisort(pars[2],pars)
	nlens=len(pars[2])
	self.x,self.y,self.z,self.v=pars
	self.r=sqrt(self.x**2+self.y**2) #radius
	self.phi=arctan2(self.y,self.x) #galaxy position angle
	self.r_e=zeros(nlens)*1.
	for i in range(nlens):
	    self.r_e[i]=einstein_radius(self.v[i],self.z[i],z_s,cosmo)
	self.kappa=self.r_e/2./self.r
	self.shear1=self.kappa*cos(2.*self.phi)
	self.shear2=self.kappa*sin(2.*self.phi)
	self.shear=sqrt(self.shear1**2+self.shear2**2)
	self.magnification=1./((1.-self.kappa)**2-self.shear**2)
	#Creating U_i
	u=[]
	for i in range(nlens):	    
	    u.append(zeros((2,2))*1.)
	    u[i][0,0]=self.kappa[i]+self.shear1[i]
	    u[i][0,1]=self.shear2[i]
	    u[i][1,0]=self.shear2[i]
	    u[i][1,1]=self.kappa[i]-self.shear1[i]
#	    print 'u[%i]\n'%i,u[i]
#	    print 'a_single[%i]\n'%i,identity(2)*1.-u[i]


#	print 'U=',u
	    
	#Creating beta
	beta=zeros((nlens,nlens))*1.
	for j in range(nlens):
	    for i in range(j):
		if self.z[i]<>self.z[j]:
		    beta[i,j]=(d_ang(self.z[i],self.z[j],cosmo)*
			       d_ang(0.,z_s,cosmo)/
			       d_ang(0.,self.z[j],cosmo)/
			       d_ang(self.z[i],z_s,cosmo))
		else:
		    beta[i,j]=0.

	#Creating A_i
	a=[]
	for j in range(nlens):
	    a.append(identity(2)*1.)
	    for i in range(j):
		a[j]=a[j]-beta[i,j]*matrixmultiply(u[i],a[i])
#	    print 'a[%i]\n'%j,a[j]
	
	#Creating the final A matrix
	self.a=identity(2)*1.

	for i in range(nlens):
#	    print 'u[%i]\n'%i,u[i]
#	    print 'a[%i]\n'%i,a[i]
#	    print 'u*a',matrixmultiply(u[i],a[i])
	    self.a=self.a-matrixmultiply(u[i],a[i])
#	    print 'i,self.a'
#	    print i,self.a

	self.total_shear1=0.5*(self.a[1,1]-self.a[0,0])#checked
	self.total_shear2=-self.a[0,1]
#	if self.a[0,1]<>self.a[1,0]:
#	    print 'Matrix non symmetric'
#	    print self.a[0,1],self.a[1,0]

	self.total_kappa=1.-0.5*(self.a[0,0]+self.a[1,1])
	self.total_shear=sqrt(self.total_shear1**2+self.total_shear2**2)
	self.total_phi=0.5*arctan2(self.total_shear2,self.total_shear1)
	self.total_magnification=1./((1.-self.total_kappa)**2-self.total_shear**2)

def ab_theta2chi(ab,theta):
    #Given a a/b ratio <1 and a position angle
    #it returns the complex ellipticity used for the class lens2.572
    #(tested)
    return (1.-ab**2)/(1.+ab**2)*exp(2j*theta)

def chi2ab_theta(chi):
    #Given the complex ellipticity chi
    #returns the axis ratio a/b and the position angle
    #(tested)
    ab=sqrt((1.-abs(chi))/(1.+abs(chi)))
    theta=0.5*arctan2(chi.imag,chi.real)
    return ab,theta

class lens:
    def __init__(self,kappa,shear):
	self.kappa=kappa
	self.shear=shear[0]**2+shear[1]**2
	self.shear=complex(shear[0],shear[1])
	self.g=self.shear/(1.-self.kappa)

    def chi(self,chi_s):
	gchi=g*conjugate(chi_s)
	return ((chi_s-2.*self.g+self.g**2*conjugate(chi_s))/
		(1.+abs(g)**2-2.*gchi.real))

    def chi_s(self,chi):
	gchi=g*conjugate(chi)
	return ((chi+2.*self.g+self.g**2*conjugate(chi))/
		(1.+abs(g)**2+2.*gchi.real))  	


class prior_mu_sn:
    #This function estimates the prior magnification distribution
    #For a SN at z=1.7
    #from weak lensing for standard candles based on the results
    # of Wang, Y. 1999
    def __init__(self):
	# We are assuming Omega=0.3,Lambda=0.7
	# Tabulated values from Wang:
	# Omega=0.4,Lambda=0.6
	#    z=1.5  z=2
	d0=(0.43849,0.46109)
	d05=(0.40897,0.4157)
	d1=(0.38069,0.37329)
	d15=(0.35362,0.33352)

	# Interpolation
	d0a=d0[0]+(d0[1]-d0[0])/(2.-1.5)*(1.7-1.5)
	d05a=d05[0]+(d05[1]-d05[0])/(2.-1.5)*(1.7-1.5)
	d1a=d1[0]+(d1[1]-d1[0])/(2.-1.5)*(1.7-1.5)
	d15a=d15[0]+(d15[1]-d15[0])/(2.-1.5)*(1.7-1.5)
	# Omega=0.2,Lambda=0.8
	#    z=1.5  z=2
	d0=(0.49007,0.52031)
	d05=(0.46655,0.48168)
	d1=(0.44374,0.44485)
	d15=(0.42160,0.40976)

	# Interpolation
	d0b=d0[0]+(d0[1]-d0[0])/(2.-1.5)*(1.7-1.5)
	d05b=d05[0]+(d05[1]-d05[0])/(2.-1.5)*(1.7-1.5)
	d1b=d1[0]+(d1[1]-d1[0])/(2.-1.5)*(1.7-1.5)
	d15b=d15[0]+(d15[1]-d15[0])/(2.-1.5)*(1.7-1.5)
	d0=0.5*(d0b+d0a)
	d05=0.5*(d05b+d05a)
	d1=0.5*(d1b+d1a)
	d15=0.5*(d15b+d15a)

	print 'd0,d05,d1,d15'
	print d0,d05,d1,d15

	a=4./3.*(-d0+3.*d05-3.*d1+d15)
	b=2.*(d0-2*d05+d1-3.*a/4.)
	c=-d0+d1-a-b

	print 'a,b,c'
	print a,b,c



	alfa=arange(0.5,3.,.01)

	mu=((d0+a+b+c)/(d0+a*alfa**3+b*alfa**2+c*alfa))**2

	z=1.7
	x=1./(5.*z)
	y=z/5.

	cnorm=0.01*(0.53239+2.79165*y-2.42315*y**2+1.13844*y**3)
	alfa_peak=1.0135-1.07857*x+2.05019*x*x-2.1452*x**3
	w=0.06375+1.75355*x-4.99383*x*x+5.95852*x**3
	q=0.75045+1.85924*y-2.91830*y**2+1.59266*y**3

	pa=cnorm*exp(-((alfa-alfa_peak)/(w*alfa**q))**2)

#	print (3.*a*alfa**2+2.*b*alfa+c)

	pmu=abs(pa/cho*(d0+a+b+c)/2./mu**1.5/(3.*a*alfa**2+2.*b*alfa+c))
	pmu=pmu/add.reduce(pmu)

	p=FramedPlot()
	#	p.add(Curve(log10(mu),log10(pmu)))
        p.xlog=1
        p.ylog=1
	p.add(Curve(mu,pmu))
	p.show()

def test():
    Testing('multiple_iso')
    z=array([.557,.557])
    v=array([238.,1.])
    #1.
    x=array([6.5,-6.5])      
    y=array([-6.5,6.5]) 
    #Unlensed galaxy
    print 'x,y',x[0],y[0],'  ',x[1],y[1]
    ab,theta=.5,pi/4.
    print 'Unlensed galaxy',ab,theta*180./pi
    chi_s=ab_theta2chi(ab,theta)
    z_s=1.66
    sn=multiple_iso(z_s,(x,y,z,v))
    print 'shear1',sum(sn.shear1)
    print 'shear2',sum(sn.shear2)
    print 'total_shear1',sn.total_shear1
    print 'total_shear2',sn.total_shear2
    el=lens(sn.total_kappa,(sn.total_shear1,sn.total_shear2))
    chi=el.chi(chi_s)
    ab_l,theta_l=chi2ab_theta(chi)
    print 'Lensed galaxy',ab_l,theta_l
    if ab_l<ab : print '                  BIEN'
    else : print '                  MAL'

    #2.
    x=array([0.,0.])      
    y=array([-6.5,6.5]) 
    #The lensed galaxy
    print 'x,y',x[0],y[0],'  ',x[1],y[1]
    ab,theta=.5,pi/4.
    print 'Unlensed galaxy',ab,theta*180./pi
    chi_s=ab_theta2chi(ab,theta)
    z_s=1.66
    sn=multiple_iso(z_s,(x,y,z,v))
    print 'shear1',sum(sn.shear1)
    print 'shear2',sum(sn.shear2)
    print 'total_shear1',sn.total_shear1
    print 'total_shear2',sn.total_shear2
    el=lens(sn.total_kappa,(sn.total_shear1,sn.total_shear2))
    chi=el.chi(chi_s)
    ab_l,theta_l=chi2ab_theta(chi)
    print 'Lensed galaxy',ab_l,theta_l
#    if ab_l<ab : print '                  BIEN'
#    else : print '                  MAL'

    #3.
    x=array([6.5,-6.5])      
    y=array([6.5,-6.5]) 
    #The lensed galaxy
    print 'x,y',x[0],y[0],'  ',x[1],y[1]
    ab,theta=.5,pi/4.
    print 'Unlensed galaxy',ab,theta*180./pi
    chi_s=ab_theta2chi(ab,theta)
    z_s=1.66
    sn=multiple_iso(z_s,(x,y,z,v))
    print 'shear1',sum(sn.shear1)
    print 'shear2',sum(sn.shear2)
    print 'total_shear1',sn.total_shear1
    print 'total_shear2',sn.total_shear2
    el=lens(sn.total_kappa,(sn.total_shear1,sn.total_shear2))
    chi=el.chi(chi_s)
    ab_l,theta_l=chi2ab_theta(chi)
    print 'Lensed galaxy',ab_l,theta_l
    if ab_l>ab : print '                  BIEN'
    else : print '                  MAL'

    #4.
    x=array([6.5,-6.5])      
    y=array([0.,0.]) 
    #The lensed galaxy
    print 'x,y',x[0],y[0],'  ',x[1],y[1]
    ab,theta=.5,pi/4.
    print 'Unlensed galaxy',ab,theta*180./pi
    chi_s=ab_theta2chi(ab,theta)
    z_s=1.66
    sn=multiple_iso(z_s,(x,y,z,v))
    print 'shear1',sum(sn.shear1)
    print 'shear2',sum(sn.shear2)
    print 'total_shear1',sn.total_shear1
    print 'total_shear2',sn.total_shear2
    el=lens(sn.total_kappa,(sn.total_shear1,sn.total_shear2))
    chi=el.chi(chi_s)
    ab_l,theta_l=chi2ab_theta(chi)
    print 'Lensed galaxy',ab_l,theta_l
#    if ab_l<ab : print '                  BIEN'
#    else : print '                  MAL'

    kappa=sn.kappa
    shear1=sn.shear1
    shear2=sn.shear2
    shear=sn.shear
    theta=shear1*0.
    for i in range(len(shear1)):theta[i]=arctan2(shear2[i],shear1[i])
    print 'shear1',shear1
    print 'shear2',shear2
    print 'theta',theta
   
    print 'Expected'
    kappa=add.reduce(kappa)
    shear=sqrt(shear[0]**2+shear[1]**2+2.*shear[0]*shear[1]*cos(theta[0]-theta[1]))
#    theta_shear=arctan2(shear[0]*cos(theta[0])+shear[1]*cos(theta[1]),
#	shear[0]*sin(theta[0])+shear[1]*sin(theta[1]))
    magnification=1./((1.-kappa)**2-shear**2)
    print 'kappa,shear,magnification'
    print kappa,shear,magnification

    print 'Calculated'
    print 'kappa,shear,magnification'
    print sn.total_kappa,sn.total_shear,sn.total_magnification



    test='Einstein radius'
    Testing(test)
    z1=0.1
    z2=2.0
    cosmo=(1.,0.,1.)
    print cosmo
    print einstein_radius(186.,z1,z2,cosmo)*d_ang(0.,z2,cosmo)/d_ang(z1,z2,cosmo)
    cosmo=(0.3,0.7,1.)
    print cosmo
    print einstein_radius(186.,z1,z2,cosmo)*d_ang(0.,z2,cosmo)/d_ang(z1,z2,cosmo)

    test='sigma einstein'
    Testing(test)
    theta=einstein_radius(186.,z1,z2,cosmo)
    print 186./sigma_einstein(theta,z1,z2,cosmo)

    test='einstein mass'
    Testing(test)
    m=1e15
    z1=.2
    z2=4.
    theta=4.*g*m*solar_mass/c**2*d_ang(z1,z2,cosmo)/d_ang(0.,z1)/d_ang(0.,z2)
    theta=sqrt(theta/mpc)*180.*3600./pi
    print m/einstein_mass(theta,z1,z2,cosmo)

    print 
    print
    print

    test='d_filled_beam'
    Testing(test)
    print 'Compare with Fig 1a of Fukugita et al. 1992'
    d1=[]
    d2=[]
    xz=arange(0.01,5.,.1)
    for z in xz:
	d1.append(d_filled_beam(0.,z,(1.0,0.0,1.))/cho)
	d2.append(d_filled_beam(0.,z,(0.1,0.9,1.))/cho)
    d1=array(d1)
    d2=array(d2)
    p=FramedPlot()
    p.xrange=0.,5.
    p.yrange=0.,1.
    p.add(Curve(xz,d1))
    p.add(Curve(xz,d2))
    p.show()
    print
    print
    print


    test='d_dyer_roeder'
    Testing(test)
    print 'Compare with Fig 1b of Fukugita et al. 1992'
    d1=[]
    d2=[]
    xz=arange(0.01,5.,.1)
    for z in xz:
	d1.append(d_dyer_roeder(0.,z,(1.0,0.0,1.))/cho)
	d2.append(d_dyer_roeder(0.,z,(0.1,0.9,1.))/cho)
    d1=array(d1)
    d2=array(d2)
    p=FramedPlot()
    p.xrange=0.,5.
    p.yrange=0.,1.
    p.add(Curve(xz,d1))
    p.add(Curve(xz,d2))
    p.show()
    print
    print
    print

if __name__ == '__main__':
    test()
else:
    pass


#def d_dr_old(z_l,z_s,cosmo=(0.3,0.7,1.),nsteps=1000):
#    #Assumes that Omega+lambda=1.
#    #Reference Fukugita et al. 1992,393,3
#    #\alpha=0, that is, all the matter is in clumps
#    if cosmo[0]+cosmo[1] <> 1.:
#        print 'Omega+Lambda must be 1'
#        sys.exit()
#    dz=(z_s-z_l)/float(nsteps)
#    z=arange(z_l,z_s+dz,dz)
#    f=1./sqrt((1.+z)**3*cosmo[0]+1.-cosmo[0])/(1.+z)**2
#    return cho/cosmo[2]*trapz(f,z)*(1.+z_l)    

#def d_fb_old(z_l,z_s,cosmo=(0.3,0.7,1.),nsteps=1000):
#    dz=(z_s-z_l)/float(nsteps)
#    z=arange(z_l,z_s+dz,dz)
#    f=1./sqrt((1.+z)**3*cosmo[0]+1.-cosmo[0])
#    return cho/cosmo[2]*trapz(f,z)/(1.+z_s)    



#def einstein_radius(mass,z_l,z_s,cosmo=(0.3,0.7,1.)):
#    #Assuming Omega=1!
#    d=d_ang(z_l,z_s,cosmo)/d_ang(0.,z_s,cosmo)/d_ang(0.,z_l,cosmo)
#    theta=sqrt(4.*g*mass*d/c**2)*sqrt(solar_mass/mpc)
#    return theta*3600.*180./pi

#def einstein_mass(e_radius,z_l,z_s,cosmo):
#    #Assuming Omega=1!
#    e_radius=e_radius/3600./180.*pi
#    mass=c**2/4./g/d_lensing(z_l,z_s,cosmo)*e_radius**2
#    return mass/solar_mass*mpc

