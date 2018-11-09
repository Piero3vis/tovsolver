import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.constants import G, c, M_sun
from math import pi
from sys import argv
global K, Gamma, rho0_c                             #define globally EoS and central density


def Eul(u , t, dt, rhs):                            #Euler Method
    n = len(u)
    k1 = np.zeros(n)
    up = np.zeros(n)
    k1 = dt*rhs(u,t)
    up = u + k1
    return up


def RK4(u, t, dt ,rhs):                         #Runge-Kutta 4 Method
    n = len(u)
    up = np.zeros(n)
    k1 = np.zeros(n)
    k2 = np.zeros(n)
    k3 = np.zeros(n)
    k4 = np.zeros(n)
    k1 = dt*rhs(u,t)
    k2 = dt*rhs(u + 0.5*k1, t + 0.5*dt)             
    k3 = dt*rhs(u + 0.5*k2, t + 0.5*dt)             
    k4 = dt*rhs(u + k3, t + dt)                     
    up = u + (k1 + 2*k2 + 2*k3 + k4) / 6.0
    return up


def rho0fp(p):                                          #Advice for negative pressure
    if p <   0.0:
        print "Negative p found: %.2e" %(p)
        p = abs(p)
    return (p/K)**(1.0/Gamma) + p/(Gamma - 1.0)         #Energy density
    

def TOVrhs(u,r):
    [m,p] = u
    rho = rho0fp(p)

    #Regularize system at the origin, Taylor series for m
    if r < 1e-3:

        m_rrr = 4.0/3.0*pi * (rho0_c**2/(Gamma - 1.0) + 4.0**(-1.0/Gamma)*(rho0_c)**(2.0/Gamma))
        m_rr = m_rrr * r
        m_r = m_rr * r
    else:
        m_r = m / r
        m_rr = m_r / r
        m_rrr = m_rr / r

    #The real TOV

    MRHS = 4.0*pi*r**2 * rho    #continuity
    PRHS = - m_rr * (rho + p)*(1.0 + 4.0*pi*p / m_rrr) * 1.0/(1.0 - 2.0*m_r)   #pressure equation

    return np.array([MRHS, PRHS])

    #Solver using RK4 or Euler for a given central density

def tovint(rho0_c, dr):
    nsteps = 500000
    sol = np.zeros((nsteps,2))      #store the solution
    n=0
    r=0.0
    p_c = K*rho0_c**Gamma           #Central Pressure
    u = np.array([0.0, p_c])        # [m, p](0,0)
    sol[0] = u
    
    #RK4 METHOD
    
    while  sol[n,1]> 1.0e-12 and  n < nsteps-1:     #P < 1.0e-12 = zeropressure
        u = RK4(u, r , dr, TOVrhs)  
        n += 1
        sol[n] = u
        r+= dr

    '''#<---------------DECOMMENT HERE TO USE EULER METHOD
    #EULER METHOD: 
    
    while sol[n,1]> 1.0e-12 and  n < nsteps-1:
        u = Eul(u, r , dr, TOVrhs)
        n += 1
        sol[n] = u
        r += dr
                  #DECOMMENT HERE TO USE EULER METHOD-------->'''
    
    return (r, sol[:n,0], sol[:n,1], rho0fp(p_c))               #(R, M, P, rho_c)
    
    

    #Plotter for p(r), m(r) and rho(r)

def tovplot(r, m, p):

    #generality of the plot

    rhost = str(rho0_c/1.6199e-18/1e14) #string for title
    
    ri = np.linspace(0.00001,r, len(m))
    fig = plt.figure(1)
    plt.xlim([0.00001,r])
    gs = gridspec.GridSpec(3,3)
    plt.rc('text', usetex= True)
    plt.rc('font', family = 'Iwona ')
    title = plt.suptitle(r'\textsc{Central density} = %s $\times$ 10$^{14}$ g cm$^{-3}$' %rhost, fontsize = 16) 
    plt.subplots_adjust(wspace=0, hspace=0.38)
    


    #presure subplot
    
    ax1 = plt.subplot(gs[0,:])
    ax1.set_xscale("log")
    line1, = ax1.plot(ri, p, 'r--')
    ax1.set_ylabel(r'P(r) [dyn cm$^2$]', fontsize = 14)
    lb1, ub1 = ax1.get_ylim()
    ax1.set_xlim(0.00001,r+5)
    ax1.set_ylim(1, ub1 + 1e34)
    ax1.set_title(r'\textit{Pressure profile}', fontsize = 14)
    
    #mass subplot
    
    ax2 = plt.subplot(gs[1,:])
    ax2.set_xscale("log")
    line2, = ax2.plot(ri, m, 'b-')
    lb2, ub2 = ax2.get_ylim()
    ax2.set_yticks( np.linspace(lb2, ub2, 5))
    plt.ylabel(r'm(r) [M$_{\odot}$]', fontsize = 14)
    ax2.set_ylim(lb2-0.2, ub2+0.2)
    ax2.set_xlim(0.00001,r+5)
    ax2.set_title(r'\textit{Mass profile}', fontsize = 14)

    #density subplot
    
    ax3 = plt.subplot(gs[2,:])
    ax3.set_xscale("log")
    line3, = ax3.plot(ri,rho2, 'g-')
    
    lb3, ub3 = ax3.get_ylim()
    ax3.set_ylim(lb3, ub3 + 1e14)
    ax3.set_ylabel(r'$\rho (r)$ [g cm$^{-3}$]', fontsize = 14)
    ax3.set_xlabel(r'r [km]',fontsize = 14)
    ax3.set_xlim(0.00001,r+5)
    ax3.set_title(r'\textit{Density profile}', fontsize = 14)
    
    
    plt.savefig('plot %s_%s_%s.png' % (rhost,str(m[-1]),str(r)))     #save figure giving central density, mass, radius

                                                        
    plt.show()  


###########################################################################################    

    #MODIFY HERE THE IMPORTANT PARAMETERS!
    
if __name__ == "__main__":

    # MODEL: EOS

    script, rho0_c0, K, Gamma = argv                #CHOOSE CENTRAL DENSITY IN CGS, K IN G=C=Msun=1 units, GAMMA when you run the program
    K = float(K)
    Gamma = float(Gamma)
    rho0_c0 = float(rho0_c0)
                                                        #K = 30000 in  G=C=Msun=1 units corresponds K = 1.98183e-6 in cgs
    #TYPICAL MODEL
    
    #rho0_c0 = 2.2   e14
    #K = 30000
    #Gamma = 2.75
    

    rho0_c = rho0_c0*1.6199e-18             #convert in G=C=Msun=1
    
    M_sun = M_sun.value  
    c = c.value
    G = G.value
    M_solar_G_over_c_sq = M_sun * G/ c**2   #Msun=G=C=1 lenght convert 

    
    dr = 1e-4                               #Stepsize
    
    (r, m , p ,rho_c) = tovint(rho0_c, dr)      #RESOLVE TOV 
    
    global rho1, rho2

    rho1 = (p / K) ** (1. / Gamma) + p/(Gamma - 1.0)
    rho2 = rho1/(1.6199e-18)  
    Press =  (r, m , p ,rho_c)[2]/1.8063e-39
    Prhom = (Press, rho2, m)
    rhoMR = ([(rho0_c0, m[-1], r*M_solar_G_over_c_sq/1000.0)])
    
    print 'stopped at n = %d, r = %.4f km, m = %.4e M_solar and p = %.4e' %(len(m), r*M_solar_G_over_c_sq / 1000.0,m[-1], p[-1])
    print 'Stellar mass: M = %.11e' %m[-1]
    print 'Stellar radius. R = %.11e [km]' %(r*M_solar_G_over_c_sq/1000.0)

    print 'schwartzRadius = %.4e [km]' %(2.**m[-1]*M_solar_G_over_c_sq/1000.0)                                                                        

    tovplot(r*M_solar_G_over_c_sq / 1000.0, m, p/1.8063e-39)            #PLOT PRESSURE, MASS, DENSITY PROFILE



    np.savetxt('prhomRK4%s.txt' %str(rho0_c0), np.transpose(Prhom), fmt= '%.8e', header = 'PRESSURE          DENSITY       MASS      MODEL=%s' %str(rho0_c0)) #save pressure and density and mass profile
    
    #SAVE DENSITY, FINAL MASS, FINAL RADIUS
    
    datafile_path = "PATH:\massradiusgraf.txt"  #instead of PATH include the folder you want to save it
    datafile_id = open(datafile_path, 'a+')
    np.savetxt(datafile_id, rhoMR , fmt='%.11f')
    datafile_id.close()    
    
    
