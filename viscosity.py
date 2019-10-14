import numpy as np
import scipy.optimize as spo
import scipy.integrate as spi
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class property_plot:
    """
    This class must be inherited from a viscosity model class.
      It is not well written and cannot stand alone, which is a problem.
    """

    def __init__(self):
        self.rate_min=.001
        self.rate_max=10000.
        
    def visc_plot(self):
        """
        Method to create log-log plot of viscosity versus shear rate.  
        50 pts are evenly spaced on log axis betwen rate_min and rate_max.
        This class expects to be inherited by a viscosity function class.
        """
        x = np.logspace(np.log10(.001),np.log10(10000.),51)
        y = list( map( lambda rate: self.calc_visc(rate), x))
        plt.loglog(x,y,'-')
        plt.xlabel('Shear rate')
        plt.ylabel('Viscosity')
        plt.title(self.name)
        
    def stress_plot(self):
        """
        Method to create log-log plot of shear stress versus shear rate.  
        50 pts are evenly spaced on log axis betwen rate_min and rate_max.
        This class expects to be inherited by a viscosity function class.
        """
        x = np.logspace(np.log10(.001),np.log10(10000.),51)
        y = list( map( lambda rate: self.calc_visc(rate)*rate, x))
        plt.loglog(x,y,'-')
        plt.xlabel('Shear rate')
        plt.ylabel('Stress')
        plt.title(self.name)    
    

class newtonian(property_plot):
    def __init__(self,name='Default',mu=1.):
        self.name = name
        self.mu = mu
        
    def __str__(self):
        return str(self.name+'\n'+
            'mu ='+str(self.mu)+'\n')
        
    def calc_visc(self,rate):
        return self.mu

class power_law(property_plot):
    def __init__(self,name='Default',k=1.,n=.5):
        self.name = name
        self.k = k
        self.n = n
        
    def __str__(self):
        return str(self.name+'\n'+
            'k ='+str(self.k)+'\n'+
            'n='+str(self.n)+'\n')
        
    def calc_visc(self,rate):
        return self.k*(rate+1.e-9)**(self.n-1.)
    
class carreau(property_plot):
    def __init__(self,name='Default',eta0=10.,etainf=.1,reltime=1.,a=2.,n=.5):
        self.name = name
        self.eta0 = eta0
        self.etainf = etainf
        self.reltime = reltime
        self.a = a
        self.n = n
        
    def __str__(self):
        return str(self.name+'\n'+
            'eta0 ='+str(self.eta0)+'\n'+
            'etainf ='+str(self.etainf)+'\n'+
            'reltime ='+str(self.reltime)+'\n'+
            'a ='+str(self.a)+'\n'+
            'n='+str(self.n)+'\n')
            
    def calc_visc(self,rate):
        return self.etainf + (self.eta0-self.etainf)/(1.0+(self.reltime*rate)**self.a)**((1.-self.n)/self.a)
    
class herschel_bulkley(property_plot):
    """
    Herschel-Bulkley viscosity model using Papanastasiou modification with m=1000.
    Also has eps=1.e-9 with shear rate in denominator of viscosity equation.
    """
    def __init__(self,name='Default',tauy=1.,k=1.,n=1.,m=1000.,m_flag=1):
        self.name = name
        self.tauy = tauy
        self.k = k
        self.n = n
        self.m=m
        self.m_flag=m_flag
        
    def __str__(self):
        return str(self.name+'\n'+
            'tauy ='+str(self.tauy)+'\n'+
            'k ='+str(self.k)+'\n'+
            'n='+str(self.n)+'\n'+
            'm=',str(self.m)+'\n'  )
        
    def calc_visc(self,rate):
        if self.m_flag==0.:
            return self.tauy/(rate+1.e-9) + self.k*(rate+1.e-9)**(self.n-1.)
        if self.m_flag==1:
            return (1.-np.exp(-self.m*rate))*self.tauy/(rate+1.e-9) + self.k*(rate+1.e-9)**(self.n-1.)
    
class three_component(property_plot):
    """
    Marco's 3-component viscosity model using Papanastasiou modification with m=1000.
    Also has eps=1.e-9 with shear rate in denominator of viscosity equation.
    """
    def __init__(self,name='Default',tauy=1.,gamma_crit=1.,eta_bg=1.,m=1000.,m_flag=1):
        self.name = name
        self.tauy = tauy
        self.gamma_crit = gamma_crit
        self.eta_bg = eta_bg
        self.m = m
        self.m_flag=m_flag
        
    def __str__(self):
        return str(self.name+'\n'+
            'tauy ='+str(self.tauy)+'\n'+
            'gamma_crit ='+str(self.gamma_crit)+'\n'+
            'eta_bg='+str(self.eta_bg)+'\n'+
            'm=',str(self.m)+'\n')
        
    def calc_visc(self,rate):
        if self.m_flag==0:
            return self.tauy/(rate+1.e-9) + self.tauy/(rate+1.e-9)*(rate/self.gamma_crit)**0.5 + (self.eta_bg)
        if self.m_flag==1:
            return (1.-np.exp(-self.m*rate))*self.tauy/(rate+1.e-9) + \
                self.tauy/(rate+1.e-9)*(rate/self.gamma_crit)**0.5 + \
                (self.eta_bg)

class bi_power_law(property_plot):
    def __init__(self,name='Default',k_low=1.,n_low=.9,k_high=1.,n_high=.5):
        self.name = name
        self.k_low = k_low
        self.k_high = k_high
        self.n_low = n_low
        self.n_high = n_high
        self.rate_switch = 10.**(np.log10(self.k_high/self.k_low)/(self.n_low-self.n_high))
    
    def __str__(self):
        return str(self.name+'\n'+
                  'k_low ='+str(self.k_low)+'\n'+
                   'k_high ='+str(self.k_high)+'\n'+
                   'n_low ='+str(self.n_low)+'\n'+
                   'n_high ='+str(self.n_high)+'\n'+'\n')
    
    def calc_visc(self,rate):
        eps = 1.e-9
        self.rate_switch = 10.**(np.log10(self.k_high/self.k_low)/(self.n_low-self.n_high))
        if (rate >= self.rate_switch):
            return self.k_high*(rate+eps)**(self.n_high-1.)
        else:
            return self.k_low*(rate+eps)**(self.n_low-1.)

