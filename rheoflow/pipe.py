import numpy as np
import scipy.optimize as spo
import scipy.integrate as spi
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from . import viscosity

class laminar:
    """
    This class contains a variety of methods for computing quantities of interest for laminar flow in a tube.
    The argument viscosity requires a function (or class with method) calc_visc with parameters already set 
    and the shear rate as it's only argument.  Default values are provided.  A default visosity function is provided
    """
    def __init__(self,name='Default',density=1000.,radius=.01,length=1.,viscosity=viscosity.newtonian(name='default',mu=1.), \
             scale=1.e+6,pressure_drop = None, q = None):
        self.name=name
        self.__density = density
        self.__radius=radius
        self.__length=length
        self._viscosity = viscosity
        self.__r = 0.
        self.scale=scale
        if pressure_drop:
            self.pressure_drop = pressure_drop
        else:
            self.pressure_drop = None
        if q:
            self.q = q
            self.__q = q
        else:
            self.__q = None
        
    def __str__(self):
        return str('Name ='+self.name+'\n'+
            'Radius ='+str(self.__radius)+'\n'+
            'Length ='+str(self.__length)+'\n'+
            'Pressure drop ='+str(self.__pressure_drop)+'\n'+
            'Flow rate ='+str(self.__q)+'\n'+
            'Shear rate wall = '+str(self._shear_rate_wall()))
                   
    def _shear_rate_equation(self,solve_rate,dp):
        return solve_rate * self._viscosity.calc_visc(solve_rate) -  \
            dp/self.__length*self.__r/2.

    def shear_rate(self,rad,dp):
        """
        This method computes the shear rate at a radial position (rad) that is the only argument.
        """
        self.__r = rad
        return spo.brentq(lambda x: self._shear_rate_equation(x,dp), 0.,1.e+6)
    
    def vz(self,rad,dp):
        """
        This method computes the axial velocity vz at a radial position, rad.
        """
        # (-) sign deals with vz<0 whehn pressuredrop>0
        #return -spi.quad(self.shear_rate,self.__radius,rad)[0]
        return -spi.quad(lambda x: self.shear_rate(x,dp),self.__radius,rad)[0]
    
    def stress_wall(self):
        """
        Computes shear stress at wall, radial position radius.
        """
        if self.__pressure_drop:
            return self.__radius/2.*self.__pressure_drop/self.__length
        else:
            return None
    
    def _shear_rate_wall(self):
        """
        Computes the true wall shear rate, or shear rate at radial position radius.
        """
        rad = self.__radius
        if self.__pressure_drop:
            dp = self.__pressure_drop
            return self.shear_rate(rad,dp)
        else:
            return None
    
    def shear_rate_plot(self):
        """
        Creates plot of shear rate versus radial position.
        """
        dp = self.__pressure_drop
        x = np.linspace(0.,self.__radius,51)
        y = [self.shear_rate(rad,dp) for rad in x]
        plt.plot(x,y)
        plt.xlabel('Radial position')
        plt.ylabel('Shear rate')
        
    def viscosity_wall(self):
        """
        Computes viscosity at wall, radial position radius.
        """
        return self._viscosity.calc_visc(self._shear_rate_wall())
    
    def re_wall(self):
        """
        Computes Reynolds number at the wall.
        """
        return self.__density*self.__radius*2.*self.__q/(3.14159*self.__radius**2)/self.viscosity_wall()
        
    def vz_plot(self):
        """
        Creates plot of axial velocity versus radial position.
        """
        dp = self.__pressure_drop
        x = np.linspace(0.,self.__radius,51)
        y = [self.vz(rad,dp) for rad in x]
        plt.plot(x,y)
        plt.xlabel('Radial position')
        plt.ylabel('Velocity')
    
    def __q_integrand(self,rad,dp):
        return 2.*3.141592653589*rad*self.vz(rad,dp)
    
    def __q_calc(self,dp):
        """
        Computes volumetric flow rate for preessure drop argument of dp_want.
        The object attribute self.pressure_drop is set to dp_want.
        """
        #return spi.quad(self.__q_integrand(),0.,self.__radius)[0]
        return spi.quad( lambda x: self.__q_integrand(x,dp),0.,self.__radius)[0]
    
    def __q_eqn(self,dp_v):
        return self.__q_calc(dp_v) - self.__q
    
    def __dp_calc(self):
        """
        Computes the pressure drop for a volumetric flow rate of q_want.
        The computation is iterative due to nature of many viscosity functiions.
        The object attribute self.pressure_drop is set to result.
        """
        self.__pressure_drop= spo.brentq(self.__q_eqn,self.scale_min,self.scale_max)
        return
    
    def q_plot(self,pressure_drop_min,pressure_drop_max):
        """
        Creates log-log plot of pressure drop versus flow rate.
        A log-spacing of pressure drops between args pressure_drop_min and pressure_drop_max
        are created.
        """
        x = np.logspace(np.log10(pressure_drop_min),np.log10(pressure_drop_max),51)
        y = [self.__q_calc(dp) for dp in x]
        plt.loglog(y,x,'-')
        plt.xlabel('Flow rate')
        plt.ylabel('Pressure drop')
        plt.title(self.name)

    @property
    def pressure_drop(self):
        return self.__pressure_drop

    @pressure_drop.setter
    def pressure_drop(self,pressure_drop):
        if pressure_drop:
            self.__pressure_drop = pressure_drop
            self.__q = self.__q_calc(pressure_drop)
        else:
            self.__pressure_drop = None

    @property
    def q(self):
        return self.__q

    @q.setter
    def q(self,q):
        if q:
            self.__q = q
            # Estimate dp_a to set scale to reasonalble value
            dp_a = 8.*self._viscosity.calc_visc(4.*self.__q/(3.14158*self.__radius**3))* \
                self.__length*self.__q/(3.14159*self.__radius**4)
            self.scale_max=2.*dp_a
            self.scale_min=dp_a/100.
            self.__dp_calc()
        else:
            self.__q = None

    @property
    def density(self):
        return self.__density

    @property
    def shear_rate_wall(self):
        return self._shear_rate_wall()

    @property
    def shear_stress_wall(self):
        return self.stress_wall()

    @property
    def radius(self):
        return self.__radius

    @radius.setter
    def radius(self,radius):
        self.__radius = radius
        if self.__pressure_drop:
            self.__q = self.__q_calc(self.__pressure_drop)
        else:
            self.__pressure_drop = None

    @property
    def length(self):
        return self.__length

    @length.setter
    def length(self,length):
        self.__length = length
        if self.__pressure_drop:
            self.__q = self.__q_calc(self.__pressure_drop)
        else:
            self.__pressure_drop = None
            
    #@property
    #def viscosity(self):
    #    return self._viscosity
    
    #@viscosity.setter
    #def viscosity(self,viscosity):
    #    self._viscosity = viscosity
    #    if self.pressure_drop:
    #        self.__q = self.__q_calc(self.pressure_drop)
    #    else:
    #        self.__pressure_drop = None


    
