import numpy as np
import scipy.optimize as spo
import scipy.integrate as spi
from scipy.integrate import odeint
import matplotlib.pyplot as plt


    # re_wall, and stuff? ow to access, vz -> change vz to vz_calc

    

class laminar_slit_flow:
    """
    This class contains a variety of methods for computing quantities of interest for laminar flow in a slit.
    The argument viscosity requires a function (or class with method) calc_visc with parameters already set 
    and the shear rate as it's only argument.
    """
    def __init__(self,name='default',height=0.01,width=0.1,length=1.,density=1000., \
        pressure_drop = None, q = None, viscosity=lambda x: 1.0):
        self.name=name
        self.__density = density
        # document 1/2H
        self.__height = height/2.
        self.__width = width
        self.__length = length
        self._viscosity = viscosity
        self.__y = 0.
        if pressure_drop:
            self.pressure_drop = pressure_drop
        else:
            self.pressure_drop = None
        if q:
            self.q = q
        else:
            self.q = None

        
    def __str__(self):
        h = 2.*self.__height
        return str('Name ='+self.name+'\n'+
            'Height ='+str(h)+'\n'+
            'Width ='+str(self.__width)+'\n'+
            'Length ='+str(self.__length)+'\n'+
            'Pressure drop ='+str(self.__pressure_drop)+'\n'+
            'Flow rate ='+str(self.__q)+'\n'+
            'Shear rate wall = '+str(self.shear_rate_wall()))
    
    def _shear_rate_equation(self,solve_rate,dp):
        return solve_rate * self._viscosity.calc_visc(solve_rate) -  \
            dp/self.__length*self.__y

    def shear_rate(self,h,dp):
        """
        This method computes the shear rate at a y position for dp.
        """
        self.__y = h
        return spo.brentq(lambda x: self._shear_rate_equation(x,dp), 0.,1.e+9)

    def shear_rate_wall(self):
        """
        Computes the true wall shear rate, or shear rate at radial position radius.
        """
        height = self.__height
        if self.__pressure_drop:
            dp = self.__pressure_drop
            return self.shear_rate(self.__height,dp)
        else:
            return None
        
    
    def vz(self,h,dp):
        """
        This method computes the axial velocity vz at a radial position, rad.
        """
        # (-) sign deals with vz<0 whehn pressuredrop>0
        #return -spi.quad(self.shear_rate,self.__radius,rad)[0]
        return -spi.quad(lambda x: self.shear_rate(x,dp),self.__height,h)[0]

    def stress_wall(self):
        """
        Computes shear stress at wall, radial position radius.
        """
        if self.__pressure_drop:
            return self.__heigth/2.*self.__pressure_drop/self.__length
        else:
            return None

    def __q_integrand(self,h,dp):
        return 2.*self.__width*self.vz(h,dp)

    def __q_calc(self,dp):
        """
        Computes volumetric flow rate for preessure drop argument of dp_want.
        The object attribute self.pressure_drop is set to dp_want.
        """
        #return spi.quad(self.__q_integrand(),0.,self.__radius)[0]
        return spi.quad( lambda x: self.__q_integrand(x,dp),0.,self.__height)[0]
    
    def __q_eqn(self,dp_v):
        return self.__q_calc(dp_v) - self.__q
    
    def __dp_calc(self):
        """
        Computes the pressure drop for a volumetric flow rate of q_want.
        The computation is iterative due to nature of many viscosity functiions.
        The object attribute self.pressure_drop is set to result.
        """
        comp_pressure_drop = spo.brentq(self.__q_eqn,0.,1.e+6)
        self.__pressure_drop = comp_pressure_drop
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

    def shear_rate_plot(self):
        """
        Creates plot of shear rate versus radial position.
        """
        dp = self.__pressure_drop
        x = np.linspace(0.,self.__height,51)
        y = [self.shear_rate(h,dp) for h in x]
        plt.plot(x,y)
        plt.xlabel('Height')
        plt.ylabel('Shear rate')
        
    def viscosity_wall(self):
        """
        Computes viscosity at wall, radial position radius.
        """
        return self._viscosity.calc_visc(self.shear_rate_wall())
    
    def re_wall(self):
        """
        Computes Reynolds number at the wall.
        """
        return self.__density*self.__height*2.*self.__q/(self.__width*2.*self.__height)/self.viscosity_wall()
        
    def vz_plot(self):
        """
        Creates plot of axial velocity versus radial position.
        """
        dp = self.__pressure_drop
        x = np.linspace(0.,self.__height,51)
        y = [self.vz(height,dp) for height in x]
        plt.plot(x,y)
        plt.xlabel('Y position position')
        plt.ylabel('Velocity')

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
            self.__dp_calc()
        else:
            self.__q = None


    
                   