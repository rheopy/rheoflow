import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import scipy.optimize as spo
from scipy.optimize import fsolve


class friction_factor:
    """
    This class computes pipe flow information based on the non-Newtonian Dodge-Metzner paper.

    There is a Jupyter notebook demonstrating usage.
    """
    def __init__(self,name,rho,d,l,viscosity):
        self.name=name
        self.__rho=rho
        self.__d=d
        self.__l=l
        self._viscosity=viscosity # This is the viscosity function 
        self._friction = self._f_dm # Default friction factor for now
        self.__u = None
        self.__pressure_drop = None
        self.__f = None
        self.__re = None
        self.gammadotw = None
        self.tauw = None
    
    def __str__(self):
        return str('Name= '+self.name+'\n'+
            'Diameter = '+str(self.__d)+'\n'+
            'Length = '+str(self.__l)+'\n'+
            'Density = '+str(self.__rho)+'\n'+
            'U = '+str(self.u)+'\n'+
            'Pressure drop = '+str(self.pressure_drop)+'\n'+
            'Friction factor = '+str(self.__f)+'\n'+
            'Reynolds number = '+ str(self.__re)+'\n'+
            'Wall shear rate = '+str(self.gammadotw)+'\n'+
            'Wall shear stress = '+str(self.tauw)+'\n'
            )
    
    def _f_dm(self,re,tauw):
        """
        _f_fm returns the Fanning friction factor given Re (re) and wall stress (tauw).
        It is based on Dodge-Metzner paper.
        """
        # Step 1 - compute wall shear rate, gammadot_f
        gammadot_f = spo.brentq(lambda x: x*self._viscosity(x)-tauw,0.,1.e+9)
        gammadot_f = np.abs(np.real(gammadot_f))
        # Step 2 - compute n' = dlog(tauw)/dlog(gammadotw).  Simple central FD
        dx = .01 # arbitrary for now
        #nprime = (np.log10(self.viscosity(gammadot_f+dx)*(gammadot_f+dx))- \
        #         np.log10(self.viscosity(gammadot_f-dx)*(gammadot_f-dx)))/        \
        #        (np.log10(gammadot_f+dx)-np.log10(gammadot_f-dx))
        nprime = (np.log10(self._viscosity(gammadot_f+dx)*(gammadot_f+dx))- \
                 np.log10(self._viscosity(gammadot_f)*(gammadot_f)))/        \
                (np.log10(gammadot_f+dx)-np.log10(gammadot_f))
        if (nprime<0.):
            nprime=.01
        elif (nprime>1.):
            nprime=1.

        # Step 3 - comnpute laminar f and use f=0.008 as cutoff between laminar and turbulent.
        f_laminar = 16./(np.abs(re)+1.0e-9)
        if (f_laminar < 0.008):   # Solve turbulent case
            f_dodgemetz = lambda x: np.sqrt(1.0/(x+1.e-9)) - \
                    4.0/nprime**0.75*np.log10(np.abs(re*(x+1.e-9)**(1.-nprime/2.))) + 0.4/nprime**1.2
            f_fanning = spo.brentq( f_dodgemetz, 0., 1.e+6)
        else:
            f_fanning = f_laminar
        return f_fanning

    def _equations_u(self,u,p):
        """
        Simultanious equations for pipe flow when U is input and delta P must be calculated.
        There are four equations with tauw, re, dp, and gammadotw as unknowns.  Equations are equal to 0.
        Eqs
        0 - re = rho d U / viscosity(gammadotw)
        1 - tauw = viscosity(gammadotw)*gammadotw
        2 - f(re,tauw) = delta P * D / (2 * rho U^2 L)
        3 - tauw = D/4 Delta P / L
        """
        [tauw,re,dp,gammadotw] = p
        eqs = [re-self.__rho*self.__d*u/self._viscosity(gammadotw), \
               tauw-self._viscosity(gammadotw)*gammadotw, \
                self._friction(re,tauw)-dp*self.__d/(2.*self.__rho*u**2*self.__l), \
               tauw-self.__d/4.*dp/self.__l]
        return eqs

    def _equations_dp(self,dp,tauw,gammadotw,p):
        """
        Simultanious equations for pipe flow when delta P is input and U must be calculated.
        There are two equations with u and re as unknowns.  Equations are equal to 0. tauw and 
        gammadotw are computed directly and provided as arguments.
        Eqs
        0 - re = rho D U / viscosity(gammadotw)
        1 - f(re,tauw) = delta P * D / (2 * rho U^2 L)
        """
        [re,u] = p
        eqs = [re*self._viscosity(gammadotw)-self.__rho*self.__d*u,   \
            self._friction(re,tauw)*2.*self.__rho*u**2*self.__l-dp*self.__d]
        return eqs

    def __pipe_u(self):
        """
        """
        
        # viscosity and friction are functions viscosity(rate), friction(re,tauw,viscosity)
        # Calc apparent wall shear rate for guesses
        gammadot_a = 8.*self.__u/self.__d
        # Calculate tauw guess via app shear rate
        tau_guess = self._viscosity(gammadot_a)*gammadot_a
        re_guess = self.__rho*self.__d*self.__u/self._viscosity(gammadot_a)
        dp_guess = tau_guess*4.*self.__l/self.__d
        
        # Use different guesses for laminar or turbulent cases to help convergence.
        # guess is list of initial guesses
        if (re_guess<2000.):
            guess = [1.*tau_guess,1.*re_guess,1.*dp_guess,1.*gammadot_a]
        else:
            guess = [.1*tau_guess,1.*re_guess,.1*dp_guess,.5*gammadot_a]

        # Solve for delta P (dp), U, Re, and shear rate (gammadotw). p is list of variables.
        # guess is list of initial guesses.
        ans = spo.fsolve( lambda p: self._equations_u(self.__u,p),guess)

        self.dp=ans[2]
        self.__pressure_drop = ans[2]
        self.__re=ans[1]
        self.tauw=ans[0]
        self.gammadotw = ans[3]
        self.__f = self._friction(self.__re,self.tauw)
        return

    def __pipe_dp(self):
        """
        """
        # Compute wall stress tauw
        tauw_calc = self.__d/4.*self.__dp_target/self.__l
        # Compute wall shear rate gammadot
        gammadot_calc = spo.brentq(lambda x: x-tauw_calc/self._viscosity(x),0.,1.e+9)
        # u guess - needs to be good for high re
        u_guess=self.__d/8.*gammadot_calc
        # re guess - needs to be good for high re
        re_guess = self.__rho*self.__d*u_guess/self._viscosity(gammadot_calc)
        if (re_guess<2000.):
            guess = [re_guess,u_guess]
        elif (re_guess < 5000.):
            guess = [re_guess*.1,u_guess*.01]
        else:
            guess = [re_guess*.1,u_guess*.01]
        ans = spo.fsolve( lambda p: self._equations_dp(self.__dp_target,tauw_calc,gammadot_calc,p), \
                        guess,maxfev=100000)

        #self.pressure_drop = self.__dp_target
        #self.__pressure_drop = self.dp_target
        self.u = ans[1]
        #self.__u = self.u
        self.__re=ans[0]
        self.tauw=tauw_calc
        self.gammadotw = gammadot_calc
        self.__f = self._friction(self.__re,self.tauw)
        self.__pressure_drop = self.__dp_target
        #print(self)
        return
    

    @property
    def pressure_drop(self):
        return self.__pressure_drop

    @pressure_drop.setter
    def pressure_drop(self,pressure_drop):
        self.__dp_target = pressure_drop
        self.__pipe_dp()

    @property
    def u(self):
        return self.__u 
    
    @u.setter
    def u(self,u):
        self.__u = u
        self.__pipe_u()

    @property
    def f(self):
        return self.__f

    @property
    def re(self):
        return self.__re

    # Decide later whether diameter and length are mutable
    # For now, choose U to dribve if modified
    @property
    def d(self):
        return self.__d
    
    @d.setter
    def d(self,d):
        self.__d = d
        self.__pipe_u()

    @property
    def l(self):
        return self.__l

    @l.setter
    def l(self,l):
        self.__l = l
        self.__pipe_u()

    @property
    def rho(self,rho):
        return self.__rho

    @rho.setter
    def rho(self,rho):
        self.__rho = rho
        self.__pipe_u()








       
