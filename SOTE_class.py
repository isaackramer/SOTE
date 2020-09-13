"""Class used to run the system. Includes integration of ds/dt term.
Includes soil degradation related to declines in hydraulic conductivity."""

import numpy as np
import sympy as sym
import soil_properties

class SOTE(object):
    def __init__(self, Cirr=0.0, Eirr=0.0, C_init=0.0, E_init=0.0, ET_w = 0.0,
                 s_init = 0.0, depth = 0, soil_type = "class_1", Crain = 0.0,
                 Erain = 0.0, rain_prob = 0.0, days = 0.0, mean_height = 0.0,
                 dt=0.0, ET_ratio = 0.0, runs = 100):
        self.par = {}
        self.par.update(Cirr=Cirr, Eirr=Eirr, C_init=C_init,
                        E_init=E_init, ET_w=ET_w,
                        s_init = s_init, dt=dt, Crain = Crain,
                        Erain = Erain, rain_prob = rain_prob, days = days,
                        mean_height = mean_height, ET_ratio = ET_ratio,
                        soil_type = soil_type, runs = runs)
        # create new dictionary with soil properties
        self.soil_parms = soil_properties.soil_props(soil_type = soil_type,
                                                     depth = depth)

        # initial values for saturated hydrualic conductivity and water dynamics
        self.Ksat = self.soil_parms['Ks']
        self.soil_parms['nZr'] = self.soil_parms['n']*self.soil_parms['Zr']
        self.functions_def()

    def functions_def(self):
        # define parts needed to calculate dE/dt
        self.K1 = lambda q, E, w: (self.dgsdC(q/w, E) * self.dsdt * np.power(q,2))*(self.soil_parms['nZr'])
        self.K2 = lambda q, E, w: self.dgsdC(q/w, E) * q * self.dqdt * w
        self.K3 = lambda w: self.Cinput * self.Einput * self.s_In_rate * np.power(w,2)
        self.K4 = lambda q, E, w: self.Lw * q * w * self.gs(q/w, E)
        self.K5 = lambda q, E, w: self.dqdt * np.power(w,2) * self.gs(q/w, E)
        self.K6 = lambda w: self.soil_parms['CEC'] * self.soil_parms['Msoil'] * np.power(w,2)
        self.K7 = lambda q, E, w: self.dgsdE(q/w, E) * q * np.power(w,2)

        # define Gapon equation and partial derivatives
        self.set_g()

        # equations needed for relative hydralulic conductivity, based on Ezlit
        # Eq. [11], x0 is the adjusted effective swelling factor, f is empirical
        self.x0 = lambda E, C: self.soil_parms['p7'] * self.ESPstar(E,C) * self.dstar(C) * 3.6e-4
        # Eq. [12], adjusted ESP
        self.ESPstar = lambda E, C: 100.0*E - (self.soil_parms['p5']+self.soil_parms['p6']*np.log(C))
        # Eq. [8], adjusted interlayer spacing
        self.dstar = lambda C: np.where(C < 300, abs(356.4 / np.sqrt(C) - 20.58), 0)
        # Eq. [13], c and n are constants for a given soil within a specified range of ESP
        self.n = lambda E: np.power(E, self.soil_parms['p1']) + self.soil_parms['p2']
        # Eq. [14]
        self.c = lambda E: self.soil_parms['p3'] * np.exp(self.soil_parms['p4']*E )
        # reduction function itself
        self.cx0n = lambda E, C: self.c(E) * np.power(self.x0(E,C), self.n(E))
        self.R = lambda E, C: 1.0 - self.cx0n(E,C) / (1.0+self.cx0n(E,C))

    def derivs(self):
        # define dq/dt, dE/dt, and ds/dt
        self.dEdt = lambda q, E, w: (self.K1 (q, E, w) - self.K2(q, E, w) + self.K3(w) - self.K4(q, E, w) - self.K5(q, E, w))/(self.K6(w) + self.K7(q, E, w))
        self.dqdt = self.input_salts_rate - self.output_salts_rate
        self.dsdt = (self.s_In_rate - self.s_Out_rate)/(self.soil_parms['n']*self.soil_parms['Zr'])

        # put all derivatives into a single function
        self.rhs = lambda data: np.array([self.dqdt,
                                          self.dEdt(data[0], data[1], data[2]*self.soil_parms['n']*self.soil_parms['Zr']),
                                          self.dsdt])

    def set_g(self):
        # define gapon exchange and partials
        C, E = sym.symbols('E C')
        gs_eq = lambda C, E: 2.0 / (1.0 + (1.0 + 8.0 * (self.soil_parms['Kg']) ** 2 * C * (1.0 / E - 1.0) ** 2) ** 0.5)
        self.gs = lambda C, E: gs_eq(C, E)
        diffgsdE = sym.diff(self.gs(C, E), E)
        diffgsdC = sym.diff(self.gs(C, E), C)
        self.dgsdE = sym.lambdify((C, E), diffgsdE, "numpy")
        self.dgsdC = sym.lambdify((C, E), diffgsdC, "numpy")

    def set_K(self, C, E):
        # calculate the relative hydraulic conductivity
        X0 = self.x0(E, C)
        r_ezlit = np.where(X0<=0, 1.0, self.R(E, C)) # gets rid of negative values
        return r_ezlit

    def input_q(self):
        # calculate the salinity of the input water
        irr_salts_rate = self.irrigation_rate*self.par['Cirr']
        rain_salts_rate = self.rain_rate*self.par['Crain']
        self.input_salts_rate = irr_salts_rate + rain_salts_rate
        with np.errstate(divide='ignore', invalid='ignore'):
            self.Cinput = np.where(self.s_In_rate == 0, 0, self.input_salts_rate/(self.s_In_rate))

    def input_sodicity(self):
        # calculate the sodicity of the input water
        irrigation_sodicity = self.irrigation_rate*self.par['Eirr']
        rain_sodicity = self.rain_rate*self.par['Erain']
        self.input_sodium = irrigation_sodicity + rain_sodicity
        with np.errstate(divide='ignore', invalid='ignore'):
            self.Einput = np.where(self.s_In_rate == 0, 0, self.input_sodium/(self.s_In_rate))

    def water_input(self):
        # calculate total input water and salinity and sodicity of input water
        self.irrigation_rate = self.Irr
        self.s_In_rate = self.irrigation_rate + self.rain_rate
        self.input_q()
        self.input_sodicity()

    def rain_height(self):
        # calculate rain height based on exponential distribution
        rand = np.random.rand(self.par['runs'],)
        event_height = np.where(rand <= self.par['rain_prob']*self.par['dt'],
                                np.random.exponential(self.par['mean_height'], size = self.par['runs']), 0)
        self.rain_rate = event_height/self.par['dt']
        # account for infiltration excess runoff
        self.rain_rate = np.where(self.rain_rate>self.Ksat, self.Ksat, self.rain_rate)

    def water_loss(self, s, C):
        # calculate water loss based on soil water content/properties and ETmax
        self.ET_act = np.piecewise(s,
                                   [s <= self.soil_parms['s_h'],
                                    (s <= self.soil_parms['s_w']) & (s > self.soil_parms['s_h']),
                                    (s <= self.soil_parms['s_bal']) & (s > self.soil_parms['s_w'])],
                                    [lambda s: 0,
                                     lambda s: self.par['ET_w']*((s - self.soil_parms['s_h'])/(self.soil_parms['s_w'] - self.soil_parms['s_h'])),
                                     lambda s: self.par['ET_w'] + (self.ETmax - self.par['ET_w']) *((s - self.soil_parms['s_w'])/(self.soil_parms['s_bal'] - self.soil_parms['s_w'])),
                                     lambda s: self.ETmax])

        # calculate leaching
        self.Lw = self.Ksat*np.power(s,(self.soil_parms['c']))

        # calculate total loss rate
        self.s_Out_rate = self.ET_act + self.Lw

        # calculate the loss of salts
        self.output_salts_rate = (self.Lw)*C

    def water_net(self, s, C, E):
        # net change in water content
        self.water_input()
        self.water_loss(s, C)
        self.Ksat = self.set_K(C, E)*self.soil_parms['Ks']

    def rk4_step(self, data):
        # Integrate using Runge-Kutta
        dt = self.par['dt']
        k1 = dt * self.rhs(data)
        k2 = dt * self.rhs(data + 0.5 * k1)
        k3 = dt * self.rhs(data + 0.5 * k2)
        k4 = dt * self.rhs(data + k3)
        return np.array([data + (k1 + 2.0 * (k2 + k3) + k4) / 6.0])
