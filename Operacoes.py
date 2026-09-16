from . import config
from .utils import np, pd, minimize, brentq

class Port:
    def __init__(self, nome, entrada, corrente = None):
        self.Nome = nome
        self.Tipo = entrada # 0 = entrada, 1 = saida
        self.Corrente = corrente

class Operacao:
    __slots__ = ['Tag',
                 'FP',
                 'Ports'
                ]

    def __init__(self, nome, fp):
        self.Tag = nome
        self.FP = fp

class Mixer(Operacao):
    
    __slots__ = Operacao.__slots__

    def __init__(self, nome, fp):
        super().__init__(nome, fp)
        self.Ports = {'C_IN': Port('C_IN', 0),
                      'C_OUT': Port('C_OUT', 1)
                     }

    def Add_Corrente(self, corr, port):
        if port == 'C_IN':
            if self.Ports[port].Corrente is None:
                self.Ports[port].Corrente = corr
            else:
                i = 1
                while f'C_IN_{i}' in self.Ports.keys():
                    i += 1
                self.Ports[f'C_IN_{i}'] = Port(f'C_IN_{i}', 0)
                self.Ports[f'C_IN_{i}'].Corrente = corr
        elif port == 'C_OUT':
            self.Ports[port].Corrente = corr
        self.Calc_Op()

    def Checa_Atributos(self):
        for p in self.Ports.values():
            if p.Corrente is None:
                return False
            elif not p.Corrente.Checa_Atributos() and p.Tipo == 0:
                return False
        return True        

    def Calc_Op(self):
        if self.Checa_Atributos():
            corr_in = [v.Corrente for k,v in self.Ports.items() if 'IN' in k]
            corr_out = self.Ports['C_OUT'].Corrente
            p_out = min([c.p for c in corr_in])
            n = np.array([c.z*c.F for c in corr_in])
            n = np.sum(n, axis = 0)
            F = np.array([c.F for c in corr_in])
            z_out = n/np.sum(n)
            F_out = np.sum(F)
            corr_out.F = F_out
            corr_out.z = z_out
            corr_out.p = p_out
            T0 = np.sum(F*np.array([c.T for c in corr_in]))/np.sum(F)
            def Loop(T):
                corr_out.T = T
                return corr_out.HF - np.sum([c.HF for c in corr_in])
            T = brentq(Loop, T0 - 50, T0 + 50)
            corr_out.T = T

class Compressor(Operacao):

    __slots__ = Operacao.__slots__ + ['_Compressor__pout', '_Compressor__Tout',
                                      '_Compressor__wisent', '_Compressor__eff']

    def __init__(self, nome, fp):
        super().__init__(nome, fp)
        self.Ports = {'C_IN': Port('C_IN', 0),
                      'C_OUT': Port('C_OUT', 1)
                     }
        self.__pout = None
        self.__Tout = None
        self.__wisent = None
        self.Eff = 0.75

    def Add_Corrente(self, corr, port):
        if port in self.Ports.keys():
            self.Ports[port].Corrente = corr
        if (self.Ports['C_IN'].Corrente is not None and self.Ports['C_OUT'].Corrente is not None):
            self.Ports['C_OUT'].Corrente.z = self.Ports['C_IN'].Corrente.z
            self.Ports['C_OUT'].Corrente.F = self.Ports['C_IN'].Corrente.F
            if self.Ports['C_OUT'].Corrente.p != None:
                self.P_out = self.Ports['C_OUT'].Corrente.p
        self.Calc_Op()

    @property
    def P_in(self):
        return self.Ports['C_IN'].Corrente.p
    @P_in.setter
    def P_in(self, valor):
        self.Ports['C_IN'].Corrente.p = valor
        self.Calc_Op()

    @property
    def T_in(self):
        return self.Ports['C_IN'].Corrente.T
    @T_in.setter
    def T_in(self, valor):
        self.Ports['C_IN'].Corrente.T = valor
        self.Calc_Op()

    @property
    def P_out(self):
        return self.__pout if self.__pout is not None else self.P_in
    @P_out.setter
    def P_out(self, valor):
        self.__pout = valor
        self.Calc_Op()

    @property
    def T_out(self):
        return self.__Tout if self.__Tout is not None else self.T_in
    @T_out.setter
    def T_out(self, valor):
        pass

    @property
    def W_isent(self):
        return self.__wisent
    @W_isent.setter
    def W_isent(self, valor):
        pass

    @property
    def Eff(self):
        return self.__eff
    @Eff.setter
    def Eff(self, valor):
        if all([valor >= 0, valor <= 1]):
            self.__eff = valor
            self.Calc_Op()

    @property
    def W_real(self):
        return self.__wisent/self.__eff
    @W_real.setter
    def W_real(self, valor):
        pass

    def Checa_Atributos(self):
        for p in self.Ports.values():
            if p.Corrente is None:
                return False
            elif not p.Corrente.Checa_Atributos() and p.Tipo == 0:
                return False
        if self.P_out is None:
            return False
        else:
            return True

    def Calc_Op(self):
        if self.Checa_Atributos():
            c_in = self.Ports['C_IN'].Corrente
            gamma = c_in.Cp/(c_in.Cp - config.R)
            T_ideal = self.T_in*(self.P_out/self.P_in)**((gamma - 1)/gamma)
            self.__Tout = self.T_in + (1/self.__eff)*(T_ideal - self.T_in)
            self.__wisent = self.Ports['C_IN'].Corrente.F*self.Ports['C_IN'].Corrente.Cp*(T_ideal - self.T_in)
            self.Ports['C_OUT'].Corrente.T = self.__Tout
            self.Ports['C_OUT'].Corrente.p = self.__pout

    @property
    def Output(self):
        output = {'W': self.W_real, 'T_saida': self.T_out}
        return pd.DataFrame(output, index = [self.Tag])

class TrocadordeCalor(Operacao):
    
    __slots__ = Operacao.__slots__ + ['_TrocadordeCalor__dp',
                                      '_TrocadordeCalor__Tout',
                                      '_TrocadordeCalor__q']

    def __init__(self, nome, fp):
        super().__init__(nome, fp)
        self.Ports = {'C_IN': Port('C_IN', 0),
                      'C_OUT': Port('C_OUT', 1)}
        self.__dp = 0
        self.__Tout = None
        self.__q = 0

    def Add_Corrente(self, corr, port):
        if port in self.Ports.keys():
            self.Ports[port].Corrente = corr
        self.Calc_Op()

    @property
    def T_in(self):
        if self.Checa_Atributos():
          return self.Ports['C_IN'].Corrente.T
    @T_in.setter
    def T_in(self, valor):
        self.Ports['C_IN'].Corrente.T = valor
        self.Calc_Op()

    @property
    def p_in(self):
        return self.Ports['C_IN'].Corrente.p
    @p_in.setter
    def p_in(self, valor):
        self.Ports['C_IN'].Corrente.p = valor
        self.Calc_Op()

    @property
    def T_out(self):
        return self.__Tout if self.__Tout is not None else self.T_in
    @T_out.setter
    def T_out(self, valor):
        self.__Tout = valor
        self.Calc_Op()

    @property
    def dP(self):
        return self.__dp
    @dP.setter
    def dP(self, valor):
        self.__dp = valor
        self.Calc_Op()

    @property
    def Q(self):
        return self.__q
    @Q.setter
    def Q(self, valor):
        pass

    def Checa_Atributos(self):
        for p in self.Ports.values():
            if p.Corrente is None:
                return False
            elif not p.Corrente.Checa_Atributos() and p.Tipo == 0:
                return False
        if self.dP is None:
            return False
        else:
            return True

    def Calc_Op(self):
        if self.Checa_Atributos():
            z = self.Ports['C_IN'].Corrente.z
            self.Ports['C_OUT'].Corrente.z = self.Ports['C_IN'].Corrente.z
            self.Ports['C_OUT'].Corrente.F = self.Ports['C_IN'].Corrente.F
            self.Ports['C_OUT'].Corrente.p = self.p_in - self.__dp
            self.Ports['C_OUT'].Corrente.T = self.T_out
            self.__q = self.Ports['C_OUT'].Corrente.HF - self.Ports['C_IN'].Corrente.HF
        else:
            self.__q = 0

    @property
    def Output(self):
        output = {'Q': self.Q}
        return pd.DataFrame(output, index = [self.Tag])

class VasoFlash(Operacao):

    def __init__(self, nome, fp):
        super().__init__(nome, fp)
        self.Ports = {'C_IN': Port('C_IN', 0),
                      'L': Port('L', 1),
                      'V': Port('V', 1)
                     }

    def Add_Corrente(self, corr, port):
        if port in self.Ports.keys():
            self.Ports[port].Corrente = corr
        self.Calc_Op()

    def Checa_Atributos(self):
        for p in self.Ports.values():
            if p.Corrente is None:
                return False
            elif not p.Corrente.Checa_Atributos() and p.Tipo == 0:
                return False
        return True

    def Calc_Op(self):
        if self.Checa_Atributos():
            c_in = self.Ports['C_IN'].Corrente
            L = self.Ports['L'].Corrente
            V = self.Ports['V'].Corrente
            V.z = c_in.y
            V.p = c_in.p
            V.T = c_in.T
            V.F = c_in.F*c_in.fvap
            L.z = c_in.x
            L.p = c_in.p
            L.T = c_in.T
            L.F = c_in.F*(1 - c_in.fvap)

    @property
    def Output(self):
        output = {'F_l': self.Ports['L'].Corrente.F} | {f'x_{(c.Nome)}': self.Ports['L'].Corrente.z[i].item()
                                                        for i,c in enumerate(self.FP.Componentes)
                                                       }
        output = output.copy() |{'F_v': self.Ports['V'].Corrente.F}|{f'y_{(c.Nome)}': self.Ports['V'].Corrente.z[i].item()
                                                                     for i,c in enumerate(self.FP.Componentes)
                                                                    }
        return pd.DataFrame(output, index = [self.Tag])

class ReatorGibbs(Operacao):

    __slots__ = Operacao.__slots__ + ['_ReatorGibbs__Tout']

    def __init__(self, nome, fp):
        super().__init__(nome, fp)
        self.__Tout = None
        self.Ports = {'C_IN': Port('C_IN', 0),
                      'L': Port('L', 1),
                      'V': Port('V', 1)
                     }

    @property
    def Tout(self):
        return self.__Tout
    @Tout.setter
    def Tout(self, valor):
        self.__Tout = valor
        self.Calc_Op()

    def Add_Corrente(self, corr, port):
        if port in self.Ports.keys():
            self.Ports[port].Corrente = corr
        self.Calc_Op()

    def Checa_Atributos(self):
        for p in self.Ports.values():
            if p.Corrente is None:
                return False
            elif not p.Corrente.Checa_Atributos() and p.Tipo == 0:
                return False
        if self.Tout is None:
            return False
        return True

    def Balanco_Atomico(self, componentes):
        a_df = pd.DataFrame({c.Nome: {k: v for k,v in c.Atomos.items()} for c in componentes}).fillna(0)
        A_atom = a_df.to_numpy()
        return A_atom

    def FObjetivo(self, n, corr):
        corr.F = np.sum(n)
        corr.z = n/corr.F
        f_obj = corr.F*corr.G/(config.R*corr.T)
        return f_obj

    def Calc_Op(self):
        if self.Checa_Atributos():
            c_out = self.Ports['C_IN'].Corrente.copy()
            z0, F0 = self.Ports['C_IN'].Corrente.z, self.Ports['C_IN'].Corrente.F
            n0 = z0*F0
            c_out.T = self.__Tout
            A_atom = self.Balanco_Atomico(c_out.FP.Componentes)
            limites = [(0, None)]*len(c_out.FP.Componentes)
            restricao = {'type': 'eq', 'fun' : lambda n : A_atom @ (n - n0)}
            resultado = minimize(self.FObjetivo, n0, args = (c_out), method = 'SLSQP', bounds = limites, constraints = restricao)
            L = self.Ports['L'].Corrente
            V = self.Ports['V'].Corrente
            V.z = c_out.y if np.sum(c_out.y) > 0 else c_out.x
            V.p = c_out.p
            V.T = c_out.T
            V.F = c_out.F*c_out.fvap
            L.z = c_out.x if np.sum(c_out.x) > 0 else c_out.y
            L.p = c_out.p
            L.T = c_out.T
            L.F = c_out.F*(1 - c_out.fvap)

    @property
    def Output(self):
        return None