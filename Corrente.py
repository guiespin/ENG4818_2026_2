from utils import np

class Corrente:
    __slots__ = ['Tag',
                 '_Corrente__FP',
                 '_Corrente__MM',
                 '_Corrente__p',
                 '_Corrente__T',
                 '_Corrente__z',
                 '_Corrente__x',
                 '_Corrente__y',
                 '_Corrente__fvap',
                 '_Corrente__Vm',
                 '_Corrente__F',
                 '_Corrente__Fmas',
                 '_Corrente__Cp',
                 '_Corrente__HF',
                 '_Corrente__S',
                 '_Corrente__G'
                ]

    def __init__(self, nome, fp):
        for atr in self.__slots__:
            self.__setattr__(atr, None)
        self.Tag = nome
        self.__FP = fp

    def copy(self):
        copia = Corrente(self.Tag + '_Copia', self.__FP)
        copia.z = self.z
        copia.F = self.F
        copia.T = self.T
        copia.p = self.p
        return copia

    @property
    def FP(self):
        return self.__FP
    @FP.setter
    def FP(self, valor):
        pass

    @property
    def F(self):
        return self.__F
    @F.setter
    def F(self, valor):
        if valor is not None:
            if self.__Fmas is None:
                self.__F = valor
                self.__Fmas = valor*self.__MM if self.__MM is not None else None
                self.Calc_Flash()
            elif self.__Fmas/self.__MM != valor:
                self.__F = valor
                self.__Fmas = valor*self.__MM if self.__MM is not None else None
                self.Calc_Flash()

    @property
    def Fmas(self):
        return self.__Fmas
    @Fmas.setter
    def Fmas(self, valor):
        if valor is not None:
            if self.__F is None:
                self.__Fmas = valor
                self.__F = valor/self.__MM if self.__MM is not None else None
                self.Calc_Flash()
            elif self.__F*self.__MM != valor:
                self.__Fmas = valor
                self.__F = valor/self.__MM if self.__MM is not None else None
                self.Calc_Flash()

    @property
    def T(self):
        return self.__T
    @T.setter
    def T(self, valor):
        if self.__T != valor:
            self.__T = valor
            self.Calc_Flash()
    @T.deleter
    def T(self):
        self.__T = None

    @property
    def p(self):
        return self.__p
    @p.setter
    def p(self, valor):
        if self.__p != valor:
            self.__p = valor
            self.Calc_Flash()
    @p.deleter
    def p(self):
        self.__p = None

    @property
    def fvap(self):
        return self.__fvap
    @fvap.setter
    def fvap(self, valor):
        if all([valor >=0, valor <= 1, self.__fvap != valor]):
          self.__fvap = valor
          self.Calc_Flash()
    @fvap.deleter
    def fvap(self):
        self.__fvap = None

    @property
    def z(self):
        return self.__z
    @z.setter
    def z(self, valor):
        if isinstance(valor, np.ndarray):
            self.__z = np.nan_to_num(valor[:len(self.FP.Componentes)]/np.sum(valor[:len(self.FP.Componentes)]), 0)
            self.__MM = np.sum(self.__z*np.array([c.MM for c in self.__FP.Componentes])).item()
            self.F = self.__F
            self.Fmas = self.__Fmas
            self.Calc_Flash()

    @property
    def MM(self):
        return self.__MM
    @MM.setter
    def MM(self, valor):
        pass

    @property
    def x(self):
        return self.__x
    @x.setter
    def x(self, valor):
        pass

    @property
    def y(self):
        return self.__y
    @y.setter
    def y(self, valor):
        pass

    @property
    def Vm(self):
        return self.__Vm
    @Vm.setter
    def Vm(self, valor):
        pass

    @property
    def Cp(self):
        return self.__Cp
    @Cp.setter
    def Cp(self, valor):
        pass

    @property
    def HF(self):
        return self.__HF
    @HF.setter
    def HF(self, valor):
        pass

    @property
    def S(self):
        return self.__S
    @S.setter
    def S(self, valor):
        pass

    @property
    def G(self):
        return self.__G
    @G.setter
    def G(self, valor):
        pass

    def Checa_Flash(self):
        if self.__z is not None and self.__F is not None:
            b_l = [self.__T, self.__p, self.__fvap]
            if sum([bool(x) for x in b_l]) >= 2:
                return True
        return False

    def Calc_Flash(self):
        if self.Checa_Flash():
            self.__x, self.__y, self.__T, self.__p, self.__fvap, self.__Vm, self.__HF, self.__S, self.__G = self.FP.Calc_Flash(z = self.z, p = self.p, T = self.T, fvap = self.fvap)
            self.__HF *= self.F
            self.__Cp = self.FP.Calc_Cp(self.__x, self.__y, self.__fvap, self.__T)
            
    def Checa_Atributos(self):
        for a in self.__slots__:
            if self.__getattribute__(a) is None:
                return False
        return True