from utils import Chemical
import config

class Componente:
    
    __slots__ = ['Nome', 'MM', 'Atomos', 'phase', 'Tb', 'Tm', 'Tc', 'Pc', 'omega', 'dHf', 'dSf', 'dGf', 'S0g', 'CpG', 'CpL', 'CpS', 'dHv', 'dSv', 'Psat']

    def __init__(self, nome):
        chem = Chemical(nome, T = config.Tref, P = config.Pref)
        self.Nome = nome
        self.MM = chem.MW # g/mol
        self.Atomos = chem.atoms
        self.phase = chem.phase
        self.Tc = chem.Tc # K
        self.Pc = chem.Pc # Pa
        self.Tb = chem.Tb
        self.Tm = chem.Tm
        self.omega = chem.omega
        self.dHf = chem.Hfm
        self.dSf = chem.Sfm
        self.dGf = chem.Gfm
        self.S0g = chem.S0m
        self.CpG = chem.HeatCapacityGas
        self.CpL = chem.HeatCapacityLiquid
        self.CpS = chem.HeatCapacitySolid
        self.dHv = chem.EnthalpyVaporization
        if self.phase == 'g':
            self.dSv = chem.S_int_Tb_to_T_ref_g + self.dHv(self.Tb)/self.Tb + self.CpL.T_dependent_property_integral_over_T(config.Tref, self.Tb)
        elif self.phase == 'l':
            self.dSv = chem.S_int_l_T_ref_l_to_Tb + self.dHv(self.Tb)/self.Tb + self.CpG.T_dependent_property_integral_over_T(self.Tb, config.Tref)
        else:
            self.dSv = 0
        self.Psat = chem.VaporPressure
        