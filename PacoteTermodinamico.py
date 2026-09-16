from .utils import Flash_Loop_TP, Flash_Loop_Tf, Flash_Loop_pf, Hres_Ideal, Hres_PengRobinson, Hres_SRK, EOS_Ideal, EOS_PengRobinson, EOS_SRK, np
from . import config

class PacoteTermodinamico:

    __slots__ = ['Nome', 'Componentes', 'Tc', 'Pc', 'w', 'kij']

    def __init__(self, Nome, Comp_list):
        self.Nome = Nome.replace('-', '')
        self.Componentes = Comp_list
        self.Tc = np.array([c.Tc for c in self.Componentes])
        self.Pc = np.array([c.Pc for c in self.Componentes])
        self.w = np.array([c.omega for c in self.Componentes])
        self.kij = np.zeros_like((len(self.Componentes), len(self.Componentes)))

    def Calc_Cp(self, x, y, fvap, T):
        if T is not None and fvap is not None:
            CpL = np.sum(x*np.array([c.CpL(T) for c in self.Componentes]))
            CpG = np.sum(y*np.array([c.CpG(T) for c in self.Componentes]))
            return (CpL*(1 - fvap) + CpG*fvap).item()
        else:
            return None
        
    def Flash_TP(self, T, p, z):
        fvap, x, y = Flash_Loop_TP(T, p, z, self.Tc, self.Pc, self.w, self.kij, self.Nome)
        return x, y, T, p, fvap

    def Flash_Tf(self, z, T, fvap):
        p, x, y = Flash_Loop_Tf(T, fvap, z, self.Tc, self.Pc, self.w, self.kij, self.Nome)
        return x, y, T, p, fvap

    def Flash_pf(self, z, p, fvap):
        T, x, y = Flash_Loop_pf(p, fvap, z, self.Tc, self.Pc, self.w, self.kij, self.Nome)
        return x, y, T, p, fvap

    '''
    Método para cálculo de flash de uma corrente, dividido em três métodos para
    flash PT (cálculo da fração de vapor), flash Pf (cálculo da temperatura) e
    flash Tf (cálculo da pressão).
    '''

    def Calc_Flash(self, z, p = None, T = None, fvap = None):
        R = config.R
        x = np.zeros_like(z)
        y = x.copy()
        if all([p is not None, T is not None]):
            x, y, T, p, fvap = self.Flash_TP(T, p, z)
    
        elif all([T is not None, fvap is not None]):
            x, y, T, p, fvap = self.Flash_Tf(z, T, fvap)
    
        elif all([p is not None, fvap is not None]):
            x, y, T, p, fvap = self.Flash_pf(z, p, fvap)
            
        Z = globals()[f'EOS_{self.Nome}'](p, T, z, self.Tc, self.Pc, self.w, self.kij) if all([p is not None, T is not None]) else 1
        try:
            Z = Z[9]
        except:
            pass
        try:
            Vm = np.max(Z)*R*T/p
        except:
            Vm = 0

        HF, S = self.Calc_H_S(z, T, p, x, y, fvap)
        G = HF - T*S
        return x, y, T, p, fvap, Vm, HF, S, G

    def Calc_H_S(self, z, T, p, x, y, fvap):
        if self.Nome == 'Ideal':
            dHf_l = np.array([c.dHf if c.phase == 'l' else c.dHf - c.dHv(config.Tref) for c in self.Componentes])
            H_id_L = np.sum(x*(np.array([c.CpL.T_dependent_property_integral(config.Tref, T) for c in self.Componentes]) + dHf_l))
            dHf_v = np.array([c.dHf if c.phase == 'g' else c.dHf + c.dHv(config.Tref) for c in self.Componentes])
            H_id_V = np.sum(y*(np.array([c.CpG.T_dependent_property_integral(config.Tref, T) for c in self.Componentes]) + dHf_v))
            dSf_l = np.array([c.dSf if c.phase == 'l' else c.dSf - c.dSv for c in self.Componentes])
            S_id_L = np.sum(x*(np.array([c.CpL.T_dependent_property_integral_over_T(config.Tref, T)
                                         for c in self.Componentes]) + dSf_l)) - config.R*np.sum(x*np.log(x), where = x > 0) - config.R*np.log(p/config.Pref)
            dSf_v = np.array([c.dSf if c.phase == 'g' else c.dSf + c.dSv for c in self.Componentes])
            S_id_V = np.sum(y*(np.array([c.CpG.T_dependent_property_integral_over_T(config.Tref, T)
                                         for c in self.Componentes]) + dSf_v)) - config.R*np.sum(y*np.log(y), where = y > 0) - config.R*np.log(p/config.Pref)
            return (1 - fvap)*H_id_L + fvap*H_id_V, (1 - fvap)*S_id_L + fvap*S_id_V
        dHf_ig = np.array([c.dHf if c.phase == 'g' else c.dHf + c.dHv(config.Tref) for c in self.Componentes])
        H_ig = np.sum(z*(np.array([c.CpG.T_dependent_property_integral(config.Tref, T) for c in self.Componentes]) + dHf_ig))
        dSf = np.array([c.dSf if c.phase == 'g' else c.dSf + c.dSv for c in self.Componentes])
        S_id_L = np.sum(x*(np.array([c.CpG.T_dependent_property_integral_over_T(config.Tref, T)
                                     for c in self.Componentes]) + dSf)) - config.R*np.sum(x*np.log(x), where = x > 0) - config.R*np.log(p/config.Pref)
        S_id_V = np.sum(y*(np.array([c.CpG.T_dependent_property_integral_over_T(config.Tref, T)
                                     for c in self.Componentes]) + dSf)) - config.R*np.sum(y*np.log(y), where = y > 0) - config.R*np.log(p/config.Pref)
        H_res_L, S_res_L = globals()[f'Hres_{self.Nome}'](p, T, x, self.Tc, self.Pc, self.w, self.kij, fase = 'L')
        H_res_V, S_res_V = globals()[f'Hres_{self.Nome}'](p, T, y, self.Tc, self.Pc, self.w, self.kij, fase = 'V')
        if np.isclose(fvap, 1., rtol = 0, atol = 1e-4) and H_res_L != H_res_L:
            H_res_L = 0
            S_res_L = 0
        elif np.isclose(fvap, 0., rtol = 0, atol = 1e-4) and H_res_V != H_res_V:
            H_res_V = 0
            S_res_V = 0
        return (1 - fvap)*H_res_L + fvap*H_res_V + H_ig, (1 - fvap)*(S_res_L + S_id_L) + fvap*(S_res_V + S_id_V)        