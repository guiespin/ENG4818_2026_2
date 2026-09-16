import pandas as pd
import numpy as np
from scipy.optimize import brentq, minimize
import thermo
from thermo import Chemical
from . import config

'''
Utilidades gerais
'''

def Checa_Database(nome):
    try:
        thermo.search_chemical(nome)
        return True
    except ValueError:
        print(f'Componente \033[1m{nome}\033[0m não encontrado no banco de dados!')
        return False

'''
Equações de Estado
'''

def EOS_Ideal(p, T, z, Tc, Pc, w, kij):
    Z = 1
    return Z

def EOS_PengRobinson(p, T, z, Tc, Pc, w, kij):
    R = config.R
    m = np.where(w <= 0.49, 0.37464 + 1.54226*w - 0.26992*(w**2), 0.379642 + 1.48503*w - 0.164423*(w**2) + 0.016666*(w**3))
    Tr = T/Tc
    alpha = (1 + m*(1 - Tr**0.5))**2
    a_i = 0.457235*((R*Tc)**2/Pc)*alpha
    b_i = 0.077796*R*Tc/Pc
    a = np.sum(z[:, None]*z[None,:]*((a_i[:,None]*a_i[None,:])**0.5)*(1 - kij))
    b = np.sum(z*b_i)
    A = a*p/((R*T)**2)
    B = b*p/(R*T)
    Z = np.roots([1, -(1 - B), A - 2*B - 3*(B**2), - (A*B - B**2 - B**3)])
    Z = np.real_if_close(Z, tol = 1e10)
    Z = Z[np.isreal(Z)].real
    Z = Z[Z > B]
    Vm = np.max(Z)*R*T/p
    return m, alpha, a_i, b_i, a, b, A, B, Z, Vm

def EOS_SRK(p, T, z, Tc, Pc, w, kij):
    R = config.R
    m = 0.480 + 1.574*w - 0.176*(w**2)
    Tr = T/Tc
    alpha = (1 + m*(1 - Tr**0.5))**2
    a_i = 0.42747*((R*Tc)**2/Pc)*alpha
    b_i = 0.08664*R*Tc/Pc
    a = np.sum(z[:, None]*z[None,:]*((a_i[:,None]*a_i[None,:])**0.5)*(1 - kij))
    b = np.sum(z*b_i)
    A = a*p/((R*T)**2)
    B = b*p/(R*T)
    Z = np.roots([1, -1, A-B-B**2, -A*B])
    Z = np.real_if_close(Z, tol = 1e10)
    Z = Z[np.isreal(Z)].real
    Z = Z[Z > B]
    Vm = np.max(Z)*R*T/p
    return m, alpha, a_i, b_i, a, b, A, B, Z, Vm

'''
Cálculos de fugacidade e equilíbrio líquido-vapor
'''

def K_Wilson(p, T, Tc, Pc, w):
    return (Pc/p)*np.exp(5.373*(1 + w)*(1 - Tc/T))

def Rachford_Rice(z, K):
    if np.sum(z * (K - 1)) <= 0:
        return 0.0, z.copy(), K*z.copy()
    if np.sum(z * (K - 1) / K) >= 0:
        return 1.0, z.copy()/K, z.copy()
    f = lambda x : np.sum(z*(K - 1)/(1 + x*(K - 1)))
    fvap = brentq(f, 0.0, 1.0)
    x = z/(1 + fvap*(K - 1))
    y = K*x
    return fvap, x, y

def Phi_Ideal(T, p, z, Tc, Pc, w, kij, fase = 'V'):
    return 0

def Phi_PengRobinson(T, p, z, Tc, Pc, w, kij, fase = 'V'):
    m, alpha, a_i, b_i, a, b, A, B, Z, Vm = EOS_PengRobinson(p, T, z, Tc, Pc, w, kij)
    if fase == 'L':
        Z = np.min(Z)
    else:
        Z = np.max(Z)
    ln_phi = (b_i/b)*(Z - 1) - np.log(Z - B) - (A/(np.sqrt(8)*B))*(2*np.sum(z*np.sqrt(a_i[:, None]*a_i[None, :])*(1 - kij), axis = 1)/a - b_i/b)*np.log((Z + (1 + np.sqrt(2))*B)/(Z + (1 - np.sqrt(2))*B))
    return ln_phi
    
def Phi_SRK(T, p, z, Tc, Pc, w, kij, fase = 'V'):
    m, alpha, a_i, b_i, a, b, A, B, Z, Vm = EOS_SRK(p, T, z, Tc, Pc, w, kij)
    if fase == 'L':
        Z = np.min(Z)
    else:
        Z = np.max(Z)
    ln_phi = (b_i/b)*(Z - 1) - np.log(Z - B) - (A/B)*(2*np.sum(z*np.sqrt(a_i[:, None]*a_i[None, :])*(1 - kij), axis = 1)/a - b_i/b)*np.log(1 + B/Z)
    return ln_phi

def K_iter(z, K, Phi_FP):
    fvap, x, y = Rachford_Rice(z, K)
    ln_phi_L = Phi_FP(T, p, x, Tc, Pc, w, kij, fase = 'L')
    ln_phi_V = Phi_FP(T, p, y, Tc, Pc, w, kij, fase = 'V')
    ln_K = ln_phi_L - ln_phi_V
    return np.exp(ln_K)

def Flash_Loop_TP(T, p, z, Tc, Pc, w, kij, FP, tol = 1e-6, maxit = 100):
    Phi_FP = globals()[f'Phi_{FP}']
    K = K_Wilson(p, T, Tc, Pc, w)
    if FP == 'Ideal':
        fvap, x, y = Rachford_Rice(z, K)
        return fvap, x, y
    for _ in range(maxit):
        fvap, x, y = Rachford_Rice(z, K)
        if np.isclose(fvap, 1., rtol = 0, atol = 1e-4):
            if _ > 0:
                return 1., z/K, z
        if np.isclose(fvap, 0., rtol = 0, atol = 1e-4):
            if _ > 0:
                return 0., z, z*K
        ln_phi_L = Phi_FP(T, p, x, Tc, Pc, w, kij, fase = 'L')
        ln_phi_V = Phi_FP(T, p, y, Tc, Pc, w, kij, fase = 'V')
        ln_K = ln_phi_L - ln_phi_V
        if np.all(np.abs(ln_K - np.log(K)) < tol):
            return fvap, x, y
        K = np.exp(ln_K)
    print(f'Flash não convergiu em {maxit} iterações')
    return fvap, x, y

def Flash_Loop_Tf(T, fvap, z, Tc, Pc, w, kij, FP, tol = 1e-5, maxit = 50):
    Phi_FP = globals()[f'Phi_{FP}']
    if fvap == 0:
        x = z.copy()
        def f(p):
            K = [K_Wilson(p, T, Tc, Pc, w)]
            for _ in range(maxit):
                y = K[-1]*x
                ln_phi_L = Phi_FP(T, p, x, Tc, Pc, w, kij, fase = 'L')
                ln_phi_V = Phi_FP(T, p, y, Tc, Pc, w, kij, fase = 'V')
                ln_K = ln_phi_L - ln_phi_V
                if np.all(np.abs(ln_K - np.log(K[-1])) < tol):
                    return np.sum(K[-1]*x) - 1, K[-1]
                K.append(np.exp(ln_K))
            return np.sum(K[-1]*x) - 1, K[-1]
        g = lambda p: f(p)[0]
        p = brentq(g, 1e-5, 1e10)
        K = f(p)[1]
        return p, x, K*x
            
    elif fvap == 1:
        y = z.copy()
        def f(p):
            K = [K_Wilson(p, T, Tc, Pc, w)]
            for _ in range(maxit):
                x = y/K[-1]
                ln_phi_L = Phi_FP(T, p, x, Tc, Pc, w, kij, fase = 'L')
                ln_phi_V = Phi_FP(T, p, y, Tc, Pc, w, kij, fase = 'V')
                ln_K = ln_phi_L - ln_phi_V
                if np.all(np.abs(ln_K - np.log(K[-1])) < tol):
                    return np.sum(y/K[-1]) - 1, K[-1]
                K.append(np.exp(ln_K))
            return np.sum(y/K[-1]) - 1, K[-1]
        g = lambda p: f(p)[0]
        p = p = brentq(g, 1e-5, 1e10)
        K = f(p)[1]
        return p, y/K, y
        
    else:
        def f(p):
            K = [K_Wilson(p, T, Tc, Pc, w)]
            for _ in range(maxit):
                x = z/(1 + fvap*(K[-1] - 1))
                y = K[-1]*x
                ln_phi_L = Phi_FP(T, p, x, Tc, Pc, w, kij, fase = 'L')
                ln_phi_V = Phi_FP(T, p, y, Tc, Pc, w, kij, fase = 'V')
                ln_K = ln_phi_L - ln_phi_V
                if np.all(np.abs(ln_K - np.log(K[-1])) < tol):
                    return np.sum(z*(K[-1] - 1)/(1 + fvap*(K[-1] - 1))), K[-1]
                K.append(np.exp(ln_K))
            return np.sum(x) - 1, K[-1]
        g = lambda p : f(p)[0]
        p = p = brentq(g, 1e-5, 1e10)
        K = f(p)[1]
        x = z/(1 + fvap*(K[-1] - 1))
        y = K*x
        return p, x, y

def Flash_Loop_pf(p, fvap, z, Tc, Pc, w, kij, FP, tol = 1e-5, maxit = 50):
    Phi_FP = globals()[f'Phi_{FP}']
    if fvap == 0:
        x = z.copy()
        def f(T):
            K = [K_Wilson(p, T, Tc, Pc, w)]
            for _ in range(maxit):
                y = K[-1]*x
                ln_phi_L = Phi_FP(T, p, x, Tc, Pc, w, kij, fase = 'L')
                ln_phi_V = Phi_FP(T, p, y, Tc, Pc, w, kij, fase = 'V')
                ln_K = ln_phi_L - ln_phi_V
                if np.all(np.abs(ln_K - np.log(K[-1])) < tol):
                    return np.sum(K[-1]*x) - 1, K[-1]
                K.append(np.exp(ln_K))
            return np.sum(K[-1]*x) - 1, K[-1]
        g = lambda T: f(T)[0]
        T = brentq(g, 1e-5, 3_000)
        K = f(T)[1]
        return T, x, K*x
            
    elif fvap == 1:
        y = z.copy()
        def f(T):
            K = [K_Wilson(p, T, Tc, Pc, w)]
            for _ in range(maxit):
                x = y/K[-1]
                ln_phi_L = Phi_FP(T, p, x, Tc, Pc, w, kij, fase = 'L')
                ln_phi_V = Phi_FP(T, p, y, Tc, Pc, w, kij, fase = 'V')
                ln_K = ln_phi_L - ln_phi_V
                if np.all(np.abs(ln_K - np.log(K[-1])) < tol):
                    return np.sum(y/K[-1]) - 1, K[-1]
                K.append(np.exp(ln_K))
            return np.sum(y/K[-1]) - 1, K[-1]
        g = lambda T: f(T)[0]
        T = brentq(g, 1e-5, 3_000)
        K = f(T)[1]
        return T, y/K, y
        
    else:
        def f(T):
            K = [K_Wilson(p, T, Tc, Pc, w)]
            for _ in range(maxit):
                x = z/(1 + fvap*(K[-1] - 1))
                y = K[-1]*x
                ln_phi_L = Phi_FP(T, p, x, Tc, Pc, w, kij, fase = 'L')
                ln_phi_V = Phi_FP(T, p, y, Tc, Pc, w, kij, fase = 'V')
                ln_K = ln_phi_L - ln_phi_V
                if np.all(np.abs(ln_K - np.log(K[-1])) < tol):
                    return np.sum(z*(K[-1] - 1)/(1 + fvap*(K[-1] - 1))), K[-1]
                K.append(np.exp(ln_K))
            return np.sum(x) - 1, K[-1]
        g = lambda T : f(T)[0]
        T = brentq(g, 1e-5, 3_000)
        K = f(T)[1]
        x = z/(1 + fvap*(K[-1] - 1))
        y = K*x
        return T, x, y
        
'''
Cálculos de entalpia e entropia residuais
'''

def Hres_Ideal(p, T, z, c, Tc, Pc, w, kij, fase = 'V'):
    return 0, 0

def Hres_PengRobinson(p, T, c, Tc, Pc, w, kij, fase = 'V'):
    R = config.R
    m, alpha, a_i, b_i, a, b, A, B, Z, Vm = EOS_PengRobinson(p, T, c, Tc, Pc, w, kij)
    if fase == 'L':
        Z = np.min(Z)
    else:
        Z = np.max(Z)
    Tr = T/Tc
    dai_dT = -a_i*m/(Tc*np.sqrt(Tr*alpha))
    da_dT = np.sum(c[:, None]*c[None, :]*np.sqrt(a_i[:, None]*a_i[None, :])*(1 - kij)*(1/2)*(dai_dT[:, None]/a_i[:, None] + dai_dT[None, :]/a_i[None, :]))
    H_res = R*T*(Z - 1) + ((T*da_dT - a)/(b*np.sqrt(8)))*np.log((Z + (1 + np.sqrt(2))*B)/(Z + (1 - np.sqrt(2))*B))
    S_res = R*np.log(Z - B) - np.log(p/config.Pref) - (1/(np.sqrt(8)*b))*da_dT*np.log((Z + (1 + np.sqrt(2))*B)/(Z + (1 - np.sqrt(2))*B))
    return H_res, S_res

def Hres_SRK(p, T, c, Tc, Pc, w, kij, fase = 'V'):
    R = config.R
    m, alpha, a_i, b_i, a, b, A, B, Z, Vm = EOS_SRK(p, T, c, Tc, Pc, w, kij)
    if fase == 'L':
        Z = np.min(Z)
    else:
        Z = np.max(Z)
    Tr = T/Tc
    dai_dT = -a_i*m/(Tc*np.sqrt(Tr*alpha))
    da_dT = np.sum(c[:, None]*c[None, :]*np.sqrt(a_i[:, None]*a_i[None, :])*(1 - kij)*(1/2)*(dai_dT[:, None]/a_i[:, None] + dai_dT[None, :]/a_i[None, :]))
    H_res = R*T*(Z - 1) + ((T*da_dT - a)/b)*np.log(1 + B/Z)
    S_res = R*np.log(Z - B) - np.log(p/config.Pref) + da_dT*(1/b)*np.log(1 + B/Z)
    return H_res, S_res