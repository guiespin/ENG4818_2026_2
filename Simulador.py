from .utils import Chemical, Checa_Database, pd
from .Componente import Componente
from .PacoteTermodinamico import PacoteTermodinamico
from .Corrente import Corrente
from .Operacoes import Compressor, Mixer, ReatorGibbs, TrocadordeCalor, VasoFlash
#from .Convergencia import Reciclo
from . import config
import warnings
warnings.filterwarnings("ignore")

class Simulacao:
    __slots__ = ['_Simulacao__componentes',
                '_Simulacao__pacotetermodinamico',
                '_Simulacao__correntes',
                '_Simulacao__operacoes',
                '_Simulacao__sequencia'
               ]

    def __init__(self):
        self.__componentes = []
        self.__pacotetermodinamico = None
        self.__correntes = []
        self.__operacoes = []
        self.__sequencia = []

    @property
    def Componentes(self):
        return [c.Nome for c in self.__componentes]
    @Componentes.setter
    def Componentes(self, valor):
        pass

    @property
    def PacoteTermodinamico(self):
        return self.__pacotetermodinamico
    @PacoteTermodinamico.setter
    def PacoteTermodinamico(self, fp):
        self.__pacotetermodinamico = fp

    @property
    def Correntes(self):
        self.Calc_Fluxograma()
        return {c.Tag: c for c in self.__correntes}
    @Correntes.setter
    def Correntes(self, valor):
        pass

    @property
    def Operacoes(self):
        self.Calc_Fluxograma()
        return {o.Tag: o for o in self.__operacoes}
    @Operacoes.setter
    def Operacoes(self, valor):
        pass

    @property
    def Sequencia(self):
        return [o.Tag for o in self.__sequencia]
    @Sequencia.setter
    def Sequencia(self, valor):
        pass

    def Add_Componentes(self, *nome):
        for n in nome:
            if Checa_Database(n) and n not in self.Componentes:
                c = Componente(n)
                self.__componentes.append(c)

    def Add_PacoteTermodinamico(self, fp):
        fp_list = ['Ideal', 'Peng-Robinson', 'SRK']
        self.__pacotetermodinamico = PacoteTermodinamico(fp if fp in fp_list else 'Ideal', self.__componentes)

    def Set_Binarios(self, kij):
        self.__pacotetermodinamico.kij = kij

    def Add_Corrente(self, nome):
        self.__correntes.append(Corrente(nome, self.__pacotetermodinamico))

    def Add_Operacao(self, nome, tipo):
        op_dict = {'Compressor': Compressor,
                   'TrocadordeCalor': TrocadordeCalor,
                   'VasoFlash': VasoFlash,
                   'ReatorGibbs': ReatorGibbs,
                   'Mixer': Mixer
                   }
        if tipo in op_dict.keys():
            self.__operacoes.append(op_dict[tipo](nome, self.__pacotetermodinamico))

    #def Add_Reciclo(self, nome):
    #    self.__operacoes.append(Reciclo(nome, self.__pacotetermodinamico))
    #    self.Calc_Fluxograma()

    def Acopla_Corrente(self, nome_corrente, nome_operacao, port):
        if nome_operacao in self.Operacoes.keys():
            op = self.Operacoes[nome_operacao]
        if nome_corrente not in self.Correntes.keys():
            self.Add_Corrente(nome_corrente)
        corr = self.Correntes[nome_corrente]
        op.Add_Corrente(corr, port)
        self.Calc_Fluxograma()

    def Calc_Fluxograma(self):
        for op in self.__operacoes:
            op.Calc_Op()

    def Print_Correntes(self):
        return pd.DataFrame({c.Tag: {'Temperatura (K)': c.T,
                                     'Pressão (Pa)': c.p,
                                     'Fração de vapor': c.fvap,
                                     'Fluxo molar (kmol/h)': c.F,
                                     'Fluxo mássico (kg/h)': c.Fmas
                                    } 
                             for c in self.__correntes
                            }
                           ).transpose()

    def Print_Composiçoes(self):
        return pd.DataFrame({c.Tag: {x.Nome: c.z[i] for i,x in enumerate(self.__componentes)} 
                             for c in self.__correntes
                            }
                           ).transpose()
            