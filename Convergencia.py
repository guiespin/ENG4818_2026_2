from .utils import np, pd
from .Operacoes import Operacao

class Reciclo(Operacao):

    def __init__(self, nome, fp):
        super().__init__(nome, fp)