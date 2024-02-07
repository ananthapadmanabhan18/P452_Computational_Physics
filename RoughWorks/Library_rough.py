import numpy as np
import math as m

class Numerical_integration:
    def __init__(self, f, a, b,n):
        self.a = a
        self.b = b
        self.n = n
        self.f = f

    def printf(self):
        print(self.a, self.b, self.n, self.f)

    def simpsons(self):
        h=(self.b-self.a)/self.n
        I=0
        x=self.a
        while x<=self.b:
            I+=h/3*(self.f(x)+4*self.f(x+h/2)+self.f(x+h))*0.5
            x+=h
        return I 