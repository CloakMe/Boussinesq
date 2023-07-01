import numpy as np
from mpmath import *

BoussinesqSoliton(x,t,c,alpha,beta1,beta2)

def BoussinesqSoliton(x,t,c,alpha,beta1,beta2):
    beta=beta1/beta2
    sech( sqrt( beta1*(c^2-1)/(beta1*c^2 - beta2) )*(x-c*t*sqrt(beta))/2 )
    u=3*(c^2 - 1)* /(2*alpha)
