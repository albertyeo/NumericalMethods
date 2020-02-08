#Greeks

'''
Calculation of Greeks (Delta, Gamma, Vega, Theta) of European call/put and
American call and put based on:
1. Cox, Ross, and Rubinstein Model (CRR)
2. Jarrow-Rudd Risk Neutral Model (JRRN)
3. Jarrow-Rudd Equal Probability Model (JREQ)
4. Tian's Model
'''
#------------------------------------------------------------------------------
from enum import Enum
import math
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------
class PayoffType(Enum):
    Call = 0
    Put = 1

class EuropeanOption():
    def __init__(self, expiry, strike, payoffType):
        self.expiry = expiry
        self.strike = strike
        self.payoffType = payoffType
    def payoff(self, S):
        if self.payoffType == PayoffType.Call:
            return max(S - self.strike, 0)
        elif self.payoffType == PayoffType.Put:
            return max(self.strike - S, 0)
        else:
            raise Exception("payoffType not supported: ", self.payoffType)
    def valueAtNode(self, t, S, continuation):
        return continuation

class AmericanOption():
    def __init__(self, expiry, strike, payoffType):
        self.expiry = expiry
        self.strike = strike
        self.payoffType = payoffType
    def payoff(self, S):
        if self.payoffType == PayoffType.Call:
            return max(S - self.strike, 0)
        elif self.payoffType == PayoffType.Put:
            return max(self.strike - S, 0)
        else:
            raise Exception("payoffType not supported: ", self.payoffType)
    def valueAtNode(self, t, S, continuation):
        return max(self.payoff(S), continuation)

def crrCalib(r, vol, t):
    b = math.exp(vol * vol * t + r * t) + math.exp(-r * t)
    u = (b + math.sqrt(b * b - 4)) / 2
    p = (math.exp(r * t) - (1 / u)) / (u - 1 / u)
    return (u, 1/u, p)

def jrrnCalib(r, vol, t):
    u = math.exp((r - vol * vol / 2) * t + vol * math.sqrt(t))
    d = math.exp((r - vol * vol / 2) * t - vol * math.sqrt(t))
    p = (math.exp(r * t) - d) / (u - d)
    return (u, d, p)

def jreqCalib(r, vol, t):
    u = math.exp((r - vol * vol / 2) * t + vol * math.sqrt(t))
    d = math.exp((r - vol * vol / 2) * t - vol * math.sqrt(t))
    return (u, d, 1/2)

def tianCalib(r, vol, t):
    v = math.exp(vol * vol * t)
    u = 0.5 * math.exp(r * t) * v * (v + 1 + math.sqrt(v*v + 2*v - 3))
    d = 0.5 * math.exp(r * t) * v * (v + 1 - math.sqrt(v*v + 2*v - 3))
    p = (math.exp(r * t) - d) / (u - d)
    return (u, d, p)

def binomialPricer(S, r, vol, trade, n, calib):
    t = trade.expiry / n
    (u, d, p) = calib(r, vol, t)
    # set up the last time slice, there are n+1 nodes at the last time slice
    vs = [trade.payoff(S * u ** (n - i) * d ** i) for i in range(n + 1)]
    # iterate backward
    for i in range(n - 1, -1, -1):
        # calculate the value of each node at time slide i, there are i nodes
        for j in range(i + 1):
            nodeS = S * u ** (i - j) * d ** j
            continuation = math.exp(-r * t) * (vs[j] * p + vs[j + 1] * (1 - p))
            vs[j] = trade.valueAtNode(t * i, nodeS, continuation)
    return vs[0]

class Greeks(Enum):
    Delta = 0
    Gamma = 1
    Vega = 2
    Theta = 3

class GreekType():
    def __init__(self, n, greeks, calib, trade):
        self.n = n
        self.greeks = greeks
        self.calib = calib
        self.trade = trade

def binomialGreeks(S, r, vol, T, strike, greekType):
    dS = 0.05*S
    dvol = 0.05*vol
    dt = 0.05*T
    VS = lambda S: binomialPricer(S, r, vol, greekType.trade, greekType.n, greekType.calib)
    Vvol = lambda vol: binomialPricer(S, r, vol, greekType.trade, greekType.n, greekType.calib)
    if greekType.greeks == Greeks.Delta:
        return (VS(S+dS) - VS(S-dS))/(2*dS)
    elif greekType.greeks == Greeks.Gamma:
        return (VS(S+dS) - 2*VS(S) + VS(S-dS))/(dS**2)
    elif greekType.greeks == Greeks.Vega:
        return (Vvol(vol+dvol) - Vvol(vol-dvol))/(2*dvol)
    elif greekType.greeks == Greeks.Theta:
        return VS(S)/dt
    else:
        raise Exception("greekType not supported: ", greekType.greeks)
#------------------------------------------------------------------------------
S, r, vol, T = 100, 0.03, 0.2, 1
n = 300
ks = range(50, 150)
#------------------------------------------------------------------------------
#European Option
#------------------------------------------------------------------------------
#Delta
#------------------------------------------------------------------------------
crrDeltaCall, crrDeltaPut = [], []
tianDeltaCall, tianDeltaPut = [], []
jrrnDeltaCall, jrrnDeltaPut = [], []
jreqDeltaCall, jreqDeltaPut = [], []

for k in ks:
    crrDeltaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, crrCalib, EuropeanOption(T, k, PayoffType.Call))))
    crrDeltaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, crrCalib, EuropeanOption(T, k, PayoffType.Put))))    
    tianDeltaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, tianCalib, EuropeanOption(T, k, PayoffType.Call))))
    tianDeltaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, tianCalib, EuropeanOption(T, k, PayoffType.Put))))
    jrrnDeltaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, jrrnCalib, EuropeanOption(T, k, PayoffType.Call))))
    jrrnDeltaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, jrrnCalib, EuropeanOption(T, k, PayoffType.Put))))
    jreqDeltaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, jreqCalib, EuropeanOption(T, k, PayoffType.Call))))
    jreqDeltaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, jreqCalib, EuropeanOption(T, k, PayoffType.Put))))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

fig.suptitle('Delta')

ax1.plot(ks, crrDeltaCall, 'b', label='European Call')
ax1.plot(ks, crrDeltaPut, 'r', label='European Put')
ax1.set_ylabel('CRR Delta')

ax2.plot(ks, tianDeltaCall, 'b', label='European Call')
ax2.plot(ks, tianDeltaPut, 'r', label='European Put')
ax2.set_ylabel('Tian Delta')

ax3.plot(ks, jrrnDeltaCall, 'b', label='European Call')
ax3.plot(ks, jrrnDeltaPut, 'r', label='European Put')
ax3.set_ylabel('JRRN Delta')

ax4.plot(ks, jreqDeltaCall, 'b', label='European Call')
ax4.plot(ks, jreqDeltaPut, 'r', label='European Put')
ax4.set_ylabel('JREQ Delta')

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

plt.tight_layout()
plt.show()
#------------------------------------------------------------------------------
#Gamma
#------------------------------------------------------------------------------
crrGammaCall, crrGammaPut = [], []
tianGammaCall, tianGammaPut = [], []
jrrnGammaCall, jrrnGammaPut = [], []
jreqGammaCall, jreqGammaPut = [], []

for k in ks:
    crrGammaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, crrCalib, EuropeanOption(T, k, PayoffType.Call))))
    crrGammaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, crrCalib, EuropeanOption(T, k, PayoffType.Put))))    
    tianGammaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, tianCalib, EuropeanOption(T, k, PayoffType.Call))))
    tianGammaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, tianCalib, EuropeanOption(T, k, PayoffType.Put))))
    jrrnGammaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, jrrnCalib, EuropeanOption(T, k, PayoffType.Call))))
    jrrnGammaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, jrrnCalib, EuropeanOption(T, k, PayoffType.Put))))
    jreqGammaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, jreqCalib, EuropeanOption(T, k, PayoffType.Call))))
    jreqGammaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, jreqCalib, EuropeanOption(T, k, PayoffType.Put))))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

fig.suptitle('Gamma')

ax1.plot(ks, crrGammaCall, 'b', label='European Call')
ax1.plot(ks, crrGammaPut, 'r', label='European Put')
ax1.set_ylabel('CRR Gamma')

ax2.plot(ks, tianGammaCall, 'b', label='European Call')
ax2.plot(ks, tianGammaPut, 'r', label='European Put')
ax2.set_ylabel('Tian Gamma')

ax3.plot(ks, jrrnGammaCall, 'b', label='European Call')
ax3.plot(ks, jrrnGammaPut, 'r', label='European Put')
ax3.set_ylabel('JRRN Gamma')

ax4.plot(ks, jreqGammaCall, 'b', label='European Call')
ax4.plot(ks, jreqGammaPut, 'r', label='European Put')
ax4.set_ylabel('JREQ Gamma')

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

plt.tight_layout()
plt.show()
#------------------------------------------------------------------------------
#Vega
#------------------------------------------------------------------------------
crrVegaCall, crrVegaPut = [], []
tianVegaCall, tianVegaPut = [], []
jrrnVegaCall, jrrnVegaPut = [], []
jreqVegaCall, jreqVegaPut = [], []

for k in ks:
    crrVegaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, crrCalib, EuropeanOption(T, k, PayoffType.Call))))
    crrVegaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, crrCalib, EuropeanOption(T, k, PayoffType.Put))))    
    tianVegaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, tianCalib, EuropeanOption(T, k, PayoffType.Call))))
    tianVegaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, tianCalib, EuropeanOption(T, k, PayoffType.Put))))
    jrrnVegaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, jrrnCalib, EuropeanOption(T, k, PayoffType.Call))))
    jrrnVegaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, jrrnCalib, EuropeanOption(T, k, PayoffType.Put))))
    jreqVegaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, jreqCalib, EuropeanOption(T, k, PayoffType.Call))))
    jreqVegaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, jreqCalib, EuropeanOption(T, k, PayoffType.Put))))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

fig.suptitle('Vega')

ax1.plot(ks, crrVegaCall, 'b', label='European Call')
ax1.plot(ks, crrVegaPut, 'r', label='European Put')
ax1.set_ylabel('CRR Vega')

ax2.plot(ks, tianVegaCall, 'b', label='European Call')
ax2.plot(ks, tianVegaPut, 'r', label='European Put')
ax2.set_ylabel('Tian Vega')

ax3.plot(ks, jrrnVegaCall, 'b', label='European Call')
ax3.plot(ks, jrrnVegaPut, 'r', label='European Put')
ax3.set_ylabel('JRRN Vega')

ax4.plot(ks, jreqVegaCall, 'b', label='European Call')
ax4.plot(ks, jreqVegaPut, 'r', label='European Put')
ax4.set_ylabel('JREQ Vega')

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

plt.tight_layout()
plt.show()
#------------------------------------------------------------------------------
#Theta
#------------------------------------------------------------------------------
crrThetaCall, crrThetaPut = [], []
tianThetaCall, tianThetaPut = [], []
jrrnThetaCall, jrrnThetaPut = [], []
jreqThetaCall, jreqThetaPut = [], []

for k in ks:
    crrThetaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, EuropeanOption(
            0.95*T, k, PayoffType.Call)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, EuropeanOption(
            T, k, PayoffType.Call)))/(0.05*T))

    crrThetaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, EuropeanOption(
            0.95*T, k, PayoffType.Put)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, EuropeanOption(
            T, k, PayoffType.Put)))/(0.05*T))  

    tianThetaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, tianCalib, EuropeanOption(
            0.95*T, k, PayoffType.Call)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, tianCalib, EuropeanOption(
            T, k, PayoffType.Call)))/(0.05*T))

    tianThetaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, EuropeanOption(
            0.95*T, k, PayoffType.Put)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, EuropeanOption(
            T, k, PayoffType.Put)))/(0.05*T))  

    jrrnThetaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jrrnCalib, EuropeanOption(
            0.95*T, k, PayoffType.Call)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jrrnCalib, EuropeanOption(
            T, k, PayoffType.Call)))/(0.05*T))

    jrrnThetaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jrrnCalib, EuropeanOption(
            0.95*T, k, PayoffType.Put)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jrrnCalib, EuropeanOption(
            T, k, PayoffType.Put)))/(0.05*T)) 
        
    jreqThetaCall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jreqCalib, EuropeanOption(
            0.95*T, k, PayoffType.Call)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jreqCalib, EuropeanOption(
            T, k, PayoffType.Call)))/(0.05*T))

    jreqThetaPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jreqCalib, EuropeanOption(
            0.95*T, k, PayoffType.Put)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jreqCalib, EuropeanOption(
            T, k, PayoffType.Put)))/(0.05*T)) 

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

fig.suptitle('Theta')

ax1.plot(ks, crrThetaCall, 'b', label='European Call')
ax1.plot(ks, crrThetaPut, 'r', label='European Put')
ax1.set_ylabel('CRR Theta')

ax2.plot(ks, tianThetaCall, 'b', label='European Call')
ax2.plot(ks, tianThetaPut, 'r', label='European Put')
ax2.set_ylabel('Tian Theta')

ax3.plot(ks, jrrnThetaCall, 'b', label='European Call')
ax3.plot(ks, jrrnThetaPut, 'r', label='European Put')
ax3.set_ylabel('JRRN Theta')

ax4.plot(ks, jreqThetaCall, 'b', label='European Call')
ax4.plot(ks, jreqThetaPut, 'r', label='European Put')
ax4.set_ylabel('JREQ Theta')

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

plt.tight_layout()
plt.show()
#------------------------------------------------------------------------------
#American Option
#------------------------------------------------------------------------------
#Delta
#------------------------------------------------------------------------------
crrDeltaACall, crrDeltaAPut = [], []
tianDeltaACall, tianDeltaAPut = [], []
jrrnDeltaACall, jrrnDeltaAPut = [], []
jreqDeltaACall, jreqDeltaAPut = [], []

for k in ks:
    crrDeltaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, crrCalib, AmericanOption(T, k, PayoffType.Call))))
    crrDeltaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, crrCalib, AmericanOption(T, k, PayoffType.Put))))    
    tianDeltaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, tianCalib, AmericanOption(T, k, PayoffType.Call))))
    tianDeltaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, tianCalib, AmericanOption(T, k, PayoffType.Put))))
    jrrnDeltaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, jrrnCalib, AmericanOption(T, k, PayoffType.Call))))
    jrrnDeltaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, jrrnCalib, AmericanOption(T, k, PayoffType.Put))))
    jreqDeltaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, jreqCalib, AmericanOption(T, k, PayoffType.Call))))
    jreqDeltaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Delta, jreqCalib, AmericanOption(T, k, PayoffType.Put))))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

fig.suptitle('Delta')

ax1.plot(ks, crrDeltaACall, 'b', label='American Call')
ax1.plot(ks, crrDeltaAPut, 'r', label='American Put')
ax1.set_ylabel('CRR Delta')

ax2.plot(ks, tianDeltaACall, 'b', label='American Call')
ax2.plot(ks, tianDeltaAPut, 'r', label='American Put')
ax2.set_ylabel('Tian Delta')

ax3.plot(ks, jrrnDeltaACall, 'b', label='American Call')
ax3.plot(ks, jrrnDeltaAPut, 'r', label='American Put')
ax3.set_ylabel('JRRN Delta')

ax4.plot(ks, jreqDeltaACall, 'b', label='American Call')
ax4.plot(ks, jreqDeltaAPut, 'r', label='American Put')
ax4.set_ylabel('JREQ Delta')

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

plt.tight_layout()
plt.show()
#------------------------------------------------------------------------------
#Gamma
#------------------------------------------------------------------------------
crrGammaACall, crrGammaAPut = [], []
tianGammaACall, tianGammaAPut = [], []
jrrnGammaACall, jrrnGammaAPut = [], []
jreqGammaACall, jreqGammaAPut = [], []

for k in ks:
    crrGammaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, crrCalib, AmericanOption(T, k, PayoffType.Call))))
    crrGammaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, crrCalib, AmericanOption(T, k, PayoffType.Put))))    
    tianGammaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, tianCalib, AmericanOption(T, k, PayoffType.Call))))
    tianGammaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, tianCalib, AmericanOption(T, k, PayoffType.Put))))
    jrrnGammaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, jrrnCalib, AmericanOption(T, k, PayoffType.Call))))
    jrrnGammaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, jrrnCalib, AmericanOption(T, k, PayoffType.Put))))
    jreqGammaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, jreqCalib, AmericanOption(T, k, PayoffType.Call))))
    jreqGammaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Gamma, jreqCalib, AmericanOption(T, k, PayoffType.Put))))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

fig.suptitle('Gamma')

ax1.plot(ks, crrGammaACall, 'b', label='American Call')
ax1.plot(ks, crrGammaAPut, 'r', label='American Put')
ax1.set_ylabel('CRR Gamma')

ax2.plot(ks, tianGammaACall, 'b', label='American Call')
ax2.plot(ks, tianGammaAPut, 'r', label='American Put')
ax2.set_ylabel('Tian Gamma')

ax3.plot(ks, jrrnGammaACall, 'b', label='American Call')
ax3.plot(ks, jrrnGammaAPut, 'r', label='American Put')
ax3.set_ylabel('JRRN Gamma')

ax4.plot(ks, jreqGammaACall, 'b', label='American Call')
ax4.plot(ks, jreqGammaAPut, 'r', label='American Put')
ax4.set_ylabel('JREQ Gamma')

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

plt.tight_layout()
plt.show()
#------------------------------------------------------------------------------
#Vega
#------------------------------------------------------------------------------
crrVegaACall, crrVegaAPut = [], []
tianVegaACall, tianVegaAPut = [], []
jrrnVegaACall, jrrnVegaAPut = [], []
jreqVegaACall, jreqVegaAPut = [], []

for k in ks:
    crrVegaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, crrCalib, AmericanOption(T, k, PayoffType.Call))))
    crrVegaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, crrCalib, AmericanOption(T, k, PayoffType.Put))))    
    tianVegaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, tianCalib, AmericanOption(T, k, PayoffType.Call))))
    tianVegaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, tianCalib, AmericanOption(T, k, PayoffType.Put))))
    jrrnVegaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, jrrnCalib, AmericanOption(T, k, PayoffType.Call))))
    jrrnVegaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, jrrnCalib, AmericanOption(T, k, PayoffType.Put))))
    jreqVegaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, jreqCalib, AmericanOption(T, k, PayoffType.Call))))
    jreqVegaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Vega, jreqCalib, AmericanOption(T, k, PayoffType.Put))))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

fig.suptitle('Vega')

ax1.plot(ks, crrVegaACall, 'b', label='American Call')
ax1.plot(ks, crrVegaAPut, 'r', label='American Put')
ax1.set_ylabel('CRR Vega')

ax2.plot(ks, tianVegaACall, 'b', label='American Call')
ax2.plot(ks, tianVegaAPut, 'r', label='American Put')
ax2.set_ylabel('Tian Vega')

ax3.plot(ks, jrrnVegaACall, 'b', label='American Call')
ax3.plot(ks, jrrnVegaAPut, 'r', label='American Put')
ax3.set_ylabel('JRRN Vega')

ax4.plot(ks, jreqVegaACall, 'b', label='American Call')
ax4.plot(ks, jreqVegaAPut, 'r', label='American Put')
ax4.set_ylabel('JREQ Vega')

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

plt.tight_layout()
plt.show()
#------------------------------------------------------------------------------
#Theta
#------------------------------------------------------------------------------
crrThetaACall, crrThetaAPut = [], []
tianThetaACall, tianThetaAPut = [], []
jrrnThetaACall, jrrnThetaAPut = [], []
jreqThetaACall, jreqThetaAPut = [], []

for k in ks:
    crrThetaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, AmericanOption(
            0.95*T, k, PayoffType.Call)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, AmericanOption(
            T, k, PayoffType.Call)))/(0.05*T))

    crrThetaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, AmericanOption(
            0.95*T, k, PayoffType.Put)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, AmericanOption(
            T, k, PayoffType.Put)))/(0.05*T))  

    tianThetaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, tianCalib, AmericanOption(
            0.95*T, k, PayoffType.Call)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, tianCalib, AmericanOption(
            T, k, PayoffType.Call)))/(0.05*T))

    tianThetaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, AmericanOption(
            0.95*T, k, PayoffType.Put)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, crrCalib, AmericanOption(
            T, k, PayoffType.Put)))/(0.05*T))  

    jrrnThetaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jrrnCalib, AmericanOption(
            0.95*T, k, PayoffType.Call)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jrrnCalib, AmericanOption(
            T, k, PayoffType.Call)))/(0.05*T))

    jrrnThetaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jrrnCalib, AmericanOption(
            0.95*T, k, PayoffType.Put)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jrrnCalib, AmericanOption(
            T, k, PayoffType.Put)))/(0.05*T)) 
        
    jreqThetaACall.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jreqCalib, AmericanOption(
            0.95*T, k, PayoffType.Call)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jreqCalib, AmericanOption(
            T, k, PayoffType.Call)))/(0.05*T))

    jreqThetaAPut.append(binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jreqCalib, AmericanOption(
            0.95*T, k, PayoffType.Put)))/(0.05*T) -
            binomialGreeks(S, r, vol, T, k, GreekType(
            n, Greeks.Theta, jreqCalib, AmericanOption(
            T, k, PayoffType.Put)))/(0.05*T)) 

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

fig.suptitle('Theta')

ax1.plot(ks, crrThetaACall, 'b', label='American Call')
ax1.plot(ks, crrThetaAPut, 'r', label='American Put')
ax1.set_ylabel('CRR Theta')

ax2.plot(ks, tianThetaACall, 'b', label='American Call')
ax2.plot(ks, tianThetaAPut, 'r', label='American Put')
ax2.set_ylabel('Tian Theta')

ax3.plot(ks, jrrnThetaACall, 'b', label='American Call')
ax3.plot(ks, jrrnThetaAPut, 'r', label='American Put')
ax3.set_ylabel('JRRN Theta')

ax4.plot(ks, jreqThetaACall, 'b', label='American Call')
ax4.plot(ks, jreqThetaAPut, 'r', label='American Put')
ax4.set_ylabel('JREQ Theta')

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

plt.tight_layout()
plt.show()
