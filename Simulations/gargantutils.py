import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

def airDensity(temperature=295,staticPressure=1e5):
    """
    Donne la densité de l'air pour une température et une pression statique donnée (Beranek et Mellow 2012 eq. 1.6)
    Paramètres:
        temperature: température de l'air en K
        staticPressure: pression statique en Pa
    Rend la densité de l'air en kg/m3
    Devrait retourner 1.18kg/m3 pour T=295K et P0=1e5
    """
    return staticPressure/(287*temperature)
def helmholtzParameters(S, V, L, Q=5):
    """
    Décrit un résonnateur d'Helmholtz
    Retourne la constante de rappel (s), la masse (m) et la résistance (R) et meme la frequence de resonnance (w0)
    pour une aire et longueure du tuyau (S,L), un volume de compression (V) et un facteur de qualité donné
    """
    m=airDensity()*S*L
    s= airDensity()*airSoundSpeed()**2*S**2/V
    w0=airSoundSpeed()*np.sqrt(S/(L*V))
    R=w0*m/Q
    return s, m, R,w0

def airSoundSpeed(temperature=295):
    """
    Donne la vitesse du son dans l'air pour une température donnée (Beranek et Mellow 2012 eq. 1.8)
    Paramètres:
        temperature: température de l'air en K
    Rend la vitesse du son en m/s (344.8 m/s pour T=295K)
    """
    return 331.5*np.sqrt(temperature/273)
def freq2k(freq,temperature=295):
    """
    Convertit la fréquence en Hz en nombre d'onde dans l'air
    freq: la fréquence de l'onde en Hz
    temperature: la température de l'air
    retourne le nombre d'onde en rad/m
    """
    return 2*np.pi*freq/airSoundSpeed(temperature=temperature)
def airCharacImpedance(temperature=295):
    """
    Impédance charactéristique de l'air à une température donnée
        temperature: température de l'air en K
    Rend l'impédance charactéristique de l'air en N s/m^3 (rayls) (407 rayls pour T=295K)
    """
    return airDensity(temperature=temperature)*airSoundSpeed(temperature=temperature)
def freq2lambda(freq,temperature=295):
    """
    Convertit la fréquence en Hz en longueure d'onde dans l'air
    freq: la fréquence de l'onde en Hz
    temperature: la température de l'air
    retourne la longueure d'onde en m
    """
    return airSoundSpeed(temperature=temperature)/freq

def radiativeImpedance(radius,wavevector,temperature=295):
    """
    Donne l'impédance méchanique radiative d'un piston isolé (flanged) eq. 7.5.13 et 7.5.14 KinsleyFreyCoopSanders
    radius: rayon du piston en m
    wavector: nombre d'onde en rad/m 
    retourne l'impédance radiative complexe en rayls*m2
    """
    area=np.pi*radius**2
    Z0=airCharacImpedance(temperature=temperature)
    Zr=Z0*area*(0.5*(radius*wavevector)**2)+1j*8/(3*np.pi)*wavevector*radius
    return Zr
def tubeImpedance(length,wavevector,radius,temperature=295):
    """
    Donne l'impédance mécanique regardant vers l'entrée d'un tube 1D chargé à l'autre extrémité (rayls)
    length: longueure du tube en m
    wavevector: nombre d'onde en rad/m
    radius: rayon du tube en m
    retourne impédance du tube chargé vue par la bouche en rayls*m^2
    """ 
    area=np.pi*radius**2
    Zm0=airCharacImpedance(temperature=temperature)*area
    loadImpedance=radiativeImpedance(radius,wavevector)
    inputImpedance=(loadImpedance+1j*Zm0*area*np.tan(wavevector*length))/(1+1j*(loadImpedance/(Zm0)*np.tan(wavevector*length)))
    return inputImpedance

def driveImpedance(v,R,m,s):
    """
        Rends l'impédance mécanique d'un oscillateur forcé vue de la force (eq. 10.6.3 KinsleyFreyCoopSanders)
        R: résistance mécanique de l'oscillateur (rayl.m2)
        v: fréquence de forcage (Hz)
        m: masse de l'oscillateur
        s: constante de rappel de l'oscillateur
    """    
    return R+1j*(2*np.pi*v*m-s/(2*np.pi*v))
def drivenTubeResonnances(length,radius,R,m,s,maxHarm=10,temperature=295):
    """
        Calcule les fréquences de résonnance d'un piston oscillateur monté sur un tube chargé et ouvert selon eq. 10.6.3 KinsleyFreyCoopSanders
        length: longueure du tube en m
        loadImpedance: impédance de la charge en rayls
        radius: rayon du tube en m
        R: résistance mécanique de l'oscillateur (rayl.m2)
        m: masse de l'oscillateur
        s: constante de rappel de l'oscillateur
        maxHarm: Harmonique supérieure à considérer
    """
    loadLine=lambda v: np.imag(driveImpedance(v,R,m,s)+tubeImpedance(length,freq2k(v),radius,radiativeImpedance(radius,freq2k(v))))
    roots=[]
    minfreq=1e-3
    #minfreq=(0.5)/4*airSoundSpeed()
    maxfreq=(2*(maxHarm)-1)/4*airSoundSpeed()
    #f=np.linspace(minfreq,maxfreq,300)
    #f=np.linspace(20,200,300)
    #plt.figure()
    #plt.plot(f,loadLine(f))
    sol=optimize.root_scalar(loadLine,bracket=[minfreq,maxfreq],method='brentq')
    roots.append(sol.root)
    return np.array(roots)
def averagePowerFactor(driveImpedance,loadImpedance):
    """
    Calculates the average power factor (in 1/N^2) of a loaded driver
    driverImpedance: the driver's complex mechanical impedance
    loadImpedance: the load's complex mechanical impedance

    taken from eq 10.3.9 (see 1.9.1 for ref to understand why the denominator has absolute value) of KinsleyFreyCoopSanders

    """
    Gamma=1/2*(np.real(loadImpedance))/np.absolute(driveImpedance+loadImpedance)**2
    return Gamma

if __name__ == "__main__":

    print("Densite de lair %.2e"%airDensity())
    print("Vitesse du son dans lair %.2e"%airSoundSpeed())
    print("Impédance charactéristique de lair %.2e"%airCharacImpedance())
    freq=30
    print("Une onde a %.2e Hz a une longueure d<onde de %.2e m et un nombre d<onde de %.2e rad m^-1"%(freq,freq2lambda(freq),freq2k(freq)))
    freqs=2*np.logspace(0,2,1000)
    tubeLength=6
    tubeRadius=0.25
    oscMass=1
    oscResistance=10
    oscSpring=36000
    print(drivenTubeResonnances(tubeLength,tubeRadius,oscResistance,oscMass,oscSpring))
    print("Osc. resonnance at %.2e Hz"%(1./(2*np.pi)*np.sqrt(oscSpring/oscMass)))
    Zt=lambda tubeLength: tubeImpedance(tubeLength,freq2k(freqs),tubeRadius)
    Zm=lambda oscMass,oscSpring: driveImpedance(freqs,oscResistance,oscMass,oscSpring)
    plt.figure()
    plt.plot(freqs,np.imag(Zm(1,36000)),'.',label="Heavy stiff oscillator")
    plt.plot(freqs,-np.imag(Zt(3)),'.',label="(-) Tube ($L=3$m)")
    plt.plot(freqs,np.imag(Zm(0.2,36000/5)),'.',label="Light flexible oscillator")
    plt.plot(freqs,-np.imag(Zt(6)),'.',label="(-) Tube ($L=6$m)")
    plt.plot(freqs,-np.imag(Zt(1)),'.',label="(-) Tube ($L=1$m)")
    plt.xscale('log')
    plt.grid('both')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Imag. Impedance (Rayls.m$^2$)')
    plt.legend()

    plt.figure()
    tubeLength=6
    tubeRadius=0.25
    oscMass=1
    oscResistance=10
    oscSpring=36000
    plt.plot(freqs,averagePowerFactor(driveImpedance(freqs,oscResistance,oscMass,oscSpring),tubeImpedance(tubeLength,freq2k(freqs),tubeRadius)),'-',label="(-) Tube ($L=%d$m), Heavy driver"%(tubeLength))
    tubeLength=3
    tubeRadius=0.25
    oscMass=1
    oscResistance=10
    oscSpring=36000
    plt.plot(freqs,averagePowerFactor(driveImpedance(freqs,oscResistance,oscMass,oscSpring),tubeImpedance(tubeLength,freq2k(freqs),tubeRadius)),'-',label="(-) Tube ($L=%d$m), Heavy driver"%(tubeLength))
    tubeLength=1
    tubeRadius=0.25
    oscMass=1/5
    oscResistance=10
    oscSpring=36000/5
    plt.plot(freqs,averagePowerFactor(driveImpedance(freqs,oscResistance,oscMass,oscSpring),tubeImpedance(tubeLength,freq2k(freqs),tubeRadius)),'-',label="(-) Tube ($L=%d$m), Heavy driver"%(tubeLength))
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True,which='both')
    plt.xlabel('Frequency (Hz)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power factor')
    plt.figure()
    tubeLength=6
    tubeRadius=0.25
    oscMass=1/5
    oscResistance=100
    oscSpring=36000/5
    plt.plot(freqs,averagePowerFactor(driveImpedance(freqs,oscResistance,oscMass,oscSpring),tubeImpedance(tubeLength,freq2k(freqs),tubeRadius)),'-',label="(-) Tube ($L=%d$m), Light driver"%(tubeLength))
    tubeLength=3
    plt.plot(freqs,averagePowerFactor(driveImpedance(freqs,oscResistance,oscMass,oscSpring),tubeImpedance(tubeLength,freq2k(freqs),tubeRadius)),'-',label="(-) Tube ($L=%d$m), Light driver"%(tubeLength))
    tubeLength=1
    plt.plot(freqs,averagePowerFactor(driveImpedance(freqs,oscResistance,oscMass,oscSpring),tubeImpedance(tubeLength,freq2k(freqs),tubeRadius)),'-',label="(-) Tube ($L=%d$m), Light driver"%(tubeLength))
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True,which='both')
    plt.xlabel('Frequency (Hz)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power factor')

    D=3*0.0254
    S=np.pi*(D/2)**2
    L=0.0254*5
    V=np.pi*(3/2*0.0254)**2*1
    Q=100
    s,m,R,w0=helmholtzParameters(S,V,L,Q)
    print("Spring constant %.2e, mass %.2e, resistance %.2e, Frequency %.2e"%(s,m,R,w0/(2*np.pi)))
    plt.figure()
    tubeLength=3
    plt.plot(freqs,averagePowerFactor(driveImpedance(freqs,R,m,s),tubeImpedance(tubeLength,freq2k(freqs),tubeRadius)),'-',label="(-) Tube ($L=%d$m), Helmholtz resonnator"%(tubeLength))
    plt.xscale('log')
    plt.legend()
    plt.yscale('log')
    plt.grid(True,which='both')
    plt.xlabel('Frequency (Hz)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power factor')
    plt.figure()
    Zm=lambda oscMass,oscSpring,oscResistance: driveImpedance(freqs,oscResistance,oscMass,oscSpring)
    plt.plot(freqs,np.imag(Zm(m,s,R)),'.',label="Helmholtz resonnator")
    plt.plot(freqs,-np.imag(Zt(3)),'.',label="(-) Tube ($L=3$m)")
    plt.xscale('log')
    plt.grid('both')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Imag. Impedance (Rayls.m$^2$)')
    plt.legend()
    plt.show()