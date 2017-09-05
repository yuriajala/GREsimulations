# INICIALIZAÇÃO DAS BIBLIOTECAS
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d.axes3d import Axes3D

def set_rho(ff):
    rho_fat = np.array(0.5)
    rho_water = np.array((rho_fat*((1-ff)/ff)))
    return (rho_fat,rho_water)

def discrete_timeseries(n_echoes):
    '''
    This function calculates the GRE IP-OP signal from human liver, for any number of TEs,
    considering a low flip angle (Taylor approx. around zero) and 3T (TEs adjusted for 3T).
    This is based on studies by Hamilton and Yokoo.
    parameters: n_echoes     -> integer >0, number of equally spaced echoes
    returns:    echo_times
     '''
    
    #TE acquisition times for n echoes
    echo_times = np.linspace(0.0011515, 0.0011515*n_echoes, n_echoes, endpoint=True)
    
    return echo_times

def gen_GRE(echo_times, rho_fat, rho_water):
    '''
    This function calculates the GRE IP-OP signal from human liver, for any number of TEs,
    considering a low flip angle (Taylor approx. around zero) and 3T (TEs adjusted for 3T).
    This is based on studies by Hamilton and Yokoo.
    parameters: echo times -> array with TEs
                rho_fat    -> float >0, fat proton density
                rho_water  -> float >0, water proton density
    returns:    Sabs (absolute signal intensity) (equivalent to ydata)
    '''
    rho = np.array([rho_fat*0.047,\
    rho_fat*0.039,\
    rho_fat*0.006,\
    rho_fat*0.12,\
    rho_fat*0.7,\
    rho_fat*0.088,\
    rho_water]) #rho_fat is a general variable that multiplies intrinsic fat constants (Hamilton)
    
    #standard, individual T2* values
    t2star = np.array([0.030,\
              0.040,\
              0.051,\
              0.052,\
              0.062,\
              0.083,\
              0.023] )
    
    #frequency in parts per million for each moiety
    freq_ppm = np.array([5.3,\
             4.2,\
             2.75,\
             2.1,\
             1.3,\
             0.9,\
             4.7])

    freq = np.zeros(len(freq_ppm))
    magField = np.array(3.0)      #magnetic field, Tesla (if change, remember to change echo spacing)
    hydroGyro = np.array(42.57e6) #Hydrogen gyromagnetic ratio, Hz/T
    
    for n, fppm in enumerate(freq_ppm):
        freq[n] = fppm*magField*hydroGyro/1e6
    
    #calculate the complex signal for each moiety    
    S = (len(freq_ppm),len(echo_times))
    S = np.zeros(S,dtype='complex')
    
    for n, moieties in enumerate(freq_ppm):
        for t,TE in enumerate(echo_times):
            S[n,t] = rho[n]*np.exp(-TE/t2star[n])*np.exp(2.0*np.pi*1j*freq[n]*TE)         
    
    #obtain sum & absolute value
    Sabs= (len(echo_times))
    Sabs=np.zeros(Sabs,dtype='complex')
    
    Sabs = np.sum(S, axis=0)
    Sabs = np.absolute(Sabs)
    return (Sabs)

def continuous_timeseries(endpoint=0.01):
    return np.linspace(0, endpoint, 1000, endpoint=True)

def calc_fatfrac(fitted_parameters):
    ff =  fitted_parameters[0] / (fitted_parameters[0] + fitted_parameters[1])
    return ff

#################################################################################################

x = continuous_timeseries()  # plotar função contínua
ff = np.linspace(0.01,0.5)   # fat fraction
error = np.linspace(-20,20)    # erro para valores iniciais
fitted1 = []                 # lista com parametros encontrados no fit
        
for n in ff:
    rho_fat, rho_water = set_rho(n)
    xdata = discrete_timeseries(7)
    ydata = gen_GRE(xdata, rho_fat, rho_water)
    
    for erro in error:
        p0 = [rho_fat+erro, rho_water] # FAT + ERROR
        fitted_par, pcov = curve_fit(gen_GRE, xdata, ydata, p0=p0) 
        fitted1.append(calc_fatfrac(fitted_par))

# Cada linha do fitted1 é o parâmetro encontrado no curve_fit, cada coluna é um valor de erro
fitted1 = np.reshape(np.array(fitted1), (np.size(error), np.size(ff)))

#print(fitted1)

plt.figure(1)
plt.xlabel('error')
plt.ylabel('fitted')
for i in range(int(np.size(fitted1[0])/2)):
    plt.title('0-50%')
    plt.plot(error, ff-fitted1[i,:])
plt.show()

plt.figure(2)
plt.xlabel('error')
plt.ylabel('fitted')
for i in range(int(np.size(fitted1[0])/2)):
    plt.title('50-100%')
    plt.plot(error, ff-fitted1[(int(np.size(fitted1[0])/2)+i),:])
plt.show()

#############################################################################

x = continuous_timeseries()  # plotar função contínua
ff = np.linspace(0.01,0.5)   # fat fraction
error = np.linspace(-20,20)    # erro para valores iniciais
fitted1 = []                 # lista com parametros encontrados no fit
        
for n in ff:
    rho_fat, rho_water = set_rho(n)
    xdata = discrete_timeseries(7)
    ydata = gen_GRE(xdata, rho_fat, rho_water)
    
    for erro in error:
        p0 = [rho_fat+erro, rho_water] # FAT + ERROR
        fitted_par, pcov = curve_fit(gen_GRE, xdata, ydata, p0=p0) 
        fitted1.append(calc_fatfrac(fitted_par))

# Cada linha do fitted1 é o parâmetro encontrado no curve_fit, cada coluna é um valor de erro
fitted1 = np.reshape(np.array(fitted1), (np.size(error), np.size(ff)))

#print(fitted1)

plt.figure(1)
plt.xlabel('error')
plt.ylabel('fitted')
for i in range(int(np.size(fitted1[0])/2)):
    plt.title('0-50%')
    plt.plot(error, ff-fitted1[i,:])
plt.show()
plt.tight_layout()

plt.figure(2)
plt.xlabel('error')
plt.ylabel('fitted')
for i in range(int(np.size(fitted1[0])/2)):
    plt.title('50-100%')
    plt.plot(error, ff-fitted1[(int(np.size(fitted1[0])/2)+i),:])
plt.show()
plt.tight_layout()

#############################################################################

x = continuous_timeseries()  # plotar função contínua
ff = np.linspace(0.01,0.5)   # fat fraction
error = np.linspace(-20,20)    # erro para valores iniciais
fitted1 = []                 # lista com parametros encontrados no fit
        
for n in ff:
    rho_fat, rho_water = set_rho(n)
    xdata = discrete_timeseries(7)
    ydata = gen_GRE(xdata, rho_fat, rho_water)
    
    for erro in error:
        p0 = [rho_fat+erro, rho_water] # FAT + ERROR
        fitted_par, pcov = curve_fit(gen_GRE, xdata, ydata, p0=p0) 
        fitted1.append(calc_fatfrac(fitted_par))

# Cada linha do fitted1 é o parâmetro encontrado no curve_fit, cada coluna é um valor de erro
fitted1 = np.reshape(np.array(fitted1), (np.size(error), np.size(ff)))

#print(fitted1)

plt.figure(1)
plt.xlabel('error')
plt.ylabel('fitted')
for i in range(int(np.size(fitted1[0])/2)):
    plt.title('0-50%')
    plt.plot(error, ff-fitted1[i,:])
plt.show()

plt.figure(2)
plt.xlabel('error')
plt.ylabel('fitted')
for i in range(int(np.size(fitted1[0])/2)):
    plt.title('50-100%')
    plt.plot(error, ff-fitted1[(int(np.size(fitted1[0])/2)+i),:])
plt.show()


#############################################################################

fig = plt.figure(3,  tight_layout=True)
ax = fig.add_subplot(111, projection='3d')
cs = ['r','g','b','y','r','g','b','y','r','g','b','y']

fat_fract = list(np.linspace(0.01,0.5,100))
for ff in fat_fract:    
    x = continuous_timeseries()

    fitted_parameter3 = []

    erro = list(np.linspace(-100,100,100))

    for n in erro:
        rho_fat, rho_water = set_rho(ff)
        xdata = discrete_timeseries(7)
        ydata = gen_GRE(xdata, rho_fat, rho_water)

        p0 = [rho_fat + n, rho_water + n]
        fitted_par, pcov = curve_fit(gen_GRE, xdata, ydata, p0=p0) # valores iniciais + random em 0.01
        fitted_parameter3.append(fitted_par)

    fitted3 = np.array([calc_fatfrac(x) for x in fitted_parameter3])
    
    # x axis
    X = erro

    # y axis
    Y = np.zeros(100)+ff

    # data - plotted against z axis
    Z = ff-fitted3

    # populate figure
    ax.plot(np.array(X), np.array(Y), np.array(Z),'.')

ax.set_xlabel('Erro inserido')
ax.set_ylabel('FF real')
ax.set_zlabel('Erro na estimacao')
    # set initial view of plot
ax.view_init(elev=40., azim=345.)
plt.show()