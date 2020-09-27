# Title: Muon_Cherenkov_library
# Description: Library containing all functions and parameters used in my Cherenkov Radiation study
# project.
#
# Author: Reiss Pikett
# Date: 02/10/2020
#
###########################################################################################################

# Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from  skimage.measure import EllipseModel
import matplotlib.ticker as ticker
from scipy.interpolate import UnivariateSpline

###########################################################################################################

# Constants and Parameters
# Physical Constants
alpha = 0.0072973525693 # Fine structure constant - (2 pi e^2)/(hc)
c = 299792458 # speed of light in a vacuum (m/s)

# Particle parameters
z = -1  #Z particle charge number: multiples of e 
beta = 0.994 # velocity/speed of light in vacuum

# Device and optics parameters
# Main detector
optic_depth = 0.00888 # meters travelled through detector optic
n_g = 1.45 # refractive index of optic material - fused silica
length_per_channel = 15.873/16
orig = 7.5

#Single channel detector Parameters
#

# Scintillator parameters
# Physical properties
e_muon=4e9  #energy of muon at surface /eV
radiation_length=2e6 # Radiation Length - material attenuation of energy. Energy per (density * depth) / E cm^2 g^-1
scint_density = 1.2 # g/cm^3
tile_depth = 0.5 # depth of one tile/cm
num_tiles=3 # number of tiles per scintillator set-up
scint_dist=0.15 # 
photon_yield_MeV = 5000 # photons per MeV
photon_yield = photon_yield_MeV*1e-6

# Dimensions
full_scint_dimensions = [3,0.5,10]
full_bar_width = full_scint_dimensions[0]
full_bar_height = full_scint_dimensions[1]
full_bar_length = full_scint_dimensions[2]

scint_width = full_scint_dimensions[0]
scint_height = full_scint_dimensions[1]
scint_length = full_scint_dimensions[2]/6
scint_z_pos = full_scint_dimensions[2]*2/3 

# Calculated parameters
num_photons_coef= 2*np.pi*alpha*z**2*optic_depth*(1-1/(beta*n_g)**2) # collective term to represent product of 
# constants [m]. Represents the number of photons emitted over distance travelled over all wavelengths

##############################################################################

# Vertex coordinates of whole scintillator bar
full_bar_vertices = np.array([[0,0,0],[full_bar_width,0,0],[full_bar_width,full_bar_height,0],
                                [full_bar_width,full_bar_height,full_bar_length],
                                [0,full_bar_height,0],[0,full_bar_height,full_bar_length],
                                [0,0,full_bar_length], [full_bar_width,0,full_bar_length]])

# Vertex coordinates of shared area of two scintilators
scint_vertices = np.array([[0,0,scint_z_pos],[scint_width,0,scint_z_pos],[scint_width,scint_height,scint_z_pos],
                           [scint_width,scint_height,scint_length+scint_z_pos], [0,scint_height,scint_z_pos],
                           [0,scint_height,scint_length+scint_z_pos],[0,0,scint_length+scint_z_pos],
                           [scint_width,0,scint_length+scint_z_pos]])


###############################################################################

# function definitions and descriptions

def read_and_clean(
        data,photocathode):
    """ ELIM

    """
    df=pd.read_excel(data)
    df.drop(columns=['Unnamed: 0'],inplace=True)
    df.columns=df.iloc[1]
    df=df[2:]
    df.reset_index(inplace=True,drop=True)
    lam=str(df.columns[0])
    
    s20=df[[lam,photocathode]]
    s20['QM']=s20[photocathode]*124/s20[lam]
    s20.drop(columns=[photocathode],inplace=True)
    vals=s20.values.T
    x=vals[0].astype(float)*1e-9
    y=(vals[1].astype(float))/100
    return np.array([x,y])

def plot_Cherenkov_Spect(
        x,y):
    """shows a plot of the Cherenkov emission spectrum from the quantum efficiency data
    arguments:
    x -- independant variable data array
    y -- dependant valiable data array
    
    returns:
    void, displays line chart of y vs x in console 
    """
    
    plt.plot(x,num_photons_coef/x**2)
    plt.title("Cherenkov Emission Spectrum")
    plt.xlabel('Wavelength/m')
    plt.ylabel("Number of Photons per Wavelength (dN/dLambda)")
    plt.show()

def interpolate_and_plot(
        x,y,
        num_points=1000):
    """interpolates the Quantum Effiency (QE) data such that it can be same size as a linspace
    array. The interpolated QE is then plotted with respect to wavelength, and the function
    returns the linspace array and the interpolated QE.
    
    arguments:
    x -- independant variable data array
    y -- dependant valiable data array
    
    returns:
    x_line -- 1xnum_points array of equally distanced values ranging from the minimum and maximum of x argument
    QE -- 1xnum_points array of interpolated QE data as predictions of the univariate spline on the x_line values
    plot -- plots the original data points, a basic spline and the final smoothed spline of QE vs wavelength
    """
    
    univ_spline = UnivariateSpline(
        x, y)
    wavelength_line = np.linspace(x[0], x[-1], num_points)
    plt.plot(x, y, 'ro', ms=5)
    plt.plot(wavelength_line, univ_spline(wavelength_line), 'g', lw=3)
    univ_spline.set_smoothing_factor(0)
    QuantumEfficiency=univ_spline(wavelength_line)
    plt.plot(wavelength_line, QuantumEfficiency, 'b', lw=3)
    plt.xlabel('Wavelength/m')
    plt.ylabel("Quantum Efficiency /%")
    plt.title("Interpolated Quantum Efficiency data")
    plt.show()
    return np.array([wavelength_line,QuantumEfficiency])

def Cherenkov_Calc(
        x_line):
    """ Calculation of the number of photons emmitted due to Cherenkov Radiation as a function of
    multiple constants (specified in document description) and the min-max wavelengths in question.
    """
    
    lamda1=x_line[0]
    lamda2=x_line[-1]
    reciprocal_difference_calc=1/lamda1-1/lamda2
    N=num_photons_coef*reciprocal_difference_calc
    return int(N)

def integrate_over_wavelength(
        wavelength_line,QuantumEfficiency):
    """ Takes the quantum efficiencies of the detector at each wavelength and the wavelength
    emmision spectrum of the Cherenkov radiation, and integrates over all wavelengths to give
    total number of particle detected """
    
    integrand1=1/wavelength_line**2
    diff=wavelength_line[1]-wavelength_line[0]
    integ1=integrand1.cumsum()*diff
    N1=integ1[-1]*num_photons_coef
    # integrate QuantumEfficiency/lam^2
    integrand2=QuantumEfficiency/wavelength_line**2
    integ2=integrand2.cumsum()*diff
    N2=integ2[-1]*num_photons_coef
    return int(N1),int(N2)

def detector_proportional_area(
        width,distance):
    """ Calculating the proportion of emmited scintillator photons hitting the single channel detector. This assumes 
     an exposed scintillator of a given ditance from the detector free to emit light in all direactions. Calculates
     the proportinal area of the detector surface as a proportion the the total surface area of a sphere originated
     at the scintillator, where the radius is the distance from the scintillator to the detector. THis is to calculate
     the intensity of light emmitted from the scintillator that is incident on the detector surface. 
     
     Using an approximation:
     detector_area = 2*distance*width*arctan(width/2*distance), 
     theta = arctan(width/2*distance)
     
     arguments
     width = width of detector square surface (cm)
     distance=distance from scintillator to detector (cm)
     
     returns
     Proportion of the area of a square detector surface to the total surface area
     
     """
     
    theta=np.arctan(width/(2*distance)) 
    detector_area=2*distance*width*theta
    radiation_area=4*np.pi*distance**2
    return detector_area/radiation_area;

def scintillator_absorption(
        radiation_length,scint_density,
        tile_depth,num_tiles):
    """ 
    The total enery loss of a particle travelling through a given number of scintillator tiles of given depth (m),
    density (kg*m^-3), and radiation length(J*m^2/kg^2)
    
    arguments 
    """
    energy_attenuation = radiation_length*scint_density*tile_depth*num_tiles
    return energy_attenuation

def random_photon_coordinates(
        origin,radius, numPhotons):
    """ 
    Generates random set of cartesian coordinates corresponding to Cherenkov photon detection. A set of polar coordinates
    are generated about a given origin, up to a maximum radius (minimum being zero), and a given number of data points, 
    corresponding to number of photons detected, and converted to cartesian coorindates. Used to simulate Cherenkov
    photons and the detection of Cherenkov rings.
    
    arguments
    origin = origin position, taken in form [x,y]
    radius = maximum radius of random polar coordinates
    numPhotons = number of coordinates generated/length of dataset, corresponding to number of photons hitting detector
    
    returns
    array of random coordinates 
    """
    phi=np.random.uniform(0,2*np.pi,numPhotons)
    r=np.random.uniform(0,radius,numPhotons)
    pos =polar_to_cartesian(phi, r)
    return pos[0]+origin,pos[1]+origin

def polar_to_cartesian(
        theta, r):
    """ converts array of polar coordinates (with origin= [0,0]) to array x,y coordinates
    
    arguments
    theta = polar angle of coordinate
    r = radius of polar coordinate
    """
    pos= np.asarray([r[i]*np.array([np.cos(theta[i]),np.sin(theta[i])]) for i in range(len(r))]).T
    return pos[0],pos[1]

def plot_data_with_channels(
        x,y):
    """ plot data in a 16x16 grid. x,y data displayed as inputted, the grid signifies the channel layout """
    
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1)
    major_ticks = np.arange(0, 16, 1)
    ax.set_xticks(major_ticks)
    ax.set_yticks(major_ticks)
    ax.grid(which='major', alpha=0.5)
    # Hide major tick labels
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    # Customize minor tick labels
    ax.xaxis.set_minor_locator(ticker.FixedLocator(major_ticks+0.5))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(major_ticks))
    ax.yaxis.set_minor_locator(ticker.FixedLocator(major_ticks+0.5))
    ax.yaxis.set_minor_formatter(ticker.FixedFormatter(major_ticks))
    plt.scatter(x,y)
    ax.set_xlim([0,16])
    ax.set_ylim([0,16])
    ax.set_xlabel('X axis Channels')
    ax.set_ylabel('Y axis Channels')
    
def channel_histogram(
        x,y):
    """ returns histogram array of channel hit frequency """
    
    xedges=[i-0.5 for i in range(17)]
    yedges=xedges
    channel_count_distribution,xedges,yedges=np.histogram2d(x,y,bins=[xedges,yedges])
    channel_count_distribution=channel_count_distribution.T
    return channel_count_distribution

def hits_HistogramToPositions(
        dist):
    """converts the numpy histogram of hit freq vs channel positions (16 x 16) into an array of
     channel positions of each hit (numPh x numPh, [X,Y]) """
     
    grid=[]
    for i in range(16):
        for j in range(16):
                if dist[i][j]>0:
                    for k in range(int(dist[i][j])):
                        grid.append([j,i])
    return np.asarray(grid).T*1.0
    
def channel_heat_map(
        x,y,scatter = True):
    """ Displayes 2d histogram heatmap to signify channel hit frequency, + unbinned data if scatter==True 
    params x, y, 
    
    """
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1)
    major_ticks = np.arange(0, 16, 1)
    ax.set_xticks(major_ticks)
    ax.set_yticks(major_ticks)
    if scatter == True:
        plt.scatter(x,y)
    dist = channel_histogram(x,y)
    ax.set_xlabel('X axis Channels')
    ax.set_ylabel('Y axis Channels')
    plt.imshow(dist)


def least_squares_ellipse(
        x,y):
    """ Uses Ordinary Least Squares to generate ellipse parameters for best fit ellipse of data """
    
    xy = np.array([x,y]).T
    ellipse =EllipseModel()
    ellipse.estimate(xy)
    xc,yc,a,b,theta_t=ellipse.params
    return xc,yc,a,b,theta_t

def find_max_radius(
        X,Y,xc,yc):
    """ Takes channel position data and returns the maximum radius of possible rings, in channel number and in meters """
    
    max_radius=np.sqrt(((X-xc)**2+(Y-yc)**2).max())
    return max_radius

def generate_ellipse(
        a,b,xc,yc,
        theta_t):
    """ gives coordinates of created ellipse drawn from circular linespace, based on OLS fit """
    
    anglespace = np.linspace(0, 2*np.pi, 100)
    ellipse_coords=[a*np.cos(anglespace) ,b*np.sin(anglespace)]
    rotation=np.array([[np.cos(theta_t),-np.sin(theta_t)],
                      [np.sin(theta_t),np.cos(theta_t)]])
    ellipse_coords=np.dot(rotation,ellipse_coords)
    return ellipse_coords

######################################################################################################################

# Scintillator Internal reflection

def random_light_direction():
    """
    Generates a vector representing light, traveling in a random direction. Creates a random 3D unit
    vector (direction) with a uniform spherical distribution, and multiplies by he speed of light in a vacuum (m/s).
    
    Arguments
    ---------
    void
    
    Return
    ---------
    3d velocity vector of light speed and direction
    
    """
    z=0
    while(z==0):
        phi = np.random.uniform(0,np.pi*2)
        costheta = np.random.uniform(-1,1)
    
        theta = np.arccos( costheta )
        x = np.sin( theta) * np.cos( phi ) * c
        y = np.sin( theta) * np.sin( phi ) * c
        z = np.cos( theta ) * c
    return np.array([x,y,z])

def first_border_time(pos,vel,full_scint_dimensions):
    """
    Calculates the time for the photon, traveling at vel, starting as pos, to arrive at its next border.
    Determines the direction of travel, and calculates for each dimension the time taken for the photon
    to reach each respected border. The minimum of these times is then returned.
    
    Arguments
    ----------
    pos - starting position of the photon
    vel - velocity vector of photon
    
    Returns
    -------
    min_time - time taken for photon to reach first border 
    

    """
    
    if vel[0]>0:
        t_x = (full_scint_dimensions[0]-pos[0])/np.abs(vel[0])
    else:
        t_x = pos[0]/np.abs(vel[0])
    if vel[1]>0:
        t_y  =(full_scint_dimensions[1]-pos[1])/np.abs(vel[1])
    else:
        t_y = pos[1]/np.abs(vel[1])
    if vel[2]>0:
        t_z = (full_scint_dimensions[2]-pos[2])/np.abs(vel[2])
    else:
        t_z=pos[2]/np.abs(vel[2])
    min_time= min(t_x,t_y,t_z)
    return min_time

def reflect(pos,vel,full_scint_dimensions):
    """
    Takes the position and velocity of a photon that has arrived at a border and returns a velocity
    vector of the light as reflected from that border, where the photon velocity dimension that is
    normal to the border surface. It is assumed that there is no scattering and so that the angle
    of incidence is equal to angle of reflection.
    
    Arguments
    ----------
    pos - position of photon, given that it is at a border
    vel - velocity vector of photon at arrival
    
    Returns
    -------
    vel - velocity vector of reflected photon
    
    """
    if pos[0]==0 or pos[0]==full_scint_dimensions[0]:
        vel[0] = (-1.0)*vel[0]
    if pos[1]==0 or pos[1]==full_scint_dimensions[1]:
        vel[1] = (-1.0)*vel[1]
    if pos[2]==0 or pos[2]==full_scint_dimensions[2]:
        vel[2] = (-1.0)*vel[2]
    return vel


def exit_time(position, full_scint_dimensions):
    """
    Calculates the exit time of a photon emitted from a given position. Takes the starting position,
    generates a velocity vector representing the emitted photon, calculates least time to nearest border,
    moves photon to nearest border, reflects off the border by reversing velocity direction normal to border,
    and repeats until the photon has left the scintillator (position in z-direction is less than zero). Returns -
    
    Arguments
    ----------
    pos - position of excited scintillator atom, representing the starting position of a photon.
    
    Returns
    -------
    time_accumulated - total time taken for the photon to leave the scintillator through open boundary at z=0.
    also returns -1 if number of reflections becomes too many
    
    """
    light = random_light_direction()
    time_accumulated = 1.4e-9
    num = 0
    while(position[2]>0):
        time = first_border_time(position, light, full_scint_dimensions)
        time_accumulated+=time
        position = position+time*light
        light = reflect(position,light,full_scint_dimensions)
        num+=1
        if num >100:
            return -1
    return time_accumulated 

def exit_time_and_path(position, full_scint_dimensions):
    """
    Similar function to exit_time, but intended purely for visual purposes for one emitted photon, and 
    will always return a completed path where the photon exits the scinillator. Calculates the exit time of
    a photon emitted from a given position, as well as each point the the photon hits a border. Takes the
    starting position, generates a velocity vector representing the emitted photon, calculates least time to
    nearest border, moves photon to nearest border, reflects off the border by reversing velocity direction
    normal to border, and repeats until the photon has left the scintillator (position in z-direction is less
    than zero).
    
    Arguments
    ----------
    pos - position of excited scintillator atom, representing the starting position of a photon.
    
    Returns
    -------
    time_accumulated - total time taken for the photon to leave the scintillator through open boundary at z=0.
    all_paths - positions of photon path at each point of path when hitting a border.
    
    """
    light = random_light_direction()
    time_accumulated = 1.4e-9
    num = 0
    changing_position = position
    path=[changing_position]
    while(changing_position[2]>0):
        time = first_border_time(changing_position, light, full_scint_dimensions)
        time_accumulated+=time
        path.append(changing_position+time*light)
        changing_position=changing_position+time*light
        light = reflect(changing_position,light,full_scint_dimensions)
        num+=1
        if num >200:
            changing_position=position
            time_accumulated, allpaths = exit_time_and_path(position, full_scint_dimensions)
            return  time_accumulated, allpaths
    all_paths = np.asarray(path)
    return time_accumulated, all_paths


def all_exit_times(scint_points, full_scint_dimensions):
    """
    Calculates exit times of all photons emmitted by excited scintillator atoms. Takes array of starting
    positions of emitted photons, returns an array of the exit time of each photon and plots histogram of
    all exit times.
    
    Arguments
    ----------
    scint_points - array of positions of each excited scintillator atom
    
    Returns
    -------
    exit_times - array of exit times for each photon emitted from excited atom.
    
    """
    exit_times = []
    for scint_point in scint_points:
        time = exit_time(scint_point, full_scint_dimensions)
        exit_times.append(time)
    exit_times = np.asarray(exit_times)
    exit_times = exit_times[exit_times!=(-1)]
    plt.hist(exit_times*1e9,bins = 30)
    plt.xlabel('Time of arrival after of Muon passing /ns')
    plt.ylabel('Scintillation Photons Detected')
    return exit_times

def plot_3D(data, ax, lines=False):
    """
    Takes three dimentional data and plots 3D plot of data, with option of lines joining data points. 
    
    Arguments
    ----------
    data - array of 3d input data
    lines - Boolean turn on lines to join data points
    
    Returns
    -------
    void
    
    """
    x_line = data.T[2]
    y_line = data.T[0]
    z_line = data.T[1]
    if lines:
        ax.plot3D(x_line, y_line, z_line)
    ax.scatter3D(x_line, y_line, z_line, cmap='hsv')
    
def generate_scintillation_positions(
        num_scint_photons,
        x=np.random.uniform(0,scint_width),
        z=np.random.uniform(scint_z_pos,scint_z_pos+scint_length),
        random = False):
    """
     Generate points of scintillated atoms in area of scintillator that has shared area with other
     scintillator below. Simulates a muon passing through the scintillator, exciting all atoms in the scintillator
     in that axis. Allows for input x,z coordinates or will generate random x,z coordinates. 
    
    Arguments
    ---------
    num_scint_photons - total number of scintillation photons
    x - x position of muon as it passes through scintillator
    z - z position of muon as it passes through scintillator
    random - set to true generates a random x,z position
    
    Returns
    -------
    scint_positions - position array of all scintillated atom coordinates
    """
    muonPos = np.array([x,0,z])
    x_perturb = np.random.normal(0,0.01, num_scint_photons)
    z_perturb = np.random.normal(0,0.01, num_scint_photons)

    scint_x = muonPos[0]+x_perturb
    scint_y = np.random.uniform(0, scint_height, num_scint_photons)
    scint_z = muonPos[2]+z_perturb
    scint_positions = np.array([scint_x,scint_y,scint_z]).T
    return scint_positions

