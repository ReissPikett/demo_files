# Title: Cherenkov_Radiation_Model
# Description: Main Cherenkov Radiation study script, applying all functionality from my project library.
#
# Author: Reiss Pikett
# Date: 02/10/2020
#
###########################################################################################################

#import packages 
import Muon_Cherenkov_library
import numpy as np
import matplotlib.pyplot as plt

###########################################################################################################
#Constants and Parameters

# Scintillator parameters
full_scint_dimensions = Muon_Cherenkov_library.full_scint_dimensions # length of each dimension of scintillator [x,y,z]
tile_depth = Muon_Cherenkov_library.tile_depth # depth of one tile (cm)
num_tiles=Muon_Cherenkov_library.num_tiles # number of tiles per scintillator arrangement

e_muon=Muon_Cherenkov_library.e_muon  # energy of muon at earth's surface /eV
radiation_length=Muon_Cherenkov_library.radiation_length # material attenuation of energy. Energy/(density*depth) (E cm^2 g^-1)
scint_density =Muon_Cherenkov_library.scint_density # (g/cm^3)
photon_yield = Muon_Cherenkov_library.photon_yield # photons per MeV

# Particle parameters
beta=Muon_Cherenkov_library.beta# ratio velocity/speed of light 
z=Muon_Cherenkov_library.z

# detector optic parameters
optic_depth=Muon_Cherenkov_library.optic_depth #m travelled through detector optic
n_g=Muon_Cherenkov_library.n_g # refractive index of optic material

#####################################################################################################
# Main
#####################################################################################################

# Scintillation

# Scintillator photon calculations
energy_attenuation=Muon_Cherenkov_library.scintillator_absorption(radiation_length,scint_density,tile_depth,num_tiles)
e_1=e_muon-energy_attenuation
e_2=e_1-energy_attenuation
num_scint_photons=int(energy_attenuation*photon_yield)
print("Estimated number of scintillation photons emmitted per bar per passing Muon: "+str(int(num_scint_photons)))

full_bar_vertices = Muon_Cherenkov_library.full_bar_vertices
scint_vertices = Muon_Cherenkov_library.scint_vertices

# Scintillated atom positions, visualisation of photon exit path, all exit times of scintillation photons
scint_positions = Muon_Cherenkov_library.generate_scintillation_positions(num_scint_photons, random = True)
trial_path=Muon_Cherenkov_library.exit_time_and_path(scint_positions[0],full_scint_dimensions)
exit_times=Muon_Cherenkov_library.all_exit_times(scint_positions,full_scint_dimensions)

# Scintillation 3D plots
fig = plt.figure()
ax = plt.axes(projection="3d")
Muon_Cherenkov_library.plot_3D(scint_vertices,ax)
Muon_Cherenkov_library.plot_3D(full_bar_vertices,ax)
Muon_Cherenkov_library.plot_3D(trial_path[1],ax, lines = True)
Muon_Cherenkov_library.plot_3D(scint_positions,ax)
plt.show()

################################################################################################################

#Cherenkov photon position simulation and channel binning

#Chereknov angle calculations
cos_theta=1/(n_g*beta)  
cerAngle=np.arccos(cos_theta) # Cherenkov angle
max_radius= optic_depth*np.tan(cerAngle) # maximum radius of Cherenkov radiations

# Cherenkov Photons calculation and plots
x,y=Muon_Cherenkov_library.read_and_clean('Photek_Cathodes.xlsx', 'S20Q')

Muon_Cherenkov_library.plot_Cherenkov_Spect(x,y);

wavelengths,QuantumEfficiency = Muon_Cherenkov_library.interpolate_and_plot(x,y)

NumPhotons_NoQE_Calculated = Muon_Cherenkov_library.Cherenkov_Calc(wavelengths)

NumPhotons_NoQE_Integrated,NumPhotons_detected = Muon_Cherenkov_library.integrate_over_wavelength(wavelengths,QuantumEfficiency)


lamda1 = x[0] # Minimum wavelength
lamda2 = x[-1] # Maximum wavelength
print()
print(f' At length {optic_depth*1000}mm, refractive index {n_g}, muon speed {beta}c, over wavelengths: [{int(lamda1*1e9)}nm,{int(lamda2*1e9)}nm]')
print()
print(f"Cherenkov Angle: {np.round(cerAngle/np.pi*180,1)} degrees, Cherenkov Ring Maximum radius: {np.round(max_radius*1000,2)}mm")
print("Number of Cherenkov photons emitted: " + str(int(NumPhotons_NoQE_Integrated)))
print("Number of Cherenkov photons detected: " + str(int(NumPhotons_detected)))

#####################################################################################################

# Cherenkov Detection and Ring Classification

numPhotons=int(NumPhotons_detected) # estimated number of photons from Cherenkov Ring pred
origin = 8 # mean posistion detected from channels 0-15, both x and y coordinates
max_radius_pred=max_radius*1000 # predicted maximum cherenkov ring radius /mm
length_per_channel = Muon_Cherenkov_library.length_per_channel # length of 16 x 16 channels (15.873 mm) /mm 
max_radius_in_channels = max_radius_pred/length_per_channel


# Generating sample data and plotting with channel as grid
x,y = Muon_Cherenkov_library.random_photon_coordinates(origin,max_radius_in_channels,NumPhotons_detected)
Muon_Cherenkov_library.plot_data_with_channels(x,y)
plt.show()

#Plotting heat map of channel hits
Muon_Cherenkov_library.channel_heat_map(x,y,scatter=False)

channel_count_distribution = Muon_Cherenkov_library.channel_histogram(x,y)            
X,Y = Muon_Cherenkov_library.hits_HistogramToPositions(channel_count_distribution)    

# Calculating ellipse regression for unbinned photons
x_center_unbinned,y_center_unbinned,a_unbinned,b_unbinned,theta_t_unbinned = Muon_Cherenkov_library.least_squares_ellipse(x,y)
max_radius_in_channels_unbinned= Muon_Cherenkov_library.find_max_radius(x,y,x.mean(),y.mean())
max_radius_unbinned_in_mm = max_radius_in_channels_unbinned*length_per_channel

# Calculating ellipse regression for binned photons
x_center,y_center,a,b,theta_t=Muon_Cherenkov_library.least_squares_ellipse(X,Y)
max_radius_in_channels_calc = Muon_Cherenkov_library.find_max_radius(X,Y,X.mean(),Y.mean())
max_radius_in_mm = max_radius_in_channels_calc*length_per_channel

# Drawing best fit ellipse and circle of maximum radius for Cherenkov angle calc
ellipse_coords=Muon_Cherenkov_library.generate_ellipse(a,b,x_center,y_center,theta_t)
circle_coords=Muon_Cherenkov_library.generate_ellipse(max_radius_in_channels_calc,max_radius_in_channels_calc,x_center,y_center,0)

# Plotting binned photons, best fit ellips and max cherenkov ring
Muon_Cherenkov_library.plot_data_with_channels(X+0.5,Y+0.5)
plt.plot(x_center+ellipse_coords[0]+0.5 , y_center+ellipse_coords[1]+0.5, color='r')
plt.plot(x_center+circle_coords[0]+0.5 , y_center+circle_coords[1]+0.5, color='g')
plt.show()


print("Maximum radius of unbinned particles is: "+ str(np.round(max_radius_unbinned_in_mm,2))+"mm")
print("Cherenkov angle of unbinned particles is: ",np.round(np.arctan(max_radius_unbinned_in_mm/Muon_Cherenkov_library.optic_depth)/np.pi*180,2))
print("Maximum radius is: "+ str(np.round(max_radius_in_mm,2))+"mm")
print("Cherenkov angle is: ",np.round(np.arctan(max_radius_in_mm/(1000*Muon_Cherenkov_library.optic_depth))/np.pi*180,2))