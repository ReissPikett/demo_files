Title: Cherenkov Radiation and Detection Model

ALMOST FINISHED - almost all functionality and documentation complete, just needs some tidying and fixing one minor bug

Description:

	   Section 1 - Predicting Number of Scintillation photons and their detection time.
	   
	   Calculating the number of photons predicted to be emmitted due to a muon, at given speed and charge,
	   passing through a scintillator, of given depth, radiation length and density. We then construct a Monte
	   Carlo algorithm that generates randomly directed photons emitted from all scintillated atoms due to the passing
	   Muon and simulate the internal reflections of the photons, until exiting through an open boundary at one end
	   of the scintillator. We then take the exit times of these photons, plot the data in a histogram.
	   ________________________________________________________

	   Section 2 - Predicting number of Cherenkov photons emitted and detected by single photon
           detector.

           The particles are assumed to travel roughly normal to the optical surface of the
           detector surface at a given speed. The model is based on the classical theory of Cherenkov 
           Radiation developed by I.M. Frank and I. E. Tamm 1934, defining the number of Cherenkov photons
           emmitted per wavelength interval and the Cherenkov Angle, and are given by the equations:
                  dN_emit/dλ =2παd*z^2/λ^2 (1-1/(nβ)^2)
                  theta_c=arccos(1/(nβ))

            where N_emit- number of photons emitted, λ - wavelength, d - optical depth travelled,
                  z- particle charge, n- refractive index of optic, β particle speed/c,
                  α- fine structure constant, theta_c- Cherenkov angle 

            This model also takes into account the quantum efficiency (QE), in relation to wavelength,
            of the photocathode material used in the detector. The final result is a 
            numerically calculated integral of the above differential equation multiplied by the quantum
            efficiency over each wavelength, between the minimum and maximum wavelengths which the photocathode
            is sensitive to, and is specified as:
                  N_det=2παd*z^2*sin^2(theta_c) ∫ QE(λ)/λ^2 dλ

	    where N_det - number of photons detected                   
            ________________________________________________________
                                     
            Section 3 - Simulating Cherenkov rings inside the detector.
           
            Simulate Cherenkov photons hitting single photon detector, binning the positions into 16x16 detector 
            channel grid, and generating best fits to characterise Cherenkov ring. When Cherenkov radiation
            occurs over short distances, such as those concerned here, when a particle passes through a medium,
            the number of emitted Cherenkov photons is constant in time. As the Cherenkov angle remains the same
            throughought the particle's travel, the perpendicular distance that the photon travels relative to
            the muon increases with increased distance after emission, or decreasing radial distance from the
            detector surface, as shown by the trigonometric relationship shown below:
 
                θ_c=arccos(ring_radius/muon_distance)
 
                where θ_c - Cherenkov angle,
                      ring_radius - radial distance of the photon from the muon after hitting detector,
                      ring_radius - distance travelled through optic window by muon
 
            This means that as the muon gets closer to the detector, the Cherenkov ring gets smaller, though the
            number of photons emitted remains the same, therefore we would expect a higher density of photons
            closer to position of the muon passting through the detector, and that density decreasing as we move
            outwards from that position.
 
            ________________________________________________________
 
            Section 4 - Characterising the Cherenkov Rings.
    
            After we have our total number of photons detected, and simulated their positions, we then take the
            coordinate data and bin them by a 16x16 channel grid to simulate the detector channels triggering as
            photon hits it. We and contructs a best fit ellipse over all channel hit data to define an 'average'
            Cherenkov ring, and constructs the maximum possible Cherenkov ring based on the data point with the
            largest distance from the mean data position. This maximum ring is then used to calculate the Cherenkov
            angle based on the radius and the depth travelled through the optic window, given by:
                        
                θ_c=arccos(max_radius/optic_depth)
 
                where θ_c - Cherenkov angle,
                      max_radius - maximum ring radius or maximum radial distance of the
                      photon from the muon after hitting detector,
                      optic_depth - depth of optical window


## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install relevant packages.

''' command line(windows/Linux)
pip install numpy
pip install matplotib
pip install pandas
pip install skimage
pip install scipy
'''

## Usage

Brief example and descriptions of script and contained functions. More information is available in the function docstrings.

'''python
import Muon_Cherenkov_library
YET TO BE COMPLETED

'''

