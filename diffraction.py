import numpy as np 
from FormFactor import FormFactor
import math
from crystals import Crystal, Atom, Element, LatticeSystem
import matplotlib.pyplot as plt
import mplcursors

# lattice_vectors = 3.35*np.eye(3)
# unitcell = [Atom('Po', coords = [0,0,0])]
# polonium = Crystal(unitcell, lattice_vectors)
# for atom in polonium:
#     print(atom)
FormFactor = FormFactor()
def calculate_structure_factor(atoms, h, k, l, wavelength, theta):
    sf = 0
    for atom in atoms:
        # element = Element(atom.element)
        element = atom.element
        form_factor = calculate_formFactor(theta, element, wavelength)
        #form_factor = 1
        x, y, z = atom.coords_fractional
        #print(x, y, z)
        sf += form_factor*np.exp(2*np.pi*1j*(h*x + k*y + l*z))
    # sf = 1 + (-1)**(h+k) + (-1)**(k+l) + (-1)**(h+l)
    return sf

def calculate_theta(h, k, l, a, b, c, alpha, beta, gamma, volume, wavelength, lattice):
    d = calculate_spacing(h, k, l, a, b, c, alpha, beta, gamma, volume, wavelength, lattice)
    theta = 2*np.arcsin(wavelength/2/d)*180/np.pi
    return theta, d

def calculate_spacing(h, k, l, a, b, c, alpha, beta, gamma, volume, wavelength, lattice):
    if lattice == LatticeSystem.cubic:
        d = 1/np.sqrt((h**2 + k**2 + l**2)/a**2)
    elif lattice == LatticeSystem.tetragonal:
        d = 1/np.sqrt((h**2+k**2)/a**2+l**2/c**2)
    elif lattice == LatticeSystem.orthorhombic:
        d = 1/np.sqrt(h**2/a + k**2/b + l**2/c)
    elif lattice == LatticeSystem.hexagonal:
        d = 1/np.sqrt(4/3*(h**2 + h*k + k**2)/a**2) + l**2/c**2
    elif lattice == LatticeSystem.rhombohedral:
        d = 1/np.sqrt((h**2 + k**2 + l**2)*np.sin(alpha)**2 + 2*(h*k + k*l + h*l)*np.cos(alpha)**2 - np.cos(alpha)/
        (a**2*(1-3*np.cos(alpha)**2 + 2*np.cos(alpha)**3)))
    elif lattice == LatticeSystem.monoclinic:
        d = 1/np.sqrt(1/np.sin(beta)**2*(h**2/a**2 + k**2*np.sin(beta)**2/b**2+l**2/c**2-2*h*l*np.cos(beta)/(a*c)))
    elif lattice == LatticeSystem.triclinic:
        s11 = b**2*c**2*np.sin(alpha)**2
        s22 = a**2*c**2*np.sin(beta)**2
        s33 = a**2*b**2*np.sin(gamma)**2
        s12 = a*b*c**2*(np.cos(alpha)*np.cos(beta)-np.cos(gamma))
        s23 = a**2*b*c*(np.cos(beta)*np.cos(gamma)-np.cos(alpha))
        s13 = a*b**2*c*(np.cos(gamma)*np.cos(alpha)-np.cos(beta))
        d = 1/volume**2*(s11*h**2 + s22*k**2 + s33*l**2 + 2*s12*h*k + 2*s23*k*l + 2*s13*h*l)
    return d

def calculate_formFactor(theta, element, wavelength):
    ff_coeffs = FormFactor.getCoeffs(element)
    a = ff_coeffs.a_coeffs
    b = ff_coeffs.b_coeffs
    c = ff_coeffs.c
    #print(a, b, c)
    orig_theta = theta
    theta = theta*np.pi/180
    q = 4*np.pi*np.sin(theta)/wavelength
    formFactor = c
    for i in range(4):
        formFactor += a[i]*np.exp(-1*b[i]*(q/(4*np.pi))**2)
    #print("theta: ", orig_theta, "q: ", 4*np.pi*np.sin(theta)/wavelength, "formFactor: ", formFactor)
    return formFactor

def calculate_peaks(crystal, wavelength):
    atoms = list(crystal.unitcell)
    a, b, c, alpha, beta, gamma = crystal.lattice_parameters
    volume = crystal.volume
    lattice = crystal.lattice_system
    intensity_dict = {}

    for h in range(-7,7):
        for k in range(-7,7):
            for l in range(-7,7):
                if h == 0 and k == 0 and l == 0:
                    continue
                two_theta, d = calculate_theta(h, k, l, a, b, c, alpha, beta, gamma, volume, wavelength, lattice)
                if not math.isnan(two_theta):
                    if d in intensity_dict:
                        # intensity_dict[d][0] = two_theta
                        #intensity_dict[d][1] += inten
                        # sity_dict[d][1]
                        intensity_dict[d][2] = (h,k,l)
                        intensity_dict[d][3] += 1
                    else:
                        theta = two_theta/180*np.pi/2
                        structure_factor = calculate_structure_factor(atoms, h, k, l, wavelength, two_theta/2)
                        lorentz_factor = (1+(np.cos(2*theta))**2)/((np.sin(theta))**2*np.cos(theta))
                        if structure_factor > 0.0001:
                            intensity_dict[d] = []
                            intensity_dict[d].append(two_theta)
                            intensity_dict[d].append(lorentz_factor*np.absolute(structure_factor)**2)
                            intensity_dict[d].append((h,k,l))
                            intensity_dict[d].append(1)
                            intensity_dict[d].append(np.real(structure_factor))
                            intensity_dict[d].append(np.imag(structure_factor))
                            intensity_dict[d].append(lorentz_factor)
    print(intensity_dict)
    return intensity_dict
def plot_diff_pattern(intensity_dict):
    max_intensity = max(intensity_dict[key][3]*intensity_dict[key][1] for key in intensity_dict)
    for d_key in intensity_dict:
        two_theta = intensity_dict[d_key][0]
        multiplicity = intensity_dict[d_key][3]
        rel_intensity = multiplicity*intensity_dict[d_key][1]/max_intensity
        hkl = intensity_dict[d_key][2]
        plt.axvline(two_theta, 0, rel_intensity, label=hkl)
    plt.xlabel("Angle (2$\Theta$)")
    plt.ylabel("Relative Intensity (counts)")
    #plt.xlim(15, 80)
    mplcursors.cursor(hover=True)
    plt.show()
    
# Cu = Crystal.from_mp(api_key="0qweAvX1EkrqqHzp", query = "mp-30") 
Cu = Crystal.from_cif("D:/Clement Research/xrd/cif/Cu.cif")
Si = Crystal.from_database('Si')
Diamond = Crystal.from_cif("D:/Clement Research/xrd/cif/diamond.cif")
NaMnFe = Crystal.from_cif("D:/Clement Research/xrd/cif/Na2Mn2Fe(VO4)3.cif")
print(list(Cu.unitcell))
# print(Cu)
#sf = calculate_structure_factor(Cu.unitcell, 1, 1, 1)
peaks = calculate_peaks(Cu, 1.54059)
#print(peaks)
plot_diff_pattern(peaks)
#print(peaks)
#sf = structure_factor()

