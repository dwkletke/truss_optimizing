import numpy
import matplotlib.pyplot as plt

# Constants
thickness = 0.0033
density_board = 141
f_ratio_1 = 20 / 12.81
f_ratio_2 = 12.81 / 20
f_ratio_3 = 12.81 / 20
f_ratio_4 = 10.2 / 20
ten_stress = 24.86e6
length_1 = 0.25  # Compressive
length_2 = 0.096  # Compressive
length_3 = 0.1601  # Tension
length_4 = 0.204  # Tension
num_dowels = 5
dowel_mass = 0.00041 * num_dowels
hole_size = 0.0033
E = 3.39e9
g = 9.81


# Finding PV based on failure for a compressive member either with Ixx inertia for a I beam or a regular inertia
def find_pv_ixx(w, member, beam):
    if member == 1:
        b1 = find_buckling(length_1, find_ixx(w, w + thickness), beam)
        m1 = find_mass_balsa(w, length_1, True)
        w2 = calc_width_from_buck(b1 * f_ratio_2, length_2)
        m2 = find_mass_balsa(w2, length_2, False)
    else:
        b1 = find_buckling(length_2, get_inertia(w))
        m1 = find_mass_balsa(w, length_2, False)
        w2 = calc_width_from_buck(b1 * f_ratio_1, length_1)
        m2 = find_mass_balsa(w2, length_1, False)
    w3 = calc_width(b1, f_ratio_3)
    m3 = find_mass_balsa(w3, length_3, False)
    w4 = calc_width(b1, f_ratio_4)
    m4 = find_mass_balsa(w4, length_4, False)
    return (b1 / g) / total_mass(m1, m2, m3, m4)


# Finding PV based on failure for a compressive member either with Iyy inertia for a I or T beam or a regular inertia
def find_pv_iyy(w, member, beam):
    if member == 1:
        b1 = find_buckling(length_1, find_iyy(w, w + thickness), beam)
        m1 = find_mass_balsa(w, length_1, True)
        w2 = calc_width_from_buck(b1 * f_ratio_2, length_2)
        m2 = find_mass_balsa(w2, length_2, False)
    else:
        b1 = find_buckling(length_2, get_inertia(w))
        m1 = find_mass_balsa(w, length_2, False)
        w2 = calc_width_from_buck(b1 * f_ratio_1, length_1)
        m2 = find_mass_balsa(w2, length_1, False)
    w3 = calc_width(b1, f_ratio_3)
    m3 = find_mass_balsa(w3, length_3, False)
    w4 = calc_width(b1, f_ratio_4)
    m4 = find_mass_balsa(w4, length_4, False)
    return (b1 / g) / total_mass(m1, m2, m3, m4)


# Finding mass of truss based on beam type and width of failure member
def find_mass_for_plot(w, beam):
    if beam == "t" or "i":
        b1 = find_buckling(length_1, find_iyy(w, w), beam)
        m1 = find_mass_balsa(w, length_1, True)
        w2 = calc_width_from_buck(b1 * f_ratio_2, length_2)
        m2 = find_mass_balsa(w2, length_2, False)
    else:
        b1 = find_buckling(length_2, get_inertia(w))
        m1 = find_mass_balsa(w, length_2, False)
        w2 = calc_width_from_buck(b1 * f_ratio_1, length_1)
        m2 = find_mass_balsa(w2, length_1, False)
    w3 = calc_width(b1, f_ratio_3)
    m3 = find_mass_balsa(w3, length_3, False)
    w4 = calc_width(b1, f_ratio_4)
    m4 = find_mass_balsa(w4, length_4, False)
    return total_mass(m1, m2, m3, m4)


# Calculates total mass of truss
def total_mass(m1, m2, m3, m4):
    return (2 * m1) + (2 * m2) + (2 * m3) + (2 * m4) + dowel_mass


# Calculates mass of member (balsa wood) based on width and beam
def find_mass_balsa(w, l, beam):
    if beam == "i":
        return 2 * (numpy.pi * (((w + hole_size) / 2) ** 2) * thickness) - (
                numpy.pi * ((hole_size / 2) ** 2) * thickness) + (
                       w * thickness * l * density_board) + 2 * ((w + thickness) * thickness * l * density_board)
    elif beam == "t":
        return 2 * (numpy.pi * (((w + hole_size) / 2) ** 2) * thickness) - (
                numpy.pi * ((hole_size / 2) ** 2) * thickness) + (
                       w * thickness * l * density_board) + ((w + thickness) * thickness * l * density_board)
    else:
        return 2 * (numpy.pi * (((w + hole_size) / 2) ** 2) * thickness) - (
                numpy.pi * ((hole_size / 2) ** 2) * thickness) + (
                       w * thickness * l * density_board)


# Calculates the max force of the beam member based on a width
def beam_strength(w, l):
    return find_buckling(l, find_iyy(w, w + thickness))


# Calculates the max force of the regular compressive member
def comp_strength(w, l):
    return find_buckling(l, get_inertia(w))


# Calculates the width of a member based on max force and ratio
def calc_width(f, f_ratio):
    return f * f_ratio / (thickness * ten_stress)


# Calculates width of a member from max buckling force
def calc_width_from_buck(f, l):
    return (f * (l ** 2) * 12) / ((numpy.pi ** 2) * E * (thickness ** 3))


# Calculates width for beam from buckling force
def calc_width_for_beam(f, l):
    return (f * (l ** 2) * 12) / ((numpy.pi ** 2) * E * (thickness ** 3))


# Calculates inertia of regular compressive member
def get_inertia(w):
    return (w * (thickness ** 3)) / 12


# Calculates max buckling force based on length and inertia
def find_buckling(l, inertia):
    return (numpy.pi ** 2) * E * inertia / (l ** 2)


# Find Ixx will always be larger than Iyy for our purpose, so we only care about Iyy
def find_ixx(w1, w2, beam):
    if beam == "i":
        return ((w1 ** 3) * thickness) / 12 + 2 * (
            (((thickness ** 3) * w2) / 12) + (((thickness * w2) * ((w1 + thickness) ** 2)) / 4))
    elif beam == "t":
        return 0


# Find Iyy Inertia (is the same for a I and a T beam)
def find_iyy(w1, w2):
    return (((thickness ** 3) * w1) / 12) + (2 * (((w2 ** 3) * thickness) / 12))


# Finds tension force based on width and length
def get_tension(w, l):
    return ten_stress * w * l


# Plots width vs force based on input member as failure mode
def plot_width_vs_force(member):
    buckling_values_top = list()
    buckling_values_side = list()
    tension_member_3 = list()
    tension_member_4 = list()
    x = numpy.arange(0.001, 0.01, 0.0001)
    for i in x:
        buckling_values_top.append(find_buckling(length_1, find_iyy(i, i + thickness)))
        buckling_values_side.append(find_buckling(length_2, get_inertia(i)))
        tension_member_3.append(get_tension(i, length_3))
        tension_member_4.append(get_tension(i, length_4))
    if member == 1:
        plt.plot(x, buckling_values_top, 'b', label="Top Compression Member")
    elif member == 2:
        plt.plot(x, buckling_values_side, 'r', label="Side Compression Member")
    elif member == 3:
        plt.plot(x, tension_member_3, 'g', label="Top Tension Member")
    else:
        plt.plot(x, tension_member_4, 'y', label="Bottom Tension Member")
    plt.legend(loc="upper left")
    plt.xlabel('Width [m]')
    plt.ylabel('Force [N]')
    plt.title('Width vs. Max Force')
    plt.show()


# Plots width vs pv based on input member as failure mode
def plot_width_vs_pv(member):
    pv_values_ixx = list()
    pv_values_iyy = list()
    x = numpy.arange(0.001, 0.01, 0.0001)
    for i in x:
        pv_values_ixx.append(find_pv_ixx(i, member))
        pv_values_iyy.append(find_pv_iyy(i, member))
    plt.plot(x, pv_values_ixx, 'b', label="Ixx")
    plt.plot(x, pv_values_iyy, 'r', label="Iyy")
    plt.legend(loc="lower right")
    plt.xlabel('Width [m]')
    plt.ylabel('PV [kg/kg]')
    plt.title('Width vs. PV (I-Beam)')
    plt.show()


# Plots width vs mass based on input beam
def plot_width_vs_mass(beam):
    mass_values = list()
    x = numpy.arange(0.001, 0.01, 0.0001)
    for i in x:
        mass_values.append(find_mass_for_plot(i, beam))
    plt.plot(x, mass_values, 'b', label="Total Mass")
    plt.legend(loc="lower right")
    plt.xlabel('Width [m]')
    plt.ylabel('Total mass [kg]')
    plt.title('Width vs. Mass (I-Beam)')
    plt.show()


# Plots mass vs force of  based on input beam
def plot_mass_vs_ibeam_vs_comp():
    mass_values = list()
    ibeam_strength_values = list()
    comp_strength_values = list()

    x = numpy.arange(0.001, 0.01, 0.0001)
    for i, val in enumerate(x):
        mass_values.append(find_mass_balsa(val, length_1, True))
        ibeam_strength_values.append(beam_strength(val, length_1))
        comp_strength_values.append(comp_strength(val, length_1))

    plt.plot(mass_values, ibeam_strength_values, 'b', label="I Beam")
    plt.plot(mass_values, comp_strength_values, 'r', label="Comp.")
    plt.legend(loc="upper left")
    plt.xlabel('Mass [kg]')
    plt.ylabel('Force [N]')
    plt.title('Mass vs. Force (I-Beam and Comp.)')
    plt.show()


plot_width_vs_force(1)

print(calc_width_from_buck(38.5, length_2))
# Find width from width vs force for exact number
