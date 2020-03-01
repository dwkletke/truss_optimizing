from cmath import pi
import numpy
import matplotlib.pyplot as plt

thickness = 0.0033
density_board = 141
f_ratio_1 = 20 / 12.81
f_ratio_2 = 12.81 / 20
f_ratio_3 = 12.81 / 20
f_ratio_4 = 10.2 / 20
ten_stress = 24.86e6
length_1 = 0.15  # Compressive
length_2 = 0.096  # Compressive
length_3 = 0.1601  # Tension
length_4 = 0.204  # Tension
num_dowels = 5
dowel_mass = 0.00041 * num_dowels
hole_size = 0.0033
E = 3.39e9
g = 9.81


def find_pv_ixx(w, member):
    if member == 1:
        b1 = find_buckling(length_1, find_ixx(w, w + thickness))
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


def find_pv_iyy(w, member):
    if member == 1:
        b1 = find_buckling(length_1, find_iyy(w, w + thickness))
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


def find_mass_for_plot(w, ibeam):
    if ibeam:
        b1 = find_buckling(length_1, find_iyy(w, w))
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


def total_mass(m1, m2, m3, m4):
    return (2 * m1) + (2 * m2) + (2 * m3) + (2 * m4) + dowel_mass


def find_mass_balsa(w, l, is_beam):
    if is_beam:
        return 2 * (pi * (((w + hole_size) / 2) ** 2) * thickness) - (pi * ((hole_size / 2) ** 2) * thickness) + (
                w * thickness * l * density_board) + 2 * ((w + thickness) * thickness * l * density_board)
    else:
        return 2 * (pi * (((w + hole_size) / 2) ** 2) * thickness) - (pi * ((hole_size / 2) ** 2) * thickness) + (
                w * thickness * l * density_board)


def ibeam_strength(w, l):
    return find_buckling(l, find_iyy(w, w + thickness))


def comp_strength(w, l):
    return find_buckling(l, get_inertia(w))

def calc_width(f, f_ratio):
    return f * f_ratio / (thickness * ten_stress)


def calc_width_from_buck(f, l):
    return (f * (l ** 2) * 12) / ((pi ** 2) * E * (thickness ** 3))


def calc_width_for_ibeam(f, l):  # for members experiencing buckling with i-beam
    return (f * (l ** 2) * 12) / ((pi ** 2) * E * (thickness ** 3))


def get_inertia(w):
    return (w * (thickness ** 3)) / 12


def find_buckling(l, inertia):
    return (pi ** 2) * E * inertia / (l ** 2)


def find_ixx(w1, w2):
    return ((w1 ** 3) * thickness) / 12 + 2 * (
                (((thickness ** 3) * w2) / 12) + (((thickness * w2) * ((w1 + thickness) ** 2)) / 4))


def find_iyy(w1, w2):
    return (((thickness ** 3) * w1) / 12) + (2 * (((w2 ** 3) * thickness) / 12))


def buckling_check_top(w):
    return find_buckling(length_1, find_iyy(w, w))


def buckling_check_side(w):
    return find_buckling(length_2, ((w * (thickness ** 3)) / 12.0))


def get_tension(w, l):
    return ten_stress * w * l


def plot_width_vs_force():
    buckling_values_top = list()
    buckling_values_side = list()
    tension_member_3 = list()
    tension_member_4 = list()
    x = numpy.arange(0.001, 0.01, 0.0001)
    for i in x:
        buckling_values_top.append(buckling_check_top(i))
        buckling_values_side.append(buckling_check_side(i))
        tension_member_3.append(get_tension(i, length_3))
        tension_member_4.append(get_tension(i, length_4))
    plt.plot(x, buckling_values_top, 'b', label="Top Compression")
    # plot(x, buckling_values_side, 'r', label="Side Compression")
    # plt.plot(x, tension_member_3, 'g', label="Top Tension")
    # plt.plot(x, tension_member_4, 'y', label="Bottom Tension")
    plt.legend(loc="upper left")
    plt.xlabel('Width [m]')
    plt.ylabel('Force [N]')
    plt.title('Width vs. Max Force (All Members)')
    plt.show()


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


def plot_width_vs_mass():
    mass_values = list()
    x = numpy.arange(0.001, 0.01, 0.0001)
    for i in x:
        mass_values.append(find_mass_for_plot(i))
    plt.plot(x, mass_values, 'b', label="Total Mass")
    plt.legend(loc="lower right")
    plt.xlabel('Width [m]')
    plt.ylabel('Total mass [kg]')
    plt.title('Width vs. Mass (I-Beam)')
    plt.show()


def plot_mass_vs_ibeam_vs_comp():
    mass_values = list()
    ibeam_strength_values = list()
    comp_strength_values = list()

    x = numpy.arange(0.001, 0.01, 0.0001)
    for i, val in enumerate(x):
        mass_values.append(find_mass_balsa(val, length_1, True))
        ibeam_strength_values.append(ibeam_strength(val, length_1))
        comp_strength_values.append(comp_strength(val, length_1))

    plt.plot(mass_values, ibeam_strength_values, 'b', label="I Beam")
    plt.plot(mass_values, comp_strength_values, 'r', label="Comp.")
    plt.legend(loc="upper left")
    plt.xlabel('Mass [kg]')
    plt.ylabel('Force [N]')
    plt.title('Mass vs. Force (I-Beam and Comp.)')
    plt.show()



plot_mass_vs_ibeam_vs_comp()


# Find width from width vs force for exact number