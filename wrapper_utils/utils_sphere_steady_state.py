import numpy as np
from string import Template

ndim = 9

labels=["$\\left( p_1 - 1 \\right) \\log_{10}\\left(\\gamma_{min}\\right)$", \
        "$\\log_{10}\\left(\\Gamma \\cdot B \\cdot \\gamma_{break}^2\\right)$", \
        "$\\left( 3 - p_2 \\right) \\log_{10}\\left(\\gamma_{max}\\right)$", \
        "$p_1 - \\left( p_1 - 1 \\right) \\log_{10}\\left(\\gamma_{min}\\right) / 4$", \
        "$p_2 - \\log_{10}\\left(\\Gamma \\cdot B \\cdot \\gamma_{break}^2\\right) / 2$", \
        "$\\log_{10}\\left(B\\right)$", \
        "$\\log_{10}\\left(\\Gamma \\cdot B \\right)$", \
        "$\\log_{10}\\left(h \\cdot R^2 \\right)$", \
        "$\\log_{10}\\left(\\Gamma \\cdot B \\cdot R \\right)$", \
        "$\\log_{10}\\left(\\rho \cdot h \cdot \\sqrt[3]{\\frac{\\gamma_{min}^{2}} {\\Gamma \\cdot B}} \\right)$"]

initial_theta = [[-0.03,  37.25, 0.03], # nu_Fnu_opt
                 [-1.78,  -4.38, 1.73], # xray_opt
                 [-0.21, -11.98, 0.24], # gray_opt
                 [-0.28,  -0.78, 0.27], # B
                 [-0.58,   0.48, 0.52], # gamma B
                 [-0.13,   8.66, 0.16], # nu_br
                 [-2.80, -23.31, 2.63], # rho
                 [-0.89,  -4.21, 0.86], # R N
                 [-0.90,  10.68, 0.71]] # norm
initial_theta = np.array(initial_theta).T

config_template = Template("""
    [general]
    density = $dens
    magnetic_field = $mag
    dt = 1
    dt_max = $dt_max
    t_max = 1e7

    eta = 1

    gamma = $gamma # not used

    [volume]
    shape = "sphere"
    R = $radius

    [electrons]
    gamma_min = $g_min
    gamma_max = $g_max
    size      = 160
    distribution_type = "broken_power_law"
    break_point  = $g_break
    first_slope  = $p1
    second_slope = $p2

    [protons]
    size = 160
    distribution_type = "power_law"
    slope = 2

    [photons]
    size = 128
    """)

def lnprior(theta):
    gamma_min, gamma_break, gamma_max,  \
    p1, p2,                             \
    B, Gamma,                           \
    R, density = theta_to_params(theta)

    if (1e1 < gamma_min < gamma_break < gamma_max < 1e10 and

        1     < p1  < p2 < 6    and
        0.005 <  B       < 0.5  and
        2     <  Gamma   < 192  and
        1e10  <  R       < 1e26 and
        1e-26 <  density < 1e-16):
           return 0
    else:
        return -np.inf

def theta_to_config(theta):
    gamma_min, gamma_break, gamma_max,  \
    p1, p2,                             \
    B, Gamma,                           \
    R, density = theta_to_params(theta)

    config = config_template.substitute(
                g_min = gamma_min,
                g_break = gamma_break,
                g_max = gamma_max,
                p1 = p1,
                p2 = p2,
                mag = B,
                gamma = Gamma,
                radius = R,
                dens = density,
                # set dt_max to a tenth of the escape time (more or less)
                dt_max = R / 3e10 / 10)

    return config

def theta_to_params(theta):
    log_nuFnu_opt, \
    log_xray_opt,  \
    log_gray_opt,  \
    log_B,         \
    log_Gamma_B,   \
    log_nu_break,  \
    log_rho,       \
    log_R_N,       \
    log_norm = theta

    log_Gamma = log_Gamma_B - log_B
    log_g_break = (log_nu_break - log_Gamma_B) / 2

    aux = log_nuFnu_opt - 2 * log_nu_break

    log_R2 = aux - log_gray_opt - log_Gamma_B + log_g_break

    log_R = log_R2 / 2
    log_V = log_R  * 3
    log_N = log_R_N - log_R

    p1 = -(log_gray_opt - log_R_N + log_Gamma_B) / log_g_break

    p2 = p1 + log_norm      / log_g_break

    log_g_max = log_xray_opt / (3 - p2) + log_g_break
    log_g_min = (log_N - log_rho - np.log10(p1 - 1)) / (p1 - 1)

    return np.array([10**log_g_min,
                     10**log_g_break,
                     10**log_g_max,
                     p1,
                     p2,
                     10**log_B,
                     10**log_Gamma,
                     10**log_R,
                     10**log_rho])

def params_to_theta(params):
    gamma_min, gamma_break, gamma_max,  \
    p1, p2,                             \
    B, Gamma,                           \
    R, density = params

    N = density * (p1 - 1) * np.power(gamma_min, p1 - 1)
    nu_br = Gamma * B * gamma_break**2
    V = R**3

    flux_opt = N * V * nu_br**2 / np.power(gamma_break, p1 + 1)
    flux_x_ray_opt = np.power(gamma_max / gamma_break, 3 - p2)
    flux_g_ray_opt = R * N * np.power(gamma_break, -p1) / (Gamma * B)

    norm = np.power(gamma_break, p2 - p1)

    return np.array([np.log10(flux_opt),
                     np.log10(flux_x_ray_opt),
                     np.log10(flux_g_ray_opt),
                     np.log10(B),
                     np.log10(Gamma * B),
                     np.log10(nu_br),
                     np.log10(density),
                     np.log10(R * N),
                     np.log10(norm)])

def get_gamma(theta):
    log_nuFnu_opt, \
    log_xray_opt,  \
    log_gray_opt,  \
    log_B,         \
    log_Gamma_B,   \
    log_nu_break,  \
    log_rho,       \
    log_R_N,       \
    log_norm = theta

    return 10**(log_Gamma_B - log_B)

def get_R(theta):
    log_nuFnu_opt, \
    log_xray_opt,  \
    log_gray_opt,  \
    log_B,         \
    log_Gamma_B,   \
    log_nu_break,  \
    log_rho,       \
    log_R_N,       \
    log_norm = theta

    log_g_break = (log_nu_break - log_Gamma_B) / 2

    aux = log_nuFnu_opt - 2 * log_nu_break

    log_R2 = aux - log_gray_opt - log_Gamma_B + log_g_break
    log_R = log_R2 / 2

    return 10**(log_R)

def assert_consistency(theta):
    params = theta_to_params(theta)
    # Assert that theta_to_params and params_to_theta are actually inverse functions
    np.testing.assert_allclose(params_to_theta(params), theta)

    # Assert that you are not forgetting R or Gamma
    np.testing.assert_allclose(params[6], get_gamma(theta))
    np.testing.assert_allclose(params[7], get_R(theta))

    return True
