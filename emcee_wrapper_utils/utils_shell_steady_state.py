import numpy as np
from string import Template

ndim = 10

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

initial_theta = [[-0.03,  37.33, 0.03], # nu_Fnu_opt
                 [-1.78,  -5.51, 1.73], # xray_opt
                 [-0.21, -14.06, 0.24], # gray_opt_h_N/2
                 [-0.28,  -1.36, 0.27], # B
                 [-0.58,  -0.63, 0.52], # gamma B
                 [-0.13,   8.57, 0.16], # nu_br
                 [-2.80, -21.79, 2.63], # rho
                 [-2.94, -18.29, 2.69], # N
                 [-0.89,  -3.06, 0.86], # h N
                 [-0.90,   9.62, 0.71]] # norm
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
    shape = "shell"
    h = $height
    R = $radius # not used

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
    h, R, density = theta_to_params(theta)

    if (1e1 < gamma_min < gamma_break < gamma_max < 1e10 and

        1     < p1  < p2 < 6    and
        1e10  < h   < R  < 1e26 and
        0.005 < B        < 0.5  and
        2     < Gamma    < 192  and
        1e-26 < density  < 1e-16):
           return 0
    else:
        return -np.inf

def theta_to_config(theta):
    gamma_min, gamma_break, gamma_max,  \
    p1, p2,                             \
    B, Gamma,                           \
    h, R, density = theta_to_params(theta)

    config = config_template.substitute(
                g_min = gamma_min,
                g_break = gamma_break,
                g_max = gamma_max,
                p1 = p1,
                p2 = p2,
                mag = B,
                gamma = Gamma,
                height = h,
                radius = R,
                dens = density,
                # set dt_max to a tenth of the escape time (more or less)
                dt_max = h / 3e10 / 10)

    return config

def theta_to_params(theta):
    log_nuFnu_opt, \
    log_xray_opt,  \
    log_gray_opt,  \
    log_B,         \
    log_Gamma_B,   \
    log_nu_break,  \
    log_rho,       \
    log_N,         \
    log_h_N,       \
    log_norm = theta

    log_Gamma = log_Gamma_B - log_B
    log_h     = log_h_N     - log_N

    log_g_break = (log_nu_break - log_Gamma_B) / 2

    aux = log_nuFnu_opt - 2 * log_nu_break
    
    log_V = aux - log_gray_opt - log_Gamma_B + log_h + log_g_break + log_h_N/2
    # log_V = aux - log_gray_opt - log_Gamma_B + log_h + log_g_break
    log_R = (log_V - log_h) / 2

    p1 = -(aux - log_N - log_V) / log_g_break - 1
    p2 = p1 + log_norm / log_g_break

    log_g_max = log_xray_opt / (3 - p2) + log_g_break

    log_g_min = (log_N - log_rho - np.log10(p1 - 1)) / (p1 - 1)

    return np.array([10**log_g_min,
                     10**log_g_break,
                     10**log_g_max,
                     p1,
                     p2,
                     10**log_B,
                     10**log_Gamma,
                     10**log_h,
                     10**log_R,
                     10**log_rho])

def params_to_theta(params):
    gamma_min, gamma_break, gamma_max,  \
    p1, p2,          \
    B, Gamma,                           \
    h, R, density = params

    N = density * (p1 - 1) * np.power(gamma_min, p1 - 1)
    nu_br = Gamma * B * gamma_break**2
    V = h * R**2

    flux_opt = N * V * nu_br**2 / np.power(gamma_break, p1 + 1)
    flux_x_ray_opt = np.power(gamma_max / gamma_break, 3 - p2)
    flux_g_ray_opt = np.power(h * N, 1.5) * np.power(gamma_break, -p1) / (Gamma * B)
    # flux_g_ray_opt = h * N * np.power(gamma_break, -p1) / (Gamma * B)

    norm = np.power(gamma_break, p2 - p1)

    return np.array([np.log10(flux_opt),
                     np.log10(flux_x_ray_opt),
                     np.log10(flux_g_ray_opt),
                     np.log10(B),
                     np.log10(Gamma * B),
                     np.log10(nu_br),
                     np.log10(density),
                     np.log10(N),
                     np.log10(h * N),
                     np.log10(norm)])

def get_gamma(theta):
    log_nuFnu_opt, \
    log_xray_opt,  \
    log_gray_opt,  \
    log_B,         \
    log_Gamma_B,   \
    log_nu_break,  \
    log_rho,       \
    log_N,         \
    log_h_N,       \
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
    log_N,         \
    log_h_N,       \
    log_norm = theta

    log_h     = log_h_N     - log_N

    log_g_break = (log_nu_break - log_Gamma_B) / 2

    aux = log_nuFnu_opt - 2 * log_nu_break
    
    log_V = aux - log_gray_opt - log_Gamma_B + log_h + log_g_break + log_h_N/2
    log_R = (log_V - log_h) / 2

    return 10**(log_R)

def assert_consistency(theta):
    params = theta_to_params(theta)
    # Assert that theta_to_params and params_to_theta are actually inverse functions
    np.testing.assert_allclose(params_to_theta(params), theta)

    # Assert that you are not forgetting R or Gamma
    np.testing.assert_allclose(params[6], get_gamma(theta))
    np.testing.assert_allclose(params[8], get_R(theta))

    return True
