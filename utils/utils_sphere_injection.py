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
        "$\\log_{10}\\left(\\rho \cdot h \cdot \\sqrt[3]{\\frac{\\gamma_{min}^{2}} {\\Gamma \\cdot B}} \\right)$"]

initial_theta = [[-2.55,  2.46, 2.86],
                 [-0.13,  8.55, 0.15],
                 [-0.02, 49.36, 0.02],
                 [-1.42,  4.16, 1.37],
                 [-0.52, -8.89, 0.67],
                 [-0.64, 10.30, 0.49],
                 [-0.56, -1.04, 0.53],
                 [-0.24, -0.15, 0.24],
                 [-0.07, 41.35, 0.07]]
initial_theta = np.array(initial_theta).T

config_template = Template("""
    [general]
    density = 1e-25
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

    [external_injection]
        luminosity = $L
        eta        = 1

        [external_injection.electrons]
        distribution_type = "broken_power_law"
        break_point  = $g_break
        first_slope  = $p1
        second_slope = $p2

        [external_injection.protons]
        distribution_type = "power_law"
        slope = 2
    """)

def lnprior(theta):
    gamma_min, gamma_break, gamma_max,  \
    p1, p2,                             \
    B, Gamma,                           \
    R,                                  \
    luminosity = theta_to_params(theta)

    if (1e1 < gamma_min < gamma_break < gamma_max < 1e10 and

        1     < p1 < p2 < 6    and
        0.005 < B       < 0.5  and
        2     < Gamma   < 192  and
        1e10  < R       < 1e26 and

        1e38 < luminosity < 1e50):
           return 0
    else:
        return -np.inf

def theta_to_config(theta):
    gamma_min, gamma_break, gamma_max,  \
    p1, p2,                             \
    B, Gamma,                           \
    R,                                  \
    luminosity = theta_to_params(theta)

    config = config_template.substitute(
                g_min = gamma_min,
                g_break = gamma_break,
                g_max = gamma_max,
                p1 = p1,
                p2 = p2,
                mag = B,
                gamma = Gamma,
                radius = R,
                L = luminosity,
                # set dt_max to a tenth of the escape time (more or less)
                dt_max = h / 3e10 / 10)

    return config

def theta_to_params(theta):
    log_N,          \
    log_nu_br,      \
    log_nuFnu_opt,  \
    log_opt_xray,   \
    log_opt_gamma,  \
    log_norm,       \
    log_B,          \
    log_Gamma_B,    \
    log_L = theta

    log_g_break = (log_nu_br - log_Gamma_B) / 2
    p1 = 4  + log_opt_gamma / log_g_break
    p2 = p1 + log_norm      / log_g_break

    log_g_max = log_g_break - log_opt_xray / (3 - p2)

    log_V = log_nuFnu_opt - log_N - 2 * log_nu_br + (1 + p1) * log_g_break

    log_R = log_V / 3

    log_g_min = (log_N - np.log10(p1 - 1)) / (p1 - 1)

    return np.array([10**log_g_min,
                     10**log_g_break,
                     10**log_g_max,
                     p1,
                     p2,
                     10**log_B,
                     10**(log_Gamma_B - log_B),
                     10**(log_R),
                     10**(log_L)])

def params_to_theta(params):
    gamma_min, gamma_break, gamma_max,  \
    p1, p2,                             \
    B, Gamma,                           \
    R,                                  \
    luminosity = params

    N = (p1 - 1) * np.power(gamma_min, p1 - 1)
    nu_br = Gamma * B * gamma_break**2
    V = R**3

    return np.array([np.log10(N),
                     np.log10(nu_br),
                     np.log10(N * V * nu_br**2 / np.power(gamma_break, p1 + 1)),
                     np.log10(np.power(gamma_break / gamma_max, 3 - p2)),
                     np.log10(np.power(gamma_break, p1 - 4)),
                     np.log10(np.power(gamma_break, p2 - p1)),
                     np.log10(B),
                     np.log10(Gamma * B),
                     np.log10(luminosity)])

def get_gamma(theta):
    log_N,          \
    log_nu_br,      \
    log_nuFnu_opt,  \
    log_opt_xray,   \
    log_opt_gamma,  \
    log_norm,       \
    log_B,          \
    log_Gamma_B,    \
    log_L = theta

    return 10**(log_Gamma_B - log_B)

def get_R(theta):
    log_N,          \
    log_nu_br,      \
    log_nuFnu_opt,  \
    log_opt_xray,   \
    log_opt_gamma,  \
    log_norm,       \
    log_B,          \
    log_Gamma_B,    \
    log_L = theta

    log_g_break = (log_nu_br - log_Gamma_B) / 2
    p1 = 4  + log_opt_gamma / log_g_break

    log_V = log_nuFnu_opt - log_N - 2 * log_nu_br + (1 + p1) * log_g_break

    log_R = log_V / 3

    return 10**(log_R)

def assert_consistency(theta):
    params = theta_to_params(theta)
    # Assert that theta_to_params and params_to_theta are actually inverse functions
    np.testing.assert_allclose(params_to_theta(params), theta)

    # Assert that you are not forgetting R or Gamma
    np.testing.assert_allclose(params[6], get_gamma(theta))
    np.testing.assert_allclose(params[7], get_R(theta))

    return True
