import numpy as np

def obs_to_data(obs_data):
    a = (obs_data[0] - 150) / (592 - 150)
    b = (obs_data[1] - 113) / (725 - 113)

    c = 10**(-6  * (1 - a) + a * 3)
    d = 10**(-10 * (1 - b) + b * -14)

    return np.array([c,d])

# Inlined experimental data
radio_obs = np.array([
    [200, 708],
    [207, 683],
    [216, 652],
    [224, 622],
    [227, 606],
    [232, 591],
    [232, 586],
    [238, 573]
]).transpose()
radio_obs_data = obs_to_data(radio_obs)

optical_obs = np.array([
    [460, 185],
    [459, 185],
    [463, 170],
    [463, 177],
    [467, 188],
    [467, 179],
    [472, 177],
    [479, 185],
    [482, 177],
    [484, 182]
]).transpose()
optical_obs_data = obs_to_data(optical_obs)
# Assume that the errors in optical data are small
optical_obs_errors = 0.1 * optical_obs_data[1]

x_rays_obs = np.array([
    [573, 350],
    [582, 347],
    [588, 361],
    [593, 382],
    [597, 382],
    [601, 391],
    [604, 396],
    [610, 414],
    [618, 420],
    [625, 430],
    [632, 429],
    [636, 415],
    [643, 412],
    [650, 398],
    [663, 400],
    [679, 380]
]).transpose()
x_rays_obs_data = obs_to_data(x_rays_obs)
# Assume that the errors in x-ray data are small
x_rays_obs_errors = 0.2 * x_rays_obs_data[1]

gamma_rays_obs = np.array([
    [ 851, 148],
    [ 854, 121],
    [ 875, 152],
    [ 900, 153],
    [ 924, 150],
    [ 948, 164],
    [ 973, 200],
    [ 981, 286],
    [ 998, 334],
    [1015, 399]
]).transpose()
gamma_rays_obs_data = obs_to_data(gamma_rays_obs)
# Assume that the upper errors in gamma ray data are small
# but the lower errors are awful
gamma_rays_obs_errors = [0.5 * gamma_rays_obs_data[1], (0.2 + 0.5) * gamma_rays_obs_data[1]]
