def lnprior(theta):
    if 1  < gamma_min    < gamma_break  < gamma_max < 1e8 and \
       1  < first_slope  < second_slope < 6 and               \
                                                              \
       1e25  <  Sf_R2   < 1e35 and  \
       1e12  <  h       < 1e18 and  \
       0.05  <  B       < 0.3 and   \
       5     <  Gamma   < 30 and    \
       1e-24 <  density < 1e-18:
           return 0
    else:
        return -np.inf

limits = np.array([(1, 7),          # log10 g_min
                   (1, 7),          # log10 g_break
                   (1, 7),          # log10 g_max
                   (1, 5),          # p1
                   (1, 5),          # p2
                   (25, 35),        # log10 Sf·R^2
                   (12, 18),        # log10 h
                   (0.05, 0.30),    # B
                   (5, 30),         # Gamma
                   (-10, +2)])      # log10 (d·h·g_min^2/3)

labels=["$\\log_{10}\\left(\\gamma_{min}\\right)$", \
        "$\\log_{10}\\left(\\gamma_{break}\\right)$", \
        "$\\log_{10}\\left(\\gamma_{max}\\right)$", \
        "$p_1$", \
        "$p_2$", \
        "$\\log_{10}\\left(S_f \cdot R^2\\right)$", \
        "$\\log_{10}\\left(h\\right)$", \
        "$B$", \
        "$\\Gamma$", \
        "$\\log_{10}\\left(\\rho \cdot h \cdot \\gamma_{min}^{2/3}\\right)$"]
