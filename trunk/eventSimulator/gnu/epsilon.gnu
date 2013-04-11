mpi = 0.139
mka = 0.493
mpr = 0.938
mel = 0.511e-3

beta(bg) = bg/sqrt(bg*bg + 1)
gamm(bg) =    sqrt(bg*bg + 1)

Z_A = 0.49848
I   = 173.e-6     # MeV
rho = 2.33        # g cm^{-2}

depth = 450e-4      # cm
K     = 0.307075    # MeV g^-1 cm^2
me    = 0.510998918 # MeV

# Load Sternheimer84
load "../../siEnergyLoss/data/density_Sternheimer84.par"

# Density correction
delta(x) = (x > x0 ? 2*log(10)*x - C + (x < x1 ? a*(x1 - x)**k : 0) : \
            d0(x0) * 10**(2 * (x-x0)) ) 
d0(x0) = delta(x0 + 1e-6)

# Xi
xi(bg) = K/2 * Z_A * (depth * rho)/beta(bg)**2

# MeV/cm
epsilon(bg) = (xi(bg) / depth * \
  (log(2*me * bg**2 * xi(bg)/I**2) + 0.200 - beta(bg)**2 - delta(log10(bg))))
