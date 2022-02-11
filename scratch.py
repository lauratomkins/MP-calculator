import numpy as np

drop_radius_um = np.array([30, 50, 100, 150, 200, 300])
drop_radius_cm = drop_radius_um * (1e-4)

drop_mass_g = (3.8e-3) * (drop_radius_cm**2)
drop_mass_ug = drop_mass_g * (1e6)
drop_mass_mg = drop_mass_g * (1e3)