"""
The 'main.jl' is written in Geometrized unit system, where distances are expressed in centimeters [cm].
If energy_density and pressure are given in [g/cm^3] unit, 
the inputs shoulds be energy_density/unit_g and pressure/unit_g [1/cm^3].
On the other hand, the 'piecewise_polytrope_eos_cgs.jl' is written in CGS unit.
"""
# Physical constants
const c = 2.99792458e10  # [cm/sec]
const G = 6.67428e-8  # [cm^3 /(g sec^2)]
const Msun = 1.9884e33  # [g]
const hbarc = 197.32698044404107 # [MeV fm]

# normalization constants. Geometrized unit system!
const unit_g = c^2/G
const unit_l = 1.0
const eps_ref = c^4/G

# unit conversion factor
const dyncm2_to_MeVfm3 = 1.0/(1.602176634e33)
const MeVfm3_to_dyncm2 = 1.602176634e33
const gcm3_to_MeVfm3 = 1.0/(1.7826619216278976e12)
const MeVfm3_to_gcm3 = 1.7826619216278976e12  # 1e(9+39)*1.602e-19/(3e8)**2  正しくは MeV/fm3/c2 to g/cm3
const gcm3_to_Msun2 = G/c^2*(G*Msun/c^2)^2 