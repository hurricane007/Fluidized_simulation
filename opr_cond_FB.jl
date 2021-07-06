if !@isdefined g
    const g = 9.81
end
if !@isdefined Rg
    const Rg = 8.3145
end

# bed properties
Lbed = 0.04
bedDiam = 0.032 # bed diameter
Ω_r = (bedDiam^2*pi*0.25)
# particle properties
dp = 250E-6 # [m]
ρs = 1800 # kg/m³
# Operation conditions
εmf = 0.55 # voidage at minimum fluidization
pres = 1.1E5 # pressure, [Pa]
y_CO2 = 0.145 # inlet CO2 fraction
y_CO = 0.03 # 3%
y_NO = 500E-6
Temp = 650. + 273.15 # [K]
Qtot = 2.2E3 # feed gas flow rate, [Nml/min]
yA_in = [y_CO2, y_CO, y_NO, (1-y_CO2-y_CO-y_NO)]
usg0 = Qtot/60*1E-6*Temp/273.15/Ω_r # [m^3/m^2*s] superficial gas velocity
Ctot = pres/Rg/Temp
CA0 = Ctot*yA_in
ρg0 = Ctot*sum(yA_in.*[0.044, 0.028, 0.030, 0.028])
μm0 = 3.94E-5 # Viscosity of Air
Deg = 1E-4 # [m^2/s]
Des = 0.06+0.1*usg0

# Calculate initial conditions
umf0 = 1.118E-13*(dp*1E6)^1.82*(ρs - ρg0)^0.94/(ρg0^0.06*μm0^0.88)
fbub0 = usg0 - umf0
db0 = 0.54/g^0.2*(usg0-umf0)^0.4*(4*sqrt(0.25pi*0.002^2))^0.8
ubr0 = 0.711*(g*db0)^0.5
ub0 = usg0-umf0+ubr0
fb0 = (usg0-umf0)/ub0
fes0 = (1-fb0)*(1-εmf)
FAin = [fbub0*yA_in*Ctot; umf0*yA_in*Ctot]
