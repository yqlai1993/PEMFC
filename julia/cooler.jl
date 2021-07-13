using ModelingToolkit,OrdinaryDiffEq

function psat(T)
  T=T-273.15
  # a=-2.1794+0.02953*T-9.1837*10^-5*T^2+1.4454*10^-7*T^3
  a=-2.18+0.0295*T-9.18*10^-5*T^2+1.44*10^-7*T^3
  p=10^(a+5)
end

function mflow(x_O2,p_water,p_air,flow_tot)
  x_N2=1-x_O2
    
  M_O2=32*10^-3
  M_N2=28*10^-3
  M_water=18*10^-3
  M_air=x_O2*M_O2+x_N2*M_N2

  omega=M_water*p_water/M_air/p_air
  y_O2=x_O2*M_O2/M_air
  y_N2=x_N2*M_N2/M_air

  flow_air=flow_tot/(1+omega)
  flow_water=flow_tot-flow_air
  flow_O2=flow_air*y_O2
  flow_N2=flow_air*y_N2
  return F=[flow_air,flow_water,flow_O2,flow_N2]
end

function cooler(p_in,flow_tot)
  x_O2=0.21
  T_amb=298.15
  T_out=273.15+80 
  p_amb=101000
  phi_amb=0.5
  phi_out=p_in*phi_amb*psat(T_amb)/(p_amb*psat(T_out))
  p_water=phi*psat(T_out)
  p_air=p_in-p_water
  (flow_air,flow_water,flow_O2,flow_N2)=mflow(x_O2,p_water,p_air,flow_tot)
end