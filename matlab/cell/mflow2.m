function [flow_gas,flow_water,flow_H2,flow_N2]=mflow2(y_H2,p_water,p_gas,flow_tot) 
    M_H2=2*10^-3;
    M_N2=28*10^-3;
    y_N2=1-y_H2;
    M_gas=y_H2*M_H2+y_N2*M_N2;
    M_water=18*10^-3;
    
    omega=M_water*p_water/M_gas/p_gas;
    flow_gas=flow_tot/(1+omega);
    flow_water=flow_tot-flow_gas;
    
    x_H2=y_H2*M_H2/M_gas;
    x_N2=y_N2*M_N2/M_gas;
    flow_H2=flow_gas*x_H2;
    flow_N2=flow_gas*x_N2;
end