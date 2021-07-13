function [flow_air,flow_water,flow_O2,flow_N2]=mflow(y_O2,p_water,p_air,flow_tot) 
    M_O2=32*10^-3;
    M_N2=28*10^-3;
    y_N2=1-y_O2;
    M_air=y_O2*M_O2+y_N2*M_N2;
    M_water=18*10^-3;

    omega=M_water*p_water/M_air/p_air;
    flow_air=flow_tot/(1+omega);
    flow_water=flow_tot-flow_air;
    
    x_O2=y_O2*M_O2/M_air;
    x_N2=y_N2*M_N2/M_air;
    flow_O2=flow_air*x_O2;
    flow_N2=flow_air*x_N2;
end