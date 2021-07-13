function [flow]=massflow(y_O2,p_water,p_gas,flow_tot,cv) 
%% 工质摩尔质量
R=8.3145;
M_H2=2*10^-3;    
M_O2=32*10^-3;
M_N2=28*10^-3;
y_N2=1-y_O2;
M_air=y_O2*M_O2+y_N2*M_N2;
M_water=18*10^-3;
if cv=='ca'
    M_gas=M_air;
else if cv=='an'
        M_gas=M_H2;
    else
        print('cv error');
    end
end
%% 气体和水蒸气的质量流量，阴极气体中含氧气和氮气，阳极气体仅含氢气
omega=M_water*p_water/M_gas/p_gas;
flow_gas=flow_tot/(1+omega);
flow_water=flow_tot-flow_gas;
if cv=='ca'
    x_O2=y_O2*M_O2/M_gas;
    x_N2=y_N2*M_N2/M_gas;
    flow_O2=flow_gas*x_O2;
    flow_N2=flow_gas*x_N2;
    flow=[flow_water,flow_gas,flow_O2,flow_N2];
else if cv=='an'
        flow=[flow_water,flow_gas];
    else
        print('flow error'); 
    end
end    


end