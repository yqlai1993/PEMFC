function [flow_gas,flow_water,flow_g1,flow_g2]=massflow(y1,p_wv,p_gas,flow_tot,cv) 
%% 工质摩尔质量
M_H2=2*10^-3;    
M_O2=32*10^-3;
M_N2=28*10^-3;
M_water=18*10^-3;

if cv=='ca'
    M1=M_O2;
    M2=M_N2;
else if cv=='an'
        M1=M_H2;
        M2=M_N2;
    else
        print('cv error');
    end
end

y2=1-y1;
M_gas=y1*M1+y2*M2;
omega=M_water*p_wv/M_gas/p_gas;
flow_gas=flow_tot/(1+omega);
flow_water=flow_tot-flow_gas;
flow_g1=flow_gas*y1*M1/M_gas;
flow_g2=flow_gas*y2*M2/M_gas;

end    

