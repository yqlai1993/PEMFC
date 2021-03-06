clear all;
clc;


%% 定值    
M_H2=2*10^-3;
M_O2=32*10^-3;
M_N2=28*10^-3;
M_water=18*10^-3;
R=8.3145;
R_O2=R/M_O2;
R_N2=R/M_N2;
R_H2=R/M_H2;
R_water=R/M_water;
F=96485;

thick_mem=1.275*10^-4;
area_mem=232*10^-4;

T_amb=273.15+23.5;
p_amb=101000;
h_amb=17;

volu_ca=0.01;
k_ca_in=0.36*10^-3;
k_ca_out=0.22*10^-5;
h_ca=10;
volu_an=0.005;
k_an_in=0.36*10^-5;
k_an_out=0.22*10^-5;
h_an=2;
dH_react=1.196*10^8;

T_cl=273.15+25;
h_cl=0;

Cp_stack=350;
m_cell=10/120;
num_cell=120;
current=10;

T_ca_in=273.15+25;
T_an_in=273.15+25;
p_ca_in=101000*5;
p_an_in=101000*7;
p_ca_out=101000*1;
p_an_out=101000*1;
Rhum_ca_in=0.5;
Rhum_an_in=0.5;

%% 初值
T_ca=273.15+25;
T_an=273.15+25;
T_cell=273.15+25;
p_ca=101000*1;
p_an=101000*1;
Rhum_ca=0.0;
Rhum_an=0.0;
[m_water_ca,m_O2_ca,m_N2_ca]=ini(0.21,p_ca,volu_ca,T_ca,Rhum_ca,'ca');
[m_water_an,m_H2_an,m_N2_an]=ini(1,p_an,volu_an,T_an,Rhum_an,'an');
tspan=100;

for i=1:tspan
    %% 
    % 阴极入口压力
    p_wv_ca_in=Rhum_ca_in*psat(T_ca_in);
    p_air_ca_in=p_ca_in-p_wv_ca_in; 
    [rho_wv_ca_in,cv_wv_ca_in,cp_wv_ca_in]=PropWV(p_wv_ca_in,T_ca_in);
    [rho_wl_ca_in,cv_wl_ca_in,k_wl_ca_in,mu_wl_ca_in]=PropWL(T_ca_in);
    % 阴极内部压力
    p_sat_ca=psat(T_ca);
    [m_wv_ca,m_wl_ca]=mvap(m_water_ca,p_sat_ca,volu_ca,R_water,T_ca);
    p_O2_ca=stateEq(m_O2_ca,R_O2,T_ca,volu_ca);
    p_N2_ca=stateEq(m_N2_ca,R_N2,T_ca,volu_ca);
    p_wv_ca=stateEq(m_wv_ca,R_water,T_ca,volu_ca) ;
    p_air_ca=p_O2_ca+p_N2_ca;
    p_ca=p_air_ca+p_wv_ca;
    Rhum_ca=p_wv_ca/p_sat_ca;
    [rho_O2_ca,cv_O2_ca,cp_O2_ca]=PropO2(p_O2_ca,T_ca);
    [rho_N2_ca,cv_N2_ca,cp_N2_ca]=PropN2(p_N2_ca,T_ca);
    [rho_wv_ca,cv_wv_ca,cp_wv_ca]=PropWV(p_wv_ca,T_ca);    
    [rho_wl_ca,cv_wl_ca,k_wl_ca,mu_wl_ca]=PropWL(T_ca);
    p_c(i)=p_ca;
    % 阴极入口流量
    flow_ca_in=k_ca_in*(p_ca_in-p_ca)/p_amb;
    y_O2_ca_in=0.21;
    [flow_air_ca_in,flow_wv_ca_in,flow_O2_ca_in,flow_N2_ca_in]=massflow(y_O2_ca_in,p_wv_ca_in,p_air_ca_in,flow_ca_in,'ca');
    flow_O2_ci(i)=flow_O2_ca_in;
    flow_w_ci(i)=flow_wv_ca_in;
    [rho_O2_ca_in,cv_O2_ca_in,cp_O2_ca_in]=PropO2(p_air_ca_in*flow_O2_ca_in/0.002/(flow_O2_ca_in/0.002+flow_N2_ca_in/0.028),298.15);
    [rho_N2_ca_in,cv_N2_ca_in,cp_N2_ca_in]=PropN2(p_air_ca_in*flow_N2_ca_in/0.028/(flow_O2_ca_in/0.002+flow_N2_ca_in/0.028),298.15);
    %阴极出口流量
    flow_ca_out=k_ca_out*(p_ca-p_ca_out)/p_amb;
    [flow_air_ca_out,flow_wv_ca_out,flow_O2_ca_out,flow_N2_ca_out]=massflow(p_O2_ca/p_air_ca,p_wv_ca,p_air_ca,flow_ca_out,'ca');
    flow_w_co(i)=flow_wv_ca_out;
    flow_O2_co(i)=flow_O2_ca_out;
    flow_N2_co(i)=flow_N2_ca_out;
    
    %% 
    % 阳极入口压力
    p_wv_an_in=Rhum_an_in*psat(T_an_in);
    p_mix_an_in=p_an_in-p_wv_an_in; 
    [rho_wv_an_in,cv_wv_an_in,cp_wv_an_in]=PropWV(p_wv_an_in,T_an_in);
    [rho_wl_an_in,cv_wl_an_in,k_wl_an_in,mu_wl_an_in]=PropWL(T_an_in);
    % 阳极内部压力
    p_sat_an=psat(T_an);
    [m_wv_an,m_wl_an]=mvap(m_water_an,p_sat_an,volu_an,R_water,T_an);
    p_H2_an=stateEq(m_H2_an,R_H2,T_an,volu_an);
    p_N2_an=stateEq(m_N2_an,R_N2,T_an,volu_an);
    p_wv_an=stateEq(m_wv_an,R_water,T_an,volu_an);
    p_mix_an=p_H2_an+p_N2_an;
    p_an=p_mix_an+p_wv_an;
    Rhum_an=p_wv_an/p_sat_an;
    [rho_H2_an,cv_H2_an,cp_H2_an]=PropH2(p_H2_an,T_an);
    [rho_N2_an,cv_N2_an,cp_N2_an]=PropN2(p_N2_an,T_an);
    [rho_wv_an,cv_wv_an,cp_wv_an]=PropWV(p_wv_an,T_an);
    [rho_wl_an,cv_wl_an,k_wl_an,mu_wl_an]=PropWL(T_an);
    p_a(i)=p_an;
    % 阳极入口流量
    flow_an_in=k_an_in*(p_an_in-p_an)/p_amb;
    [flow_mix_an_in,flow_wv_an_in,flow_H2_an_in,flow_N2_an_in]=massflow(1,p_wv_an_in,p_mix_an_in,flow_an_in,'an');
    flow_H2_ai(i)=flow_H2_an_in;
    flow_w_ai(i)=flow_wv_an_in;
    [rho_H2_an_in,cv_H2_an_in,cp_H2_an_in]=PropH2(p_mix_an_in*flow_H2_an_in/0.002/(flow_H2_an_in/0.002+flow_N2_an_in/0.028),298.15);
    [rho_N2_an_in,cv_N2_an_in,cp_N2_an_in]=PropN2(p_mix_an_in*flow_N2_an_in/0.028/(flow_H2_an_in/0.002+flow_N2_an_in/0.028),298.15);
    % 阳极出口流量
    flow_an_out=k_an_out*(p_an-p_an_out)/p_amb;
    [flow_mix_an_out,flow_wv_an_out,flow_H2_an_out,flow_N2_an_out]=massflow(p_H2_an/p_mix_an,p_wv_an,p_mix_an,flow_an_out,'an');
    flow_w_ao(i)=flow_wv_an_out;
    flow_H2_ao(i)=flow_H2_an_out;
    flow_N2_ao(i)=flow_N2_an_out;
    
    %% 反应量
    flow_H2_react=num_cell*M_H2*current/2/F;
    flow_O2_react=num_cell*M_O2*current/4/F;
    flow_water_gen=num_cell*M_water*current/2/F;

    %% 膜水合流量
    [sigma_ca,c_ca]=memW(Rhum_ca);
    [sigma_an,c_an]=memW(Rhum_an);
    sigma_ave=(sigma_an+sigma_ca)/2;
    n_d=0.0029*sigma_ave^2+0.05*sigma_ave-3.4e-19;
    D_w=1.25e-6*exp(2416*(1/303-1/T_cell));
    currentDen=current/area_mem; 
    N_water_mem=n_d*currentDen/F-D_w*(c_ca-c_an)/thick_mem; 
    flow_water_mem=num_cell*N_water_mem*M_water*area_mem/100000;
    flow_w_m(i)=flow_water_mem;
   
    %% 电压
    volt_nst=1.229-8.5e-4*(T_cell-298.15)+R*T_cell*(log(p_H2_an)+0.5*log(p_O2_ca))/(2*F);
    volt_act=0.9514-3.12e-3*T_cell+1.87e-4*T_cell*log(current)-7.4e-5*T_cell*log((p_O2_ca/p_amb)/(5.08e6*exp(-498/T_cell)));
%     resistivity_mem=(thick_mem/area_mem)*181.6*(1+0.03*currentDen+0.062*(T_cell/303)^2*currentDen^2.5)/((sigma_an-0.634-3*currentDen)*exp(4.18*(1-303/T_cell)));
    resistivity_mem=0.01605-3.5e-5*T_cell+8e-5*current;
    volt_ohm=current*resistivity_mem;
    volt_cell=volt_nst-volt_act-volt_ohm;
    power_cell=volt_cell*current;
%     eff_cell=(flow_H2_react/flow_H2_an_in)*(volt_cell/1.48);
    
    %% 控制方程
    m_O2_ca=m_O2_ca+(flow_O2_ca_in-flow_O2_ca_out-flow_O2_react);
    m_N2_ca=m_N2_ca+(flow_N2_ca_in-flow_N2_ca_out);
    m_water_ca=m_water_ca+(flow_wv_ca_in-flow_wv_ca_out+flow_water_gen+flow_water_mem);
 
    m_H2_an=m_H2_an+(flow_H2_an_in-flow_H2_an_out-flow_H2_react);
    m_N2_an=m_N2_an+(flow_N2_an_in-flow_N2_an_out);
    m_water_an=m_water_an+(flow_wv_an_in-flow_wv_an_out-flow_water_mem);
    
%     T_ca=T_ca+(h_ca*(T_cell-T_ca)+(cp_O2_ca_in*flow_O2_ca_in+cp_N2_ca_in*flow_N2_ca_in+cp_wv_ca_in*flow_wv_ca_in)*T_ca_in-(cp_O2_ca*(flow_O2_ca_out+flow_O2_react)+cp_N2_ca*flow_N2_ca_out+cp_wv_ca*flow_wv_ca_out)*T_ca)/(cp_O2_ca*m_O2_ca+cp_N2_ca*m_N2_ca+cp_wv_ca*m_wv_ca+cv_wl_ca*m_wl_ca);    
%     T_an=T_an+(h_an*(T_cell-T_an)+(cp_H2_an_in*flow_H2_an_in+cp_N2_an_in*flow_N2_an_in+cp_wv_an_in*flow_wv_an_in)*T_an_in-(cp_H2_an*(flow_H2_an_out+flow_H2_react)+cp_N2_an*flow_N2_an_out+cp_wv_an*flow_wv_an_out)*T_an)/(cp_H2_an*m_H2_an+cp_N2_an*m_N2_an+cp_wv_an*m_wv_an+cv_wl_an*m_wl_an);
%     T_cell=T_cell+(h_ca*(T_ca-T_cell)+h_an*(T_an-T_cell)+h_amb*(T_amb-T_cell)+h_cl*(T_cl-T_cell)+flow_H2_react*dH_react-num_cell*volt_cell*current)/(Cp_stack*m_cell*num_cell);  
    
    T_cell=T_cell+(h_amb*(T_amb-T_cell)+h_cl*(T_cl-T_cell)+flow_H2_react*dH_react-num_cell*volt_cell*current+(cp_O2_ca_in*flow_O2_ca_in+cp_N2_ca_in*flow_N2_ca_in+cp_wv_ca_in*flow_wv_ca_in)*T_ca_in-(cp_O2_ca*(flow_O2_ca_out+flow_O2_react)+cp_N2_ca*flow_N2_ca_out+cp_wv_ca*flow_wv_ca_out)*T_ca+(cp_H2_an_in*flow_H2_an_in+cp_N2_an_in*flow_N2_an_in+cp_wv_an_in*flow_wv_an_in)*T_an_in-(cp_H2_an*(flow_H2_an_out+flow_H2_react)+cp_N2_an*flow_N2_an_out+cp_wv_an*flow_wv_an_out)*T_an)/(Cp_stack*m_cell*num_cell+cp_O2_ca*m_O2_ca+cp_N2_ca*m_N2_ca+cp_wv_ca*m_wv_ca+cv_wl_ca*m_wl_ca+cp_H2_an*m_H2_an+cp_N2_an*m_N2_an+cp_wv_an*m_wv_an+cv_wl_an*m_wl_an);
    T_ca=T_cell;
    T_an=T_cell;
    
    m_O2_c(i)=m_O2_ca;
    m_N2_c(i)=m_N2_ca;
    m_w_c(i)=m_wv_ca;
    m_H2_a(i)=m_H2_an;
    m_N2_a(i)=m_N2_an;
    m_w_a(i)=m_wv_an;
    T_c(i)=T_ca;
    T_a(i)=T_an;
    T_s(i)=T_cell;
    volt_c(i)=volt_cell;
end

t=[1:tspan];
plot(t,flow_O2_co,'r',t,flow_H2_ao,'g',t,flow_N2_co,'m',t,flow_N2_ao,'b');
figure;
plot(t,T_c,'r',t,T_a,'g',t,T_s,'b');
figure;
plot(t,volt_c,'k');