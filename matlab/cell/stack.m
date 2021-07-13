function [ dy ] = stack( t,y )
%UNTITLED10 此处显示有关此函数的摘要
%   此处显示详细说明

    m_O2_ca=y(1);
    m_N2_ca=y(2);
    m_water_ca=y(3);
    m_H2_an=y(4);
    m_N2_an=y(5);
    m_water_an=y(6);
    T_ca=y(7);
    T_an=y(8);
    T_cell=y(9);
    
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
    y_O2_amb=0.21;
    y_N2_amb=1-y_O2_amb;
    M_air_amb=y_O2_amb*M_O2+y_N2_amb*M_N2;
    h_amb=17;

    volu_ca=0.01;
    k_ca_in=0.36*10^-5;
    k_ca_out=1.4*10^-5;
    h_ca=10;
    volu_an=0.005;
    k_an_in=0.36*10^-5;
    k_an_out=0.22*10^-5;
    h_an=2;
    dH_react=1.196*10^8;

    T_cl=273.15+23.5;
    h_cl=50;

    Cp_stack=35;
    m_cell=45/35;
    num_cell=35;
    current=20;
    
    T_ca_in=273.15+23.5;
    T_an_in=273.15+23.5;
    p_ca_in=114000*3;
    p_an_in=114000*3;
    p_ca_out=101000;
    p_an_out=101000;
    Rhum_ca_in=0.0;
    Rhum_an_in=0.0;

    %% 阴极入口压力
    p_sat_ca_in=psat(T_ca_in);
    p_water_ca_in=Rhum_ca_in*p_sat_ca_in;
    p_air_ca_in=p_ca_in-p_water_ca_in;
    
%     [rho_air_ca_in,cv_air_ca_in,cp_air_ca_in,k_air_ca_in,mu_air_ca_in]=PropAir(p_air_ca_in,T_ca_in);
    [rho_air_ca_in,cv_air_ca_in,cp_air_ca_in]=PropO2(p_air_ca_in,T_ca_in);%%paper
    [rho_wv_ca_in,cv_wv_ca_in,cp_wv_ca_in]=PropWV(p_water_ca_in,T_ca_in);
    [rho_wl_ca_in,cv_wl_ca_in,k_wl_ca_in,mu_wl_ca_in]=PropWL(T_ca_in);
    
    %% 阴极内部压力
    p_sat_ca=psat(T_ca);
    [m_wv_ca,m_wl_ca]=mvap(m_water_ca,p_sat_ca,volu_ca,R_water,T_ca);
    p_O2_ca=stateEq(m_O2_ca,R_O2,T_ca,volu_ca);
    p_N2_ca=stateEq(m_N2_ca,R_N2,T_ca,volu_ca);
    p_water_ca=stateEq(m_wv_ca,R_water,T_ca,volu_ca) ;
    p_air_ca=p_O2_ca+p_N2_ca;
    p_ca=p_air_ca+p_water_ca;
    Rhum_ca=p_water_ca/p_sat_ca;
    
    [rho_O2_ca,cv_O2_ca,cp_O2_ca]=PropO2(p_O2_ca,T_ca);
    [rho_N2_ca,cv_N2_ca,cp_N2_ca]=PropN2(p_N2_ca,T_ca);
    [rho_wv_ca,cv_wv_ca,cp_wv_ca]=PropWV(p_water_ca,T_ca);    
    [rho_wl_ca,cv_wl_ca,k_wl_ca,mu_wl_ca]=PropWL(T_ca);
    
    %% 阴极入口流量
    flow_ca_in=128e-6;
    y_O2_ca_in=1;%%y_O2_amb
    [flow_air_ca_in,flow_water_ca_in,flow_O2_ca_in,flow_N2_ca_in]=mflow(y_O2_ca_in,p_water_ca_in,p_air_ca_in,flow_ca_in);
  
    %% 阴极出口流量
    if p_ca<p_ca_out
        flow_ca_out=0;
    else
        flow_ca_out=k_ca_out*(p_ca-p_ca_out)/p_amb;
    end
    y_O2_ca=p_O2_ca/p_air_ca;
    [flow_air_ca_out,flow_water_ca_out,flow_O2_ca_out,flow_N2_ca_out]=mflow(y_O2_ca,p_water_ca,p_air_ca,flow_ca_out);
    flow_O2_react=num_cell*M_O2*current/4/F;
    flow_water_gen=num_cell*M_water*current/2/F;
   
    %% 阳极入口压力
    p_sat_an_in=psat(T_an_in);
    p_water_an_in=Rhum_an_in*p_sat_an_in;
    p_mix_an_in=p_an_in-p_water_an_in;

    [rho_H2_an_in,cv_H2_an_in,cp_H2_an_in]=PropH2(p_mix_an_in,T_an_in); 
    [rho_wv_an_in,cv_wv_an_in,cp_wv_an_in]=PropWV(p_water_an_in,T_an_in);
    [rho_wl_an_in,cv_wl_an_in,k_wl_an_in,mu_wl_an_in]=PropWL(T_an_in);
    
    %% 阳极内部压力
    p_sat_an=psat(T_an);
    [m_wv_an,m_wl_an]=mvap(m_water_an,p_sat_an,volu_an,R_water,T_an);
    p_H2_an=stateEq(m_H2_an,R_H2,T_an,volu_an);
    p_N2_an=stateEq(m_N2_an,R_N2,T_an,volu_an);
    p_water_an=stateEq(m_wv_an,R_water,T_an,volu_an);
    p_mix_an=p_H2_an+p_N2_an;
    p_an=p_mix_an+p_water_an;
    Rhum_an=p_water_an/p_sat_an;

    [rho_H2_an,cv_H2_an,cp_H2_an]=PropH2(p_H2_an,T_an);
    [rho_N2_an,cv_N2_an,cp_N2_an]=PropN2(p_N2_an,T_an);
    [rho_wv_an,cv_wv_an,cp_wv_an]=PropWV(p_water_an,T_an);
    [rho_wl_an,cv_wl_an,k_wl_an,mu_wl_an]=PropWL(T_an);
          
    %% 阳极入口流量
    flow_an_in=16e-6;
    y_H2_an_in=1;
    [flow_mix_an_in,flow_water_an_in,flow_H2_an_in,flow_N2_an_in]=mflow2(y_H2_an_in,p_water_an_in,p_mix_an_in,flow_an_in);
 
    %% 阳极出口流量
    if p_an<p_an_out
        flow_an_out=0;
    else
        flow_an_out=k_an_out*(p_an-p_an_out)/p_amb;
    end
    y_H2_an=p_H2_an/p_mix_an;
    [flow_mix_an_out,flow_water_an_out,flow_H2_an_out,flow_N2_an_out]=mflow2(y_H2_an,p_water_an,p_mix_an,flow_an_out);
    flow_H2_react=num_cell*M_H2*current/2/F;
  
    %% 膜水合流量
    [sigma_ca,c_ca]=memW(Rhum_ca);
    [sigma_an,c_an]=memW(Rhum_an);
    sigma_ave=(sigma_an+sigma_ca)/2;
    n_d=0.0029*(sigma_ca+sigma_an)^2/4+0.05*sigma_ave-3.4e-19;
    D_w=1.25e-6*exp(2416*(1/303-1/T_cell));
    currentDen=current/area_mem; 
    N_water_mem=n_d*currentDen/F-D_w*(c_ca-c_an)/thick_mem;
    flow_water_mem=num_cell*N_water_mem*M_water*area_mem;
   
    %% 电压
    volt_nst=1.229-8.5e-4*(T_cell-298.15)+R*T_cell*(log(p_H2_an)+0.5*log(p_O2_ca))/(2*F);
    volt_act=0.9514-3.12e-3*T_cell+1.87e-4*T_cell*log(current)-7.4e-5*T_cell*log((p_O2_ca/p_amb)/(5.08e6*exp(-498/T_cell)));
%     resistivity_mem=(thick_mem/area_mem)*181.6*(1+0.03*currentDen+0.062*(T_cell/303)^2*currentDen^2.5)/((sigma_an-0.634-3*currentDen)*exp(4.18*(1-303/T_cell)));
    resistivity_mem=0.01605-3.5e-5*T_cell+8e-5*current;
    volt_ohm=current*resistivity_mem;
    volt_cell=volt_nst-volt_act-volt_ohm;
    power_cell=volt_cell*current;
    eff_cell=(flow_H2_react/flow_H2_an_in)*(volt_cell/1.48);
    volt_cell
    %% 微分方程
    dy=zeros(9,1);
    dy(1)=flow_O2_ca_in-flow_O2_ca_out-flow_O2_react;
    dy(2)=flow_N2_ca_in-flow_N2_ca_out;
    dy(3)=flow_water_ca_in-flow_water_ca_out+flow_water_gen-flow_water_mem;
    dy(4)=flow_H2_an_in-flow_H2_an_out-flow_H2_react;
    dy(5)=flow_N2_an_in-flow_N2_an_out;
    dy(6)=flow_water_an_in-flow_water_an_out-flow_water_mem;
    dy(7)=(h_ca*(T_cell-T_ca)+(cp_air_ca_in*flow_air_ca_in+cp_wv_ca_in*flow_water_ca_in)*T_ca_in-(cp_O2_ca*(flow_O2_ca_out+flow_O2_react)+cp_N2_ca*flow_N2_ca_out+cp_wv_ca*flow_water_ca_out)*T_ca)/(cp_O2_ca*m_O2_ca+cp_N2_ca*m_N2_ca+cp_wv_ca*m_wv_ca+cv_wl_ca*m_wl_ca);
    dy(8)=(h_an*(T_cell-T_an)+(cp_H2_an_in*flow_H2_an_in+cp_wv_an_in*flow_water_an_in)*T_an_in-(cp_H2_an*(flow_H2_an_out+flow_H2_react)+cp_N2_an*flow_N2_an_out+cp_wv_an*flow_water_an_out)*T_an)/(cp_H2_an*m_H2_an+cp_N2_an*m_N2_an+cp_wv_an*m_wv_an+cv_wl_an*m_wl_an);
    dy(9)=(h_ca*(T_ca-T_cell)+h_an*(T_an-T_cell)+h_amb*(T_amb-T_cell)+h_cl*(T_cl-T_cell)+cp_H2_an*flow_H2_react*(T_an-T_cell)+cp_O2_ca*flow_O2_react*(T_ca-T_cell)-cp_wv_ca*flow_water_gen*(T_ca-T_cell)+flow_H2_react*dH_react-num_cell*volt_cell*current)/(Cp_stack*m_cell*num_cell);
end

