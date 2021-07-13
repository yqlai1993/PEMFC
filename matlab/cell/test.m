clear all;
clc;


% medium=5;
% t=273.15+80;
% p=101;
% q=1;
% [C,D,H,O,L,V]=refpropm('CDHOLV','T',t,'P',p,'water');
% [ rho,cv,cp,h,mu,k ] = props( medium,t,p,q );
% 
% H_liq=refpropm('H','T',t,'Q',0,'water');
% H_vap=refpropm('H','T',t,'Q',1,'water');
% dH=H_vap-H_liq;


% p=101000;
% volu=0.005;
% T=298.15;
% Rhum=0.5;
% cv='ca';
% m=ini( p,volu,T,Rhum,cv );


% current=20;
% T_cell=273.15+60;
% p_O2_ca=2.08e4;
% volt_act=0.9514-3.12e-3*T_cell+1.87e-4*T_cell*log(current)-7.4e-5*T_cell*log((p_O2_ca/101000)/(5.08e6*exp(-498/T_cell)));


clear all;
%≥ı÷µ
volu_ca=0.01;
volu_an=0.005;
T_ca=273.15+23.5;
T_an=273.15+23.5;
T_cell=273.15+23.5;
Rhum_ca_in=0.0;
Rhum_an_in=0.0;

m_ca=ini(114000*3,volu_ca,T_ca,Rhum_ca_in,'ca');
m_O2_ca=m_ca(2);
m_N2_ca=m_ca(3);
m_water_ca=m_ca(1);

m_an=ini(114000*3,volu_an,T_an,Rhum_an_in,'an');
m_H2_an=m_an(2);
m_N2_an=m_an(3);
m_water_an=m_an(1);

%«ÛΩ‚
tspan=[1,3600];
y=[m_O2_ca,m_N2_ca,m_water_ca,m_H2_an,m_N2_an,m_water_an,T_ca,T_an,T_cell];
[t,y]=ode23('stack',tspan,y);
plot(t,y(:,1),'r',t,y(:,2),'g',t,y(:,3),'b');
figure;
plot(t,y(:,4),'r',t,y(:,5),'g',t,y(:,6),'b');
figure;
plot(t,y(:,7),'r',t,y(:,8),'g',t,y(:,9),'b');



