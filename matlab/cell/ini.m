function [ m_wv,m_g1,m_g2 ] =ini( y1,p,volu,T,Rhum,cv )
%UNTITLED3 此处显示有关此函数的摘要
%% 工质摩尔质量
R=8.3145;
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

%% 气体压力
p_sat=psat(T);
p_wv=Rhum*p_sat;
p_gas=p-p_wv;
%% 气体质量
m_wv=p_wv*volu*M_water/(R*T);
m_gas=p_gas*volu*M_gas/(R*T);
m_g1=m_gas*y1*M1/M_gas;
m_g2=m_gas*y2*M2/M_gas;


end

