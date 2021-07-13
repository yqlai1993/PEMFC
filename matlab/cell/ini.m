function [ m ] =ini( p,volu,T,Rhum,cv )
%UNTITLED3 此处显示有关此函数的摘要
%% 工质摩尔质量
R=8.3145;
M_H2=2*10^-3;
M_O2=32*10^-3;
M_N2=28*10^-3;
y_O2=0.01;
y_H2=0.01;
M_water=18*10^-3;
if cv=='ca'
    M_gas=y_O2*M_O2+(1-y_O2)*M_N2;
else if cv=='an'
        M_gas=y_H2*M_H2+(1-y_H2)*M_N2;
    else
        print('cv error');
    end
end

%% 气体和水蒸气压力
p_sat=psat(T);
p_water=Rhum*p_sat;
p_gas=p-p_water;

%% 气体和水蒸气质量，阴极气体中含氧气和氮气，阳极气体只含氢气。
m_water=p_water*volu*M_water/(R*T);
m_gas=p_gas*volu*M_gas/(R*T);
if cv=='ca'
    x_O2=y_O2*M_O2/M_gas;
    x_N2=(1-y_O2)*M_N2/M_gas;
    m_O2=m_gas*x_O2;
    m_N2=m_gas*x_N2;
    m=[m_water,m_O2,m_N2];
else if cv=='an'
        x_H2=y_H2*M_H2/M_gas;
        x_N2=(1-y_H2)*M_N2/M_gas;
        m_H2=m_gas*x_H2;
        m_N2=m_gas*x_N2;
        m=[m_water,m_H2,m_N2];
    else
       print('m error'); 
    end
end

end

