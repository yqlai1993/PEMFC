function [m_water_vap,m_water_liq]=mvap(m_water,p_sat,volu,R_water,T) 
    m_water_max=p_sat*volu/R_water/T;
    if m_water<=m_water_max
       m_water_vap=m_water;
       m_water_liq=0;
     else
       m_water_vap=m_water_max;
       m_water_liq=m_water-m_water_max;
    end
end