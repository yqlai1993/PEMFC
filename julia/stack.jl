function pipe(flow,p_in,lambda,my,rho)
  length=
  width=
  height=
  area=width*height
  peri=2*(width+height)
  diam=area/peri
  v=flow/area
  
  Re=rho*v*diam/my
  if Re<3200
     f=Re/64
  else
     roughness=0.061e-3
     f=(0.2479-0.0000947*(7-log10(Re))^4)/(log10(roughness/(3.615*diam)+7.366/Re^0.9142))^2
  end
  p_out=p_in-(f*length/diam+lambda)*rho*v^2/2
end




for i in 1:N 
   if i<N
      lambda=0.1
   else
      lambda=0
   end
   
   my=PropsSI('C','T',T_in,'P',p_in,"air")
   rho=PropsSI('D','T',T_in,'P',p_in,"air")
   p=pipe(flow,p_in,lambda,my,rho)
   flow=flow-flow_branch
   p_in=p
   
end

connection=[0~flow_branch-flow_cell]