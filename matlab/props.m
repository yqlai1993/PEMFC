function [ rho,cv,cp,h,mu,k ] = props( medium,t,p,q )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
if medium==1
        fluid='air.ppf';
        else if medium==2
                fluid='O2.fld';
            else if medium==3
                    fluid='N2.fld';
                else if medium==4
                        fluid='H2.fld';
                    else if medium==5
                            fluid='water.fld';
                        else
                            print('medium error');
                        end
                    end
                end
            end
end
if medium==5&&q==1
 [cp,rho,h,cv,k,mu]=refpropm('CDHOLV','T',t,'Q',q,fluid); 
else
 [cp,rho,h,cv,k,mu]=refpropm('CDHOLV','T',t,'P',p,fluid);
end
       

end

