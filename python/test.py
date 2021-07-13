print("hello")
import CoolProp.CoolProp as CP 
rho=CP.PropsSI('D','P',101000,'Q',0,'water')
print(rho)

