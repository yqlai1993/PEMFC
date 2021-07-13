function [sys,x0,str,ts,simStateCompliance] = drop(t,x,u,flag,l,m,d,k)
switch flag,
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;
    
  case 1,
    sys=mdlDerivatives(t,x,u,l,m,d,k);
    
  case 2,
    sys=mdlUpdate(t,x,u);
    
  case 3,
    sys=mdlOutputs(t,x,u,l,m,d,k);
    
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);
    
  case 9,
    sys=mdlTerminate(t,x,u);
    
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end



function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 2;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed
sys = simsizes(sizes);

x0  = [-1;0];
str = [];
ts  = [0 0];
simStateCompliance = 'UnknownSimState';


function sys=mdlDerivatives(t,x,u,l,m,d,k)
if x(1)>=0
    b=-k*x(1);
else
    b=0;
end
sys = [x(2);10+b/m-1*x(2)/m-1*abs(x(2))*x(2)/m];

function sys=mdlUpdate(t,x,u)
sys = [];

function sys=mdlOutputs(t,x,u,l,m,d,k)
sys = [d-l-x(1)];

function sys=mdlGetTimeOfNextVarHit(t,x,u)
sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

function sys=mdlTerminate(t,x,u)
sys = [];
