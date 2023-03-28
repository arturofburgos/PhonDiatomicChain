
%input the Keff/Meff ratio based on fs considerations
KbyMratio=3.0379e-05;

%now enter the displacement rule to be obeyed by the output
dispuniform=0.1; %previous case: a good reference

%time for which the beam is forced !caution it should be adequate enough to
%allow the beam to have atleast 1 oscillation
t=10;

%input a guess for Meff
Meff=1;
Keff=1;
%input m,k
m=1;
k=1;

tguess=Meff;
%fun2=@(tguess)((((getmaxdisp(tguess,tguess*KbyMratio,m,k,t))^2-(dispuniform)^2)));
fun2=@(tguess)((((getmaxdisp(tguess,tguess*KbyMratio,m,k,t))-(dispuniform))));
options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
finalratio2=fzero(fun2,tguess,options);
M=abs(finalratio2);
K=abs(finalratio2*KbyMratio);
disp('M=')
disp(M);
disp('K=')
disp(K);
