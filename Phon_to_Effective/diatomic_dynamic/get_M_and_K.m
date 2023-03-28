%% Routine to determine Meff and Keff. Subjecting a constant unit force at the first mass


% Input the Keff/Meff ratio based on fs considerations
Keff_by_Meff = 1.5239e+06;

% Now enter the displacement rule to be obeyed by the output
disp_uniform = 0.1; % The value used in the abstract.


% Input a guess for Meff and Keff
Meff=1/300;
%Keff=;
%input m,k
m=1;
k=1;

% Number of cells
ncell = 5; 

tguess=Meff;

%fun2=@(tguess)((((get_max_disp(tguess,Keff_by_Meff,tguess,Keff_by_Meff,ncell))^2-(disp_uniform)^2)));
% Same output as below.

fun2=@(tguess)((((get_max_disp(tguess,Keff_by_Meff,tguess,Keff_by_Meff,ncell))-(disp_uniform))));
options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
finalratio2=fzero(fun2,tguess,options);
M=abs(finalratio2);
K=abs(finalratio2*Keff_by_Meff);
disp('M=')
disp(M);
disp('K=')
disp(K);