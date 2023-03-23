%% Routine to determine Meff and Keff. Subjecting a constant unit force at the first mass


% Input the Keff/Meff ratio based on fs considerations
Keff_by_Meff = 28.7287;

% Now enter the displacement rule to be obeyed by the output
disp_uniform = 0.0505; % The value used in the abstract This value is the one we found in the displacement for maximum displacement in phononic case.


% Input a guess for Keff

Keff_guess=5;



%fun2=@(tguess)((((get_max_disp(tguess,Keff_by_Meff,tguess,Keff_by_Meff,ncell))^2-(disp_uniform)^2)));
% Same output as below.

fun2=@(Keff_guess)((static_disp(Keff_guess))^2 - (disp_uniform)^2);
options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
Keff=fzero(fun2,Keff_guess,options);
M=Keff/Keff_by_Meff;
disp('Meff=')
disp(M);
disp('Keff=')
disp(Keff);