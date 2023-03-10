%routine to compute equivalent uniform mass and stiffness




%input phononic parameters
M=1000;
K=0.04;
m=0.01;
k=0.1;
nrods=2;


Lu=2/nrods;
[phon_freq,efq, kmega, Mmega]=eigen_uniform(M,K,m,k,nrods);

phon_freq;
%now we need to find the uniform version of this phononic beam
one=1;
Meffguess=(M+(m*M))/2;
Keffguess=(K+(k*K))/2;
ratioguess=Keffguess/Meffguess;


%fun=@(ratioguess)((4.730041*4.730041*(1/(2*pi))*sqrt(ratioguess))-phon_freq);
fun1=@(ratioguess)((4.730041*4.730041*4.730041*4.730041*(1/(2*pi))*(1/(2*pi))*(ratioguess))-((phon_freq)^2));
%freq=@getuniformnaturalfrequency;


finalratio=fzero(fun1,ratioguess);
KbyM=finalratio;
disp(KbyM); 




