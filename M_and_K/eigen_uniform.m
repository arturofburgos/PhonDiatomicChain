%FEM code for a dynamic periodic euler bernoulli beam
%this is the modified code  to
%model the phononic material
%nu=[5.1824,1.2956]; %mass per unit length for both the rods
%E=[210*10^9,210*10^9]; %Young's modulus for both the rods
%I=[0.0000000347,0.0000000018847]; %moment of inertia for both the rods
%fs=1.21699;
function [phon_freq, efq, kmega, Mmega] = eigen_uniform(M,K,m,k,nrods)

M=M;
K=K;
mr=m;
kr=k;


nrods=nrods; %no of rods
Nr=2; %no. of elements per rod

Lstar=2/nrods;
newno=((1/(2*pi))*sqrt(((kr)*K)/((mr)*M)));

nu=[M,mr*M]; %mass per unit length for both the rods
E=[K,kr*K]; %Young's modulus for both the rods
I=[1,1]; %moment of inertia for both the rods

b(1)=E(1)*I(1); %flexural rigidities
b(2)=E(2)*I(2); 
%FEM specs 
Ntotal= (nrods*Nr)+1;
%Ntotal=66;
ne = Ntotal-1; %number of elements
%L=0.18*nrods; %length of the total beam
L=1;
h= (L/ne); %length of an element(same for both type of rods)
delt=0.005; %time step size
t_max= 50; %time upto which you want to run the simulation
it_max=(t_max/delt);
k_element=zeros(4,4);% a hermite cubic basis function has 4 degrees of freedom
k_element=[6 -3*h -6 -3*h; -3*h  2*h*h  3*h  h*h; -6 3*h 6 3*h; -3*h h*h 3*h 2*h*h];%local stiffness matrrix
k_element1= ((2*b(1))/(h*h*h))*k_element;%local stiff matrix for type1 rod
k_element2= ((2*b(2))/(h*h*h))*k_element; %local stiff matrix for type 2 rod

Nmega=2*Ntotal; %size of the primary variable
kmega=zeros(Nmega,Nmega);

flag=1;
    i=1;
    while((i+3)<=(Nmega))
        
        if(flag<=Nr)
    kmega(i:3+i,i:3+i)=kmega(i:3+i,i:3+i)+k_element1;
    i=i+2;
 
         flag=flag+1; 
        end
        
        if(flag>Nr && flag<=(2*Nr))
    kmega(i:3+i,i:3+i)=kmega(i:3+i,i:3+i)+k_element2;
    i=i+2;
   
         flag=flag+1; 
        end
        
       if (flag>2*Nr)
       
       flag=1;
       end
        
    end


M_element=[156 -22*h 54 13*h ; -22*h 4*h*h -13*h -3*h*h; 54 -13*h 156 22*h; 13*h -3*h*h 22*h 4*h*h];
M_element1= nu(1)*((h)/(420))*M_element;
M_element2= nu(2)*((h)/(420))*M_element;
j=1;

Mmega=zeros(Nmega,Nmega);

flag=1;
    i=1;
    while((i+3)<=(Nmega))
        if(flag<=Nr)
    Mmega(i:3+i,i:3+i)=Mmega(i:3+i,i:3+i)+M_element1;
    i=i+2;
 
         flag=flag+1; 
        end
        
        if(flag>Nr && flag<=(2*Nr))
    Mmega(i:3+i,i:3+i)=Mmega(i:3+i,i:3+i)+M_element2;
    i=i+2;
   
         flag=flag+1;  
        end
        
       if (flag>2*Nr)
       
       flag=1;
       end
        
    end
    N=Nmega/2;
    kinitialcommon=kmega;
    minitialcommon=Mmega;
 kmega=kmega((3:((2*N)-2)),(3:((2*N)-2))); %clamped conditions
Mmega=Mmega((3:((2*N)-2)),(3:((2*N)-2))); 
%now we have Mmega and kmega
n_eigs = 20;

disp(kmega)
disp(Mmega)

[V,D] = eigs( kmega, Mmega, n_eigs, 0);
efq = sqrt(diag(D))/(2*pi);
%disp(efq)
%[V,D] = polyeig( kmega,0, Mmega );
%plot(V(1:2:157,1));
%hold on
phon_freq=efq(1);
%plot(V(1:2:157,2));
%hold on
%plot(V(1:2:157,3));
%hold on
%plot(V(1:2:157,4));
%hold on
%disp(efq(4));
%set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 14);
%Leg = legend('Mode 1','Mode 2','Mode 3','Mode 4');
%set( Leg, 'interpreter', 'latex' )
%xlabel( '', 'fontsize', 14, 'interpreter', 'latex')
%ylabel( 'Frequency', 'fontsize', 14, 'interpreter', 'latex')
end