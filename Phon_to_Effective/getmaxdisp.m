function maxdisp = getmaxdisp(M,K,m,k,t)

nu=[abs(M),m*abs(M)]; %mass per unit length for both the rods
E=[abs(K),k*abs(K)]; %Young's modulus for both the rods
I=[1,1]; %moment of inertia for both the rods
%fs=3;
b(1)=E(1)*I(1); %flexural rigidities
b(2)=E(2) *I(2); 
%FEM specs
nrods=10 ; %no of rods 
Nr=8; %no. of elements per rod
Ntotal= (nrods*Nr)+1;
%Ntotal=66;
ne = Ntotal-1; %number of elements
%L=0.18*nrods; %length of the total beam
L=1;
h= (L/ne); %length of an element(same for both type of rods)
delt=0.005; %time step size
t_max=t; %time upto which you want to run the simulation
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

Ue=zeros((Nmega-4),1); %odd index has deflection and the even indices have the slope of deflection
%initial shape of the beam
%betad=4.730000000000;
%kd=-1.017809410640358;
  

%for m=3:2:((2*N)-2)
%x = ((m-1)/2)*h;
%xb0=x;


%Ue(m-2,1)=-1.00000000000*0.06*((sinh(betad*xb0(j))- sin(betad*xb0(j)))+ (kd*(cosh(betad*xb0(j))- cos(betad*xb0(j)))));
%Ue(m-1,1)=-((-1*0.06*betad*cosh(betad*xb0(j)))-(-1*0.06*betad*cos(betad*xb0(j)))+(-1*0.06*kd*betad*sinh(betad*xb0(j)))+(-1*betad*1*0.06*kd*sin(betad*xb0(j))));
%end

Uini=Ue;
Fint=kmega*Ue;
k=kmega;
M=Mmega;
Unext=zeros((2*N)-4,1);
DUnext=zeros((2*N)-4,1);
DUe=zeros((2*N)-4,1);
%forcing matrix
F = zeros((2*N),1); %free vibration case without considering BC
F0=1; %constant forcing magnitude
%Felement=[6;-h;6;h];
j=1;
%while(j+3<=2*N)
    %F(j:3+j,1)=F(j:3+j,1)+ Felement;
    %j=j+2;
%end
F=F((3:(2*N)-2),1);
F=((F0*h)/(12))*F;

F=zeros(2*N,1);
F=F((3:(2*N)-2),1);
%sanity check steady state deflection
Usteady=k\F;
Velsteady=zeros(N,1);
flag2=1;
for i=1:2:((2*N)-4)
    Velsteady(flag2,1)= Usteady(i,1);
flag2=flag2+1;
end

%newmark time stepping strategy

 DUUe= M\(F-k*Ue);
 acc=DUUe;
ind=1;
khat=zeros((2*N)-4,(2*N)-4);
khat=M+(((delt*delt)/4)*k);
Vesave=zeros(N-2,1000);


while(ind<=it_max)
   
    %F(299)=0.05*sin(50*delt*ind);
%F(1597)=5*cos(220*2*pi*delt*ind);
%F(2199)=0.05*sin(64*delt*ind);
t=delt*ind;

F(79)=1;   
A1=zeros((2*N)-4,1);
A1=(Ue) + (delt*DUe) + (((delt*delt)/4)*DUUe);
Fhat= (((delt*delt)/4)*F);
Mhat=M*A1;
Rhshat=Fhat+Mhat;
Unext= khat\Rhshat ;
    
DUUnext= M\(F-k*Unext);
DUUehalf=0.5*(DUUnext+DUUe);
DUenext=DUe+((delt)*DUUehalf);
    
    
Ue=Unext;
    DUe=DUenext;
DUUe=DUUnext;
   flag=1;
   for i=1:2:((2*N)-4)
   Ve(flag,ind)= Ue(i,1);
   flag=flag+1;
   end
  %reaction force of the first element
    Xe=[0;0;Ue(1,1);Ue(2,1)];
    Xedd=[0;0;DUUe(1,1);DUUe(2,1)];
   Freac(:,ind)=(k_element1*Xe)+(M_element1*Xedd);
   
   
   Fr1(ind,1)=Freac(1,ind);

    ind=ind+1;
%disp(ind)
end
Ve(80,:)=0;


%maximum displacement
maxdisp=max((Ve(40,:)));
%disp(maxdisp)


end