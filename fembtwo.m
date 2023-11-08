clear all; clc;



%Define our geometry
geometry = @circleg ;
hmax=[1/5,1/20,1/40];

%Initial Values
alpha=4;
delta=0.0005;

for i =1:length(hmax)
    [p ,e , t ] = initmesh( geometry , 'hmax' , hmax(i) );
    
    

     %%%%%% FEM SOLVER %%%%%%
     xi_old=1+20*rand(length(p(1,:)),1); % solution (t=0)
     xi_final=xi_old;
     S=assemblerb(xi_old);
     kl=0.01; %time step
     time=0;
     T=2;

     %Initial values for pop.rate
     pop_counter=1;
     time_matrix=[];
     
     
     [A,M]=assembler(p,t);
     A=delta*A;



     

       while time<T %Crank-Nicholson
           S=assemblerb(xi_old);
           xi_final=(M/kl-M/2+A)\(M*xi_old/kl+M*xi_old/2-M*S-A*xi_old);
           
               
           xi_old=xi_final;             
           time=time+kl;


       

       %%%%%%% Population Rate %%%%%%%%%%
                 integral=0;
                
                pop(pop_counter)=assemblerc(xi_old,integral,t,p);
                time_matrix(pop_counter)=time;
                
                pop_counter=pop_counter+1;
       end

               
     

       


        figure()
        pdesurf(p,t,xi_old)
        
        figure()
        plot(time_matrix,pop)
        xlabel("time")
        ylabel("population")

        







  
end

 



%%%%%%%%%%%%%%%%%%%%%%%%% Functions used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%Computation of aplha and hat-gradients

function [b,c] = gradients(x,y,KK)
    
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/KK;
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/KK;
end





%Stiffness matrix assembler and Load vector assembler

function [A,M]=assembler(p,t)

    n=size(p,2); %number of nodes
    m=size(t,2); %number of elements
    A=sparse(n,n); %Initialize stiffness matrix
    M=sparse(n,n); %Initialize mass matrix
   
    

    for K=1:m
        nodes=t(1:3,K);
        x=p(1,nodes); % x coord of nodes
        y=p(2,nodes); % y coord of nodes
        KK=polyarea(x,y);
       [b,c] = gradients(x,y,KK); %computation of gradients of hat functions  
       
       AK = (b*b'+c*c')*KK; % element stiffness matrix
       A(nodes,nodes) = A(nodes,nodes)+AK;
       
       MK=KK/12*[2,1,1; 1,2,1;1,1,2];% element mass matrix
       M(nodes,nodes)=M(nodes,nodes)+MK;
    end
       


    end


function S=assemblerb(xi)
    alpha=4;
    S=xi.^2+(xi./(xi+alpha));
 end


 
 function pop=assemblerc(xi,integral,t,p)

    
    m=size(t,2);
    for K=1:m
      nodes=t(1:3,K);
      x=p(1,nodes); % x coord of nodes
      y=p(2,nodes); % y coord of nodes
      KK=polyarea(x,y);
      integral=integral+(KK/3)*sum(xi(nodes));
    end
    pop=integral;

end