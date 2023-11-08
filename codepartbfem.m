
clear all;   clc



%Define our geometry
geometry = @circleg ;
hmax=[1/2,1/4,1/8,1/16,1/32];%1/4,1/8,1/16,1/32%;


EnE=zeros(1,length(hmax));
pp=zeros(1,length(hmax));


for i =1:length(hmax)
    [p ,e , t ] = initmesh( geometry , 'hmax' , hmax(i) );
    m=size(t,2); %number of elements
    [A,B,uex]=assembler(p,t);



% %Boundary conditions
      I=eye(length(p));
      A(e(1,:),:) = I(e(1,:),:);
      B(e(1,:))=uex(e(1,:));


    xi=A\B;

      %Plot solutions
%       figure();
%       pdeplot(p,[],t,'XYData',xi,'ZData',xi,'ColorBar','off')

    

    %Error
    err=uex-xi;
    EnE(i)=sqrt(err'*A*err);
    



 
  


end
  figure()
  hold on
  xlabel('hmax')
  loglog(hmax,EnE,'DisplayName','energy norm of error')

  pp=polyfit(log10(hmax),log10(EnE),1);   %Convergence Rate
%   figure(1)
  loglog(hmax,hmax.^pp(1),'DisplayName','hmax^p ,p=1.23')
  legend show

 



%%%%%%%%%%%%%%%%%%%%%%%%% Functions used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%Computation of aplha and hat-gradients

function [b,c] = gradients(x,y,KK)
    
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/KK;
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/KK;
end


function[f]=f(x,y)
    f=8*((pi)^2).*sin(2.*pi.*x).*sin(2.*pi.*y);    
end



%Stiffness matrix assembler and Load vector assembler

function [A,B,uex]=assembler(p,t)

    n=size(p,2); %number of nodes
    m=size(t,2); %number of elements
    A=sparse(n,n);
    B=zeros(n,1);
    uex=zeros(n,1);
    

    for K=1:m
        nodes=t(1:3,K);
        x=p(1,nodes); % x coord of nodes
        y=p(2,nodes); % y coord of nodes
        KK=polyarea(x,y);
       [b,c] = gradients(x,y,KK); %computation of gradients of hat functions
       xc=mean(x); 
       yc = mean(y);
       bk=f(xc,yc)*KK/3;
       AK = (b*b'+c*c')*KK; % element stiffness matrix
       B(nodes)=B(nodes)+bk;
       A(nodes,nodes) = A(nodes,nodes)+AK;
       uex(nodes)=f(x,y)/(8*((pi)^2));
    end
end







