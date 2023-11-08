clear all; clc ; close all


%Boundary conditions
a=-1;
b=1;

%Initial values
delta=0.1;
TOL=10^(-5);
N=12;
lamda=0.9; %for refinement



% spatial discretization:
xvec=linspace(a,b,N+1);
NN=[];%nodes at each iteration
NN(1)=N;

%Initialize solution vectors
xi=delta*A(xvec)\B(xvec);
zeta=M(xvec)\(-A(xvec)*xi);

%Initialize eta2
eta2=eta2_as(xvec,delta,zeta);

%Initialize sum of errors
summa_eta=[];
summa_eta(1)=sum(sqrt(eta2));


 
 while sum(eta2)>TOL          
    %Refinement of the elements with biggest contribution error
    for i = 1:length(eta2)
        if eta2(i) > lamda*max(eta2) % if large residual
            xvec = [xvec (xvec(i+1)+xvec(i))/2]; % insert new node point
        end
    end
     xvec = sort(xvec); % sort node points accendingly
     
     NN=[NN,length(xvec)];
     

  %Compute uh ,laplace
   xi=delta*A(xvec)\B(xvec);
   zeta=M(xvec)\(-A(xvec)*xi);

  %Compute eta2(uh)
   eta2=eta2_as(xvec,delta,zeta);
   summa_eta=[summa_eta,sum(sqrt(eta2))];

 end


   xi=delta.*A(xvec)\B(xvec);
   zeta=M(xvec)\(-A(xvec)*xi);
   %f(xvec)
   F= arrayfun( @(x) f(x), xvec);

 
 
% %Plot the solution uh
     figure(1);
     subplot(2,2,1)
     xlabel('x domain')
     plot(xvec,xi,'DisplayName','uh')
     legend show
% 
% % % %PLot R(uh)
%     %Compute Residual
     R_uh=F+delta*zeta';
     subplot(2,2,2)
     xlabel('x domain')
     plot(xvec,R_uh,'DisplayName','R(uh)')
     legend show
% 
% %Plot Grid Size distribution
   subplot(2,2,3)
   plot(xvec(2:end),1./diff(xvec),'DisplayName','Grid size Distribution')
   legend show
% 
% %Plot error indicator
    subplot(2,2,4)
   plot(xvec(2:end),(eta2)','DisplayName','eta(uh)')
    legend show

 %Plot nodes-sum of errors at each iteration
     figure(6);
     
     plot(log10(NN),log10((summa_eta)),'DisplayName','eta_v')
     hold on
     plot(log10(NN),log10(1./(NN)),'DisplayName','N^-1 ')
     hold on
    
     legend show



 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%% All functions used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function f(x)
function f = f(x)
 g=10.*x.*sin(7*pi.*x);
 if g>abs(x) 
     f=abs(x);
 
 elseif g<-1*abs(x)
     f=-1*abs(x);

 else
     f=g;
 end
end
 
 function B=B(x)
% Input is a vector x of node coords.
N = length(x) - 1; 
B = zeros(N+1, 1); 
    for i = 1:N
        h = x(i+1) - x(i);
        n = [i i+1];
        B(n) = B(n) + [f(x(i)); f(x(i+1))]*h/2;
    end
 end

function M=M(x)
    %
    % Returns the assembled mass matrix A.
    % Input is a vector x of node coords.
    %
    N = length(x) - 1;
    M = zeros(N+1, N+1); % Initialize matrix to zero
    for i = 1:N % loop over elements
        h = x(i+1) - x(i);
        n = [i i+1];
        M(n,n) = M(n,n) + [1/3 1/6; 1/6 1/3].*h; % Our diagonal consists of 2h/3, our others consist of h/6
    end
    
    M(1,1)=1.e+6;
    M(N+1,N+1)=1.e+6;
end
% Returns the assembled stiffness matrix A.
function A=A(x)
    % Input is a vector x of node coords.
    N = length(x) - 1; % number of elements
    A = zeros(N+1, N+1); % initialize stiffnes matrix to zero
    for i = 1:N % loop over elements
        h = x(i+1) - x(i); % element length
        n = [i i+1]; % nodes
        A(n,n) = A(n,n) + [1 -1; -1 1]/h; % assemble element stiffness
    end
    A(1,1)=1.e+6;

    A(N+1,N+1)=1.e+6;
end


%Calculation of element residuals
function eta2=eta2_as(x,delta,zeta)

N=length(x)-1;
eta2= zeros(N,1); % allocate element residuals

for i = 1:N  % loop over elements
    h = x(i+1)-x(i);
    a2 = f(x(i));% temporary variables
    b2 = f(x(i+1));
    t = (a2^2+b2^2)*h/2; % integrate f^2. Trapezoidal rule
    eta2(i) = (h*sqrt(t))^2; % element residual
end          
 
end
