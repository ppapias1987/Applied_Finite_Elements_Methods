%Aparaithta gia ka8e swsth ektelesh tou script>  
clear all; close all; clc % ka8arizei variables, kleinei plot figures, ka8arizei command window

%Boundary conditions
left = -1;
right = 1;

%Initial values
delta = 0.1;
TOL=1e-3; % epishs to TOL gia 0.001 grafetai sthn matlam ws: 1e-3 dioti 1e-3 == 0.001
N=12;
% h=1/(N-1); % eisai sigouros oti to xreiasomaste edw to h? dioti orizetai mesa sthn loopa -gia apofugh variable conflict sto lew- , an to kaneis gia na oriseis to x-domain ennalaktika einai: x=linspace(a,b,N) 
eta2 = ones(N,1);

% spatial discretization:
% x = (a:h:b); %to metefera sthn grammh 18 giati to x-domain to xrhsimopoieis gia thn g sthn grammh 25 kai gia thn f sthn grammh 25
           % SYNEPWS: PRWTA ORIZEIS to x-domain KAI META to kaleis px sthn
           % grammh 28 alliws sou leei san la8os: "Unrecognized function or variable 'x'."
% enallaktika:
x = linspace(left,right,N);

%function g(x)
g= @(x) 10.*x.*sin(7*pi.*x); % evala ".x" anti gia "x" gia elementwise pol/smo mhn xehnas to x einai vector

%function f(x)
f= function_f(g(x),x); % dior8wsa to g me g(x)

while sum(eta2)>TOL
 eta = zeros(N,1); % allocate element residuals
 
for i=1:length(x) % PROSEXE: edw trexw thn loopa KATA MHKOS olou tou x-domain sunolo 12 shmeia (alliws gia N=12 sthn sxesh x=a:h:b pairnw 23 x-points) 
%  for i = 1:N % loop over elements
    h = x(i+1) - x(i); % element length
    %a = f(x(i)); % temporary variables %PROSEXE thn onomasia twn metavlhtwn soy, edw exeis variable conflict dld to a einai kai left boundary kai temp variable sthn loopa , gamietai o DIAS to antistoixo gia b
    a = f(i) % to f einai HDH VECTOR me tis times f(x(i)) ara gia na exeis prosvash se auth les f(i)
    %b = f(x(i+1));
    b = f(i+1)
    t = (a.^2+b.^2)*h/2 % integrate f^2. Trapezoidal rule
    eta2(i) = h*sqrt(t); % element residual
    
    if i == length(x)-1 %dioti alliws exeis to exhs la8os: %Error in femprojecta (line 33)                                             %h = x(i+1) - x(i); % element length
        break       
    end
    
end

%Refinement procedure
lamda = 0.9; % refinement parameter
for i = 1:length(eta)
    if eta(i) > alpha*max(eta) % if large residual
        x=[x (x(i+1)+x(i))/2]; % insert new node point
    end
end
x = sort(x); % sort node points accendingly

end