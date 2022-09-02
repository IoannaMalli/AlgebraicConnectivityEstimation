 %-----------average_consensus_estimator_and_new_update_law-----------------%
clear all
%-------------------------DECENTRALIZED-VERSION----------------------------%
%----------CONSTANTS-------------------------------------------------------%
N= 10 ; 
alpha = 1/(2*(N-1));
r1 = 10 ; 
%-----------------------INITIALIZATION-------------------------------------%

%conds = [1;1;1;32;57;4;2;2;1;14;10;24;1;2;1;25;16;20;2;3;1;4;14;28;1;3;1;40;5;14;1;3;1;31;50;27;1;2;1;33;26;
 %   28;1;1;1;6;53;18;2;1;1;25;22;16;1;1;1;5;12;8;1;1;1;38;30;15;1;3;1;5;47;12;1;2;1;6;57;29;2;1;1;15;
  %  50;1;1;1;1;30;39;14];

Cond = zeros(6,N);

for i=1:N
Cond(1:6,i)=[randi(2),randi([9,30]),randi(1),randi(4),randi(25),randi(23)];
end 
%---------------CHECKING-INITIAL-SYS-FOR-CONNECTION------------------------%
xi = Cond(5,:) ; 
yi = Cond(6,:) ; 
a3 = -2/((r1 - 20)*(r1^2 - 40*r1 + 400)) ; 
a2 = (3*(r1 + 20))/((r1 - 20)*(r1^2 - 40*r1 + 400)) ;
a1 =-(120*r1)/((r1 - 20)*(r1^2 - 40*r1 + 400)) ;
a0 = (400*(3*r1 - 20))/((r1 - 20)*(r1^2 - 40*r1 + 400));
for i=1:N
    for j=1:N
        r=sqrt((xi(i)-xi(j))^2+(yi(i)-yi(j))^2) ;
        if i == j 
            continue ; 
       end 
       if r<r1 
          A_n(i,j) = 1 ;  
       elseif r<20
           A_n(i,j)=a3*r^3+a2*r^2+a1*r+a0 ;  
      end
    end
end 
D = diag(sum(A_n(1:N,:))) ; 
L = D-A_n ; %Laplacian matrix%
e = eig(L) ; 
if e(2)<3 
    disp('not connected')
else 
%-----------------------SOLVE-THE-ODE--------------------------------------%
conds = [] ;
for i=1:N
conds = vertcat(conds,Cond(1:6,i)) ;
end 

ode = @(t1,a)  [agent(t1,X(1,a),neig(1,a,A(a,t1),conds));agent(t1,X(2,a),neig(2,a,A(a,t1),conds));agent(t1,X(3,a),neig(3,a,A(a,t1),conds));agent(t1,X(4,a),neig(4,a,A(a,t1),conds));
    agent(t1,X(5,a),neig(5,a,A(a,t1),conds));agent(t1,X(6,a),neig(6,a,A(a,t1),conds));agent(t1,X(7,a),neig(7,a,A(a,t1),conds));agent(t1,X(8,a),neig(8,a,A(a,t1),conds));agent(t1,X(9,a),neig(9,a,A(a,t1),conds));
    agent(t1,X(10,a),neig(10,a,A(a,t1),conds))];

opts = odeset('RelTol',1e-9,'AbsTol',1e-6);
[t1,a] = ode15s(ode,[0,5],conds,opts) ;

conds = [transpose(a(end,1:10));-150;-150;transpose(a(end,13:22));150;150;transpose(a(end,25:60))];
ode = @(t2,b)   [agent(t2,X(1,b),neig(1,b,A(b,t2),conds));agent(t2,X(2,b),neig(2,b,A(b,t2),conds));agent(t2,X(3,b),neig(3,b,A(b,t2),conds));agent(t2,X(4,b),neig(4,b,A(b,t2),conds));
    agent(t2,X(5,b),neig(5,b,A(b,t2),conds));agent(t2,X(6,b),neig(6,b,A(b,t2),conds));agent(t2,X(7,b),neig(7,b,A(b,t2),conds));agent(t2,X(8,b),neig(8,b,A(b,t2),conds));agent(t2,X(9,b),neig(9,b,A(b,t2),conds));
    agent(t2,X(10,b),neig(10,b,A(b,t2),conds))];

opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t2,b] = ode15s(ode,[0,5],conds,opts) ;

conds = [transpose(b(end,1:10));150;150;transpose(b(end,13:18));conds(19:22);10;10;transpose(b(end,25:60))];

ode = @(t3,c)  [agent(t3,X(1,c),neig(1,c,A(c,t3),conds));agent(t3,X(2,c),neig(2,c,A(c,t3),conds));agent(t3,X(3,c),neig(3,c,A(c,t3),conds));agent(t3,X(4,c),neig(4,c,A(c,t3),conds));
    agent(t3,X(5,c),neig(5,c,A(c,t3),conds));agent(t3,X(6,c),neig(6,c,A(c,t3),conds));agent(t3,X(7,c),neig(7,c,A(c,t3),conds));agent(t3,X(8,c),neig(8,c,A(c,t3),conds));agent(t3,X(9,c),neig(9,c,A(c,t3),conds));
    agent(t3,X(10,c),neig(10,c,A(c,t3),conds))];

opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t3,c] = ode15s(ode,[0,5],conds,opts) ;

%--------------------------EIGENVALUE-ESTIMATION---------------------------%
t2 = t2+5;
t3 = t3+10;
t = vertcat(t1(1:(end-1)),t2(1:(end-1)),t3);
af = vertcat(a(1:(end-1),:),b(1:(end-1),:),c) ; 
for n=1:N
l2n(:,n) = af(:,(n-1)*6+3)./af(:,(n-1)*6+2) ;
end 

%--------------------------ACCURATE-EIGENVALUE-----------------------------%
a3 = -2/((r1 - 20)*(r1^2 - 40*r1 + 400)) ; 
a2 = (3*(r1 + 20))/((r1 - 20)*(r1^2 - 40*r1 + 400)) ;
a1 =-(120*r1)/((r1 - 20)*(r1^2 - 40*r1 + 400)) ;
a0 = (400*(3*r1 - 20))/((r1 - 20)*(r1^2 - 40*r1 + 400));
B = zeros(N) ;
for k = 1:length(t) 
for i=1:N
    if i==2 && t(k)>5
        continue;
    end 
    if i==4 && t(k)>5 && t(k)<10
        continue;
    end 
    for j=1:N
        r=sqrt((af(k,(i-1)*6+5)-af(k,(j-1)*6+5))^2+(af(k,6*i)-af(k,6*j))^2) ;
        if i == j  
            continue ; 
        end 
       if r<r1 
          B(i,j) = 1 ;  
       elseif r<20
           B(i,j)= a3*r^3+a2*r^2+a1*r+a0 ;  
       else 
           B(i,j)=0 ;
       end 
    end 
end
D = diag(sum(B(1:N,:))) ; 
L = D-B ; 
e = eig(L) ;
if t(k)<=5
   e1_N(k) = e(2) ; 
elseif t(k)<=10
   e1_N(k) = e(4) ; 
else 
   e1_N(k) = e(3) ; 
end 
end


%plot%

%ploting eigenvalues%
figure


plot(t,l2n(:,1)) ;
hold on
plot(t1,l2n(1:length(t1),2)) ;
hold on
plot(t,l2n(:,3));
hold on
plot(t1,l2n(1:length(t1),4)) ;
hold on
plot(t,l2n(:,5:N));
hold on
plot(t,e1_N,':');
xlabel('time')
ylabel('eigenvalue estimation')

hold off

%-----------------------GRAPHS-PLOT----------------------------------------%
k=length(t1) ;
for i=1:N
    if i==2 && t(k)>5
        continue;
    end 
    if i==4 && t(k)>5 && t(k)<10
        continue;
    end 
    for j=1:N
        r=sqrt((af(k,(i-1)*6+5)-af(k,(j-1)*6+5))^2+(af(k,6*i)-af(k,6*j))^2) ;
        if i == j  
            continue ; 
        end 
       if r<r1 
          C(i,j) = 1 ;  
       elseif r<20
           C(i,j)= a3*r^3+a2*r^2+a1*r+a0 ;  
       else 
           C(i,j)=0 ;
       end 
    end 
end

figure
subplot(1,3,1)
plot(graph(A_n))
subplot(1,3,3)
plot(graph(B))
subplot(1,3,2)
plot(graph(C))
end
%----------------------FUNCTIONS-------------------------------------------%
%-----------------------A-MATRIX----------------------------------------%
function Aij = A(a,t) 
N=10 ; 
r1 = 10 ;
xi = a(((1:N)-1)*6+5) ;
yi = a(6*(1:N));
a3 = -2/((r1 - 20)*(r1^2 - 40*r1 + 400)) ; 
a2 = (3*(r1 + 20))/((r1 - 20)*(r1^2 - 40*r1 + 400)) ;
a1 =-(120*r1)/((r1 - 20)*(r1^2 - 40*r1 + 400)) ;
a0 = (400*(3*r1 - 20))/((r1 - 20)*(r1^2 - 40*r1 + 400));
Aij = zeros(N) ; 
for i=1:N
    if i==2 && t>5
        continue;
    end 
    if i==4 && t>5 && t<10
        continue;
    end 
    for j=1:N
        r=sqrt((xi(i)-xi(j))^2+(yi(i)-yi(j))^2) ;
        if i == j 
            continue ; 
        end 
       if r<r1 
          Aij(i,j) = 1 ;  
       elseif r<20
           Aij(i,j)= a3*r^3+a2*r^2+a1*r+a0 ;    
       end 
    end 
end 
end 
%------------------------X-FUNCTION----------------------------------------%
function x = X(n,an) 
s = (n-1)*6+1 ;
f = 6*n ; 
x(1:6)=an(s:f) ;
end 
%--------------------CREATION-OF-NEIG-MATRIX-------------------------------%
function m = neig(n,an,A,C) 
N = length(A) ;
%find the neighbors %
neg = [] ; 
for i=1:N 
    if A(i,n) > 10^(-5) 
        neg = horzcat(neg,i) ; 
    end 
end 
lneg = length(neg) ; 
%create m matrix %
m = zeros(10,lneg) ;
for i=1:lneg 
    m(1:6,i) = an(((neg(i)-1)*6+1):neg(i)*6) ; 
    m(7,i) = A(n,neg(i)) ; 
    m(8,i)= 1.5*abs(C((n-1)*6+1)-C((neg(i)-1)*6+1))+0.2;
    m(9,i) = 1.5*abs(C((n-1)*6+2)-C((neg(i)-1)*6+2))+0.2 ; 
    m(10,i) = 1.5*abs(C((n-1)*6+3)-C((neg(i)-1)*6+3))+0.2 ; 
end 
end 
%-------------------------ODE-FUNCTION------------------------------%
function da = agent(t,s,z)

%constants%
N=10 ;
k = 2 ; kr = 150 ; l=100 ; pf = 0.0001 ; 
alpha = 1/(2*(N-1));
beta = 350; 
gamma = 80; 

%variable-assignment%
c_i = horzcat(s(1),z(1,:)) ;
b_i = horzcat(s(2),z(2,:)) ;
a_i = horzcat(s(3),z(3,:)) ; 
x = horzcat(s(4),z(4,:)) ; 
xi = horzcat(s(5),z(5,:)) ; %syntetagmeni x axona %
yi = horzcat(s(6),z(6,:)) ; %syntetagmeni y axona % 

 
for i=1:length(c_i)-1
    %performance-functions%
    rij_1(i) = (z(8,i) - pf)*exp(-l*t)+pf ; %c_i consensus estimator%
    rij_2(i) = (z(9,i) - pf)*exp(-l*t)+pf ; %b_i consensus estimator%
    rij_3(i) = (z(10,i) - pf)*exp(-l*t)+pf ; %a_i consensus estimator%
    %consensus-a_i% 
    ji = (a_i(1)-a_i(i+1))/rij_3(i) ; 
    T3(i) = 1/2*log((1+ji)/(1-ji)) ; 
    J3(i) = (1-ji^2)^(-1);
    %consensus-c_i% 
    ji = (c_i(1)-c_i(i+1))/rij_1(i) ; 
    T1(i) = 1/2*log((1+ji)/(1-ji)) ; 
    J1(i) = (1-ji^2)^(-1);
    %consensus-b_i%
    ji = (b_i(1)-b_i(i+1))/rij_2(i) ; 
    T2(i)=1/2*log((1+ji)/(1-ji)) ; 
    J2(i) = 1/(1-ji^2);
end 
%creation-of-diff-equations-%
dx = alpha*beta*a_i(1)/b_i(1)*x(1)-beta*c_i(1)+beta*c_i(1)^2/b_i(1)*x(1)+gamma*(1-b_i(1))*x(1) ; 
for v=1:length(c_i)-1
    dx = dx-alpha*beta*z(7,v)*(x(1)-x(v+1)) ; 
end 
dc_i = -kr*(c_i(1)-x(1))+dx;
db_i = -kr*(b_i(1)-x(1)^2)+2*x(1)*dx ;
da_i = -kr*a_i(1);
dxi = 0.2;
dyi = 0.5*cos(t) ; 
for v=1:length(c_i)-1
    dc_i = dc_i-k*(J1(v)*T1(v)/rij_1(v));
    db_i = db_i-k*(J2(v)*T2(v)/rij_2(v)) ; 
    da_i = da_i+(kr*x(1)+2*dx)*(z(7,v)*(x(1)-x(v+1)))-k*(J3(v)*T3(v)/rij_3(v)) ;     
end 
da =  [dc_i ; db_i ;da_i; dx ; dxi ; dyi] ; 
end 


