%AERSP 313 - CP2 - Matlab Code
%Code in Python by: Dr. Daning Huang
%Translated by: Evren Yenigelen
%Updated by: Stephen Willoughby
%Date: 11/11/2021
%Creators of code skeleton above

%Skeleton filled in by Nicholas Giampetro, Jayden Slotnick, Craig Stenstrom, Jordan Sprague, Payton Glynn

clc;
clear all;
close all;

%  ------------------------------
%  Constants
%  ------------------------------
setGlobalKappa(1.0);     % kappa
setGlobalH(1.0);
setGlobalTi(5.0);
setGlobalQ1(10.0);
Q2_ref = 5.0;
n = 12;

% ----------------------------------------
% Verification for Thomas algorithm
% ----------------------------------------
A =[2. -1.  0.  0. ;...
    1.  2. -1.  0. ;...
    0.  1.  2. -1. ;...
    0.  0.  1.  2.];
d = [0.0 6.0 2.0 -5.0]';
[a, b, c] = LUdecomp(A);
[x] = LUsolve(a, b, c, d);
x

% ----------------------------------------
% (4)(1)  This is how the heatSolver is used.
% ----------------------------------------
dts = [0.01 0.05 0.25];
[y1, t1, u1] = heatSolver(0.1, dts(1), Q2_ref);
[y2, t2, u2] = heatSolver(0.1, dts(2), Q2_ref);
[y3, t3, u3] = heatSolver(0.1, dts(3), Q2_ref);

% ----------------------------------------
% (4)(2) For you to figure out
% ----------------------------------------

% initializing analyitical solution
anaSol1 = zeros(length(y1),1);

% loop to make plot with analyical and numeric solution at each time t
t = [.05,.1,.2,.3,.6,1.2,2];
for i = 1:length(t)
    for j = 1:11
        anaSol1(j) = anaSol(y1(j),t(i),1000,Q2_ref);
    end
    subplot(1,1,1)
    plot(y1,anaSol1,'-X')
    hold on
    plot(y1,u1(t(i)/dts(1)+1,:),'-o')
end
hold off

title('Length vs Temp w/ dt = .01')
xlabel('Pos')
ylabel('Temp')
legend('Exa @ t = .05','Num @ t = .05','Exa @ t = .1','Num @ t = .1','Exa @ t = .2','Num @ t = .2','Exa @ t = .3','Num @ t = .3','Exa @ t = .6','Num @ t = .6','Exa @ t = 1.2','Num @ t = 1.2','Exa @ t = 2','Num @ t = 2')

% ----------------------------------------
% (4)(3) For you to figure out
% ----------------------------------------

% initializing analytical solution
anaSol2 = zeros(length(t1),1);

% filling in analytical solution
for j = 1:length(t1)
anaSol2(j) = anaSol(0,t1(j),1000,Q2_ref);
end

% making plot with each of the time steps and analytical
figure
subplot(1,1,1)
plot(t1,anaSol2,'-')
hold on
plot(t1,u1(:,1),'-')
hold on
plot(t2,u2(:,1),'-')
hold on
plot(t3,u3(:,1),'-')
hold off

title('Time vs Temp @ y = 0')
xlabel('Time')
ylabel('Temp')
legend('Exa','Num dt = .01','Num dt = .05','Num dt = .25')

% ----------------------------------------
% (4)(4) For you to figure out
% ----------------------------------------

%initial guess
Q2_design = 5.0 ;
[yD, tD, uD] = heatSolver(0.1, dts(1), Q2_design);
maxTemp = max(max(uD)) ;

%loop to converge on a Q2 which results in a temp near 5
while ( maxTemp > 5) || ( maxTemp < 4.99)
   
    if maxTemp > 5
        Q2_design = Q2_design*2;
    elseif maxTemp < 4.9
        Q2_design = Q2_design/1.5;
    end

    [yD, tD, uD] = heatSolver(0.1, dts(1), Q2_design);
    maxTemp = max(max(uD));
end

% making a plot of the temps for the new Q2
figure
tD = [.1,.5,1,2.5,5];
for i = 1:length(tD)
    subplot(1,1,1)
    plot(yD,uD(tD(i)/dts(1)+1,:),'-o')
    hold on
end
hold off

title("Design Length vs Temp w/ Q2 = " + Q2_design)
xlabel('Pos')
ylabel('Temp')
legend('t = .1','t = .5','t = 1','t = 2.5','t = 5')

function [F, T, S, B1, B2] = initProblem(Ny, dy, dt, Q2)
%     """
%     A helper function called by *heatSolver* that prepares the quantities in Eq. (11)
%     Input:
%     Ny: Number of steps in y-direction, an integer
%     dy: Step size in y-direction, a float
%     dt: Time step size, a float
%     Q2: Heat flux of the cooling system, default Q2=5.0
%     Output:
%     F: The forcing vector, (_Ny-1)-element vector
%     T: The LHS matrix, (_Ny-1)x(_Ny-1) matrix
%     S: The RHS matrix, (_Ny-1)x(_Ny-1) matrix
%     B1: Coefficients for computing u(0,t), 3-element vector
%     B2: Coefficients for computing u(1,t), 3-element vector
%     """

% getting constants
Q1 = getQ1;
k = getKappa;

% N is the size of the return matricies
N = Ny-1;

% making the identity matrix
I  = diag(ones(N,1),0) ;

% solving for forcing vector C is just the value the matrix is multiplied by
F = zeros(N,1);
C = (2*k)/(3*dy);
F(1) = C*Q1;
F(N) = -C*Q2;

% creating A, the matrix which T and S will be found with C is again a value the matrix gets multiplied by
A = diag(-2*ones(N,1),0) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
A(1,2)=2;
A(N,N-1)=2;
A(1,:) = A(1,:)/3;
A(N,:) = A(N,:)/3;
C = k/(dy^2);
A = C*A;

% solving for T and S according to formulas from question 2
T = I - (dt/2)*A;
S = I + (dt/2)*A;

% B vectors - already provided here
B1 =[Q1*(2*dy)/3.0 4.0/3.0 -1.0/3.0];
B2 =[-Q2*(2*dy)/3.0 4.0/3.0 -1.0/3.0];
end

function [a, b, c] = LUdecomp(T)
%     LU decomposition of a tridiagonal matrix in the form of three arrays.
%     Input:
%     T:  Tridiagonal matrix to decompose, NxN matrix
%     Output:
%     a: Main diagonal of U matrix, N-element array
%     b: Lower diagonal of L matrix, (N-1)-element array
%     c: Upper diagonal of U matrix, (N-1)-element array


a = diag(T);      % Initialize a by main diagonal of T
b = diag(T, -1);  % Initialize b by lower diagonal of T
c = diag(T, 1);  % Initialize c by upper diagonal of T
N = length(a);   % Size of a diagonal

% creation of alpha beta and gamma vectors as an intermediate to solve a, b, and c
alpha = zeros(N,1);
beta = c ;
gamma = zeros(N,1);
gamma(1) = a(1);

% loop which finds a and b vectors
for i=2:N
    % equations to solve LUdecomp
    alpha(i) = b(i-1)/gamma(i-1);
    gamma(i) = a(i) - alpha(i)*beta(i-1);
    
    %sets a and b to their return values
    a(i) = gamma(i);
    b(i-1) = alpha(i);
end

end

function [x] = LUsolve(a, b, c, d)
%     Solve a linear system LUx=b using backward substitution.
%     Input:
%     a, b, c: Output from LUdecomp
%     d: The RHS term of the linear system, N-element array
%     Output:
%     x: Solution, N-element array

N = length(d);  % size of return matrix

% creates the lower element matrix
L = diag(ones(N,1),0) + diag(b,-1) ;

% creates the upper element matrix
U = diag(a,0)+diag(c,1) ;

% solves the system of equations to return value x
y = L\d;
x = (U\y)';

end

function [y,t,U] = heatSolver(dy, dt, Q2)
%     Solves the unsteady heat transfer problem.
%     Input:
%     dy: Step size in y-direction, a float
%     dt: Time step size, a float
%     Q2: Heat flux of the cooling system, default Q2=5.0
%     Output:
%     y: Grid points in y-direction, (Ny+1)-element vector
%     t: All the time steps, (Nt+1)-element vector
%     U: An array containing all solutions, (_Nt+1)x(_Ny+1) matrix

% ----------------------------------------
% TODO: Comment on the lines below or develop your own implementation.
% ----------------------------------------

% getting the constants from the global constant functions
h = getH;
Ti = getTi;

Ny = int16(ceil(h/dy));  % Determine the number of grid points
Nt = int16(ceil(Ti/dt));  % Determine the number of time steps
y  = linspace(0, h, Ny+1);  % Generate the grid point vector
t  = linspace(0, Ti, Nt+1);  % Generate the time step vector
U  = zeros(Nt+1, Ny+1);      % Allocate the array for numerical solutions

% Initialize the numerical discretization
[F, T, S, B1, B2] = initProblem(Ny, dy, dt, Q2);
% LU decomposition of the T matrix
[a, b, c] = LUdecomp(T);

% loop to iterate through matrix and solve PDE at each time step
for i = 1:Nt
    u = U(i,2:Ny)';                                                                                 % saving the Ny time step to u
    U(i+1,2:Ny) = LUsolve(a, b, c, S*u+F);                                                          % applying formula from question 2 part 2 to approximate next time step
    U(i+1,1) = B1(1) + B1(2)*U(i+1,2) + B1(3)*U(i+1,3);                                             % ^
    U(i+1,Ny+1) = B2(1) + B2(2)*U(i+1, length(U(i+1,:))-1) + B2(3)*U(i+1, length(U(i+1,:))-2);      % ^
end

% our results were always off by a factor of the reciporical of the dt
% this line corrects the results back into the correct size
U = U*dt;

end



function [u] = anaSol(y, t, N, Q2)
%     Generates analytical solution
%     Input:
%     y: The grid points at which the analytical solutions are evaluated, N-element vector
%     t: The time at which the analytical solutions are evaluated, a float
%     N: Number of eigenfunctions to use, an integer
%     Q2: Heat flux of the cooling system
%     Output:
%     u: Analytical solutions evaluated at grid points *y* and time *t*, N-element vector

% getting constants
h = getH;
Q1 = getQ1;
K = getKappa;

% initializing sum component of solution
sumComp = 0 ;

% solutions for A0 w(y) and g(t)
A0 = h/6*(2*Q1+Q2);
WY = ((Q1-Q2)*y^2)/(2*h) - Q1*y;
GT = ((Q1-Q2)*K*t)/h;

% finding an approximate value for the summation sum(An*cos(pny)exp(-k*pn^2*t)
for n=1:N
    Pn = ((n*pi)/h);
    An = ((2*h)/(n^2*pi^2))*((-1^n*Q2)-Q1);
    sumComp = sumComp + An*cos(Pn*y)*exp(-K*Pn^2*t);
end

% returning solution to u
u = WY + GT + A0 + sumComp;

end




% Global variable functions to make the consts easier to work with
% they make MATLAB angry but make me happy so MATLAB can cry about it
 
function setGlobalKappa(val)
global kp
kp = val;
end

function r = getKappa
global kp
r = kp;
end

function setGlobalH(val)
global h
h = val;
end

function r = getH
global h
r = h;
end

function setGlobalTi(val)
global Ti;
Ti = val;
end

function r = getTi
global Ti
r = Ti;
end

function setGlobalQ1(val)
global Q1
Q1 = val;
end

function r = getQ1
global Q1
r = Q1;
end