%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name : Sagar Jagina
% ME 402 Assignment 2: Cure Kinetics Project-2
% Problem: To generate a coputer code to simulate Cure-Kinetics
% situation during a resin curing process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
% Material is PEEK(polyetheretherketone)
k = 11.4; % conductivity
rho = 1100; % Density
cp = 935; % Specific heat in solid state
R = 8.314;
Hr = 479;

B = 0.054;
A1 = 2.101e09;
A2 = -2.014e09;
A3 = 1.960e05;
E1 = 8.07e04;
E2 = 7.78e04;
E3 = 5.66e04;

Tinitial = 300; % Initial temperature of the molten material
Twall = 420; % Temperature of the casting mold

% Making nodes in X and Y direction
x = linspace(0,2,2000);
y = linspace(0,0.1,100);
z = linspace(0,0.005,5);

delta_x = 0.001;
delta_y = 0.001;
delta_z = 0.001;


% Time step
delta_t = 0.001;

% Satbility Factor
r = k * delta_t/(rho*cp*delta_x^2);

C1 = (k*delta_t)/(rho*cp);
C2 = (Hr*delta_t)/cp;

% Initial temperature condition matrix
for i = 1:length(y)
    for j = 1:length(x)
        for k = 1:length(z)
        if i==1 || i==length(y) || j==1 || j==length(x) || k==1 || k==length(z)
            T(i,j,k) = Twall;
        else
            T(i,j,k) = Tinitial;
        end
        end
    end
end


alpha = zeros(length(y),length(x),length(z));
for i = 1:length(y)
    for j = 1:length(x)
        for k = 1:length(z)
            alpha(i,j,k) = 0.05;
        end
    end
end


Tnew  = zeros(length(y),length(x),length(z));
alpha_new = zeros(length(y),length(x),length(z));

n = 0;
MIN = 0;
ALPHA = 0;
while 0.981 -  ALPHA >= 5e-03

    n = n + 1;
        Tm(n) = T(50,1000,3);    
        alpham(n) = alpha(50,1000,3);   
        ALPHA = alpham(n);
        
for k = 1:length(z)
    for j = 1:length(x)
        for i = 1:length(y)
            
             if i==1 || i==length(y) || j==1 || j==length(x) || k==1 || k==length(z)
                Tnew(i,j,k) = Twall;

                if  le(alpha(i,j,k),0.3)
                    dalpha_dt = (A1*exp(-E1/(R*T(i,j,k))) - A2*exp(-E2/(R*T(i,j,k)))*alpha(i,j,k)) * (1 - alpha(i,j,k)) * (B - alpha(i,j,k));
                elseif alpha(i,j,k) > 0.3
                    dalpha_dt = A3*exp(-E3/(R*T(i,j,k))) * (1 - alpha(i,j,k));
                end
                        alpha_new(i,j,k) = alpha(i,j,k) + (dalpha_dt* delta_t);

             else
            
                if le(alpha(i,j,k),0.3)
                    dalpha_dt = (A1*exp(-E1/(R*T(i,j,k))) - A2*exp(-E2/(R*T(i,j,k)))*alpha(i,j,k)) * (1 - alpha(i,j,k)) * (B - alpha(i,j,k));
                elseif alpha(i,j,k) > 0.3
                    dalpha_dt = A3*exp(-E3/(R*T(i,j,k))) * (1 - alpha(i,j,k));
                end
                        alpha_new(i,j,k) = alpha(i,j,k) + (dalpha_dt* delta_t);
                        
                 Tnew(i,j,k) = T(i,j,k) + (C1/delta_z^2)*(T(i+1,j,k) + T(i-1,j,k) + T(i,j+1,k) + T(i,j-1,k) + T(i,j,k+1) + T(i,j,k-1) - 6*T(i,j,k)) + (C2 *  dalpha_dt);

                
             end
        end
    end
end
T = Tnew;
Tnew = zeros(length(y),length(x),length(z));

alpha = alpha_new;
alpha_new = zeros(length(y),length(x),length(z));
tt1 = alpha(:,:,2);
tt2 = T(:,:,3);
end



