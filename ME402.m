clear all;
close all;
clc;

k = 0.251;
rho = 1262;
cps = 1339;
cpl = 1339;
Ts = 603;
Tl = 618;
hs = cps*Ts;
hl = cpl*Tl;

Tinitial = 640;
Twall = 275;

x = linspace(0,200,20000);
y = linspace(0,100,10000);

xinitial = 0;
delta_x = 0.0001;

yinitial = 0;
delta_y = 0.0001;

delta_t = 0.1;
i = 1; j = 1; k = 1;
T = zeros(length(x),length(y),k);
Tr = Twall*ones(length(x),length(y));
for i = 1:size(x)
    for j = 1:size(y)
        if i==1 || i==n || j==1 || j==m
            T(i,j,k) = Twall;
        else
            T(i,j,k) = Tinitial;
        end
    end
end


while T(i,j,k) - Tr <= 0.000005
k = k+1;
for i = 1:size(x)
    for j = 1:size(y)
        if i==1 || i==n || j==1 || j==m
            T(i,j,k) = Twall;
        else
            if T(i,j,k-1) > Tl
                T(i,j,k) = T(i,j,k-1) + (T(i+1,j,k-1) + T(i-1,j,k-1) + T(i,j+1,k-1) + T(i,j-1,k-1))*...
                    (delta_t/cpl);
            else
                if T(i,j) <= Tl || T(i,j) >= Ts
                T(i,j,k) = T(i,j,k-1) + (T(i+1,j,k-1) + T(i-1,j,k-1) + T(i,j+1,k-1) + T(i,j-1,k-1))*...
                    (delta_t*(Tl-Ts)/(hl-hs));
                else 
                    if T(i,j) < Ts
                T(i,j,k) = T(i,j,k-1) + (T(i+1,j,k-1) + T(i-1,j,k-1) + T(i,j+1,k-1) + T(i,j-1,k-1))*...
                    (delta_t/cps);
                    end
                end
            end
        end
    end
end
end


            
                
            
            