% ME 550 -- Winter 2021 -- Project 2
% Austin Warren and Genevieve Gaudin
clear; clc; close all;


% finite diff method to determine time and heat stored for brick of PCM
% melting

%% domain set up
width = 1*0.0254; % width in meters
pcm_length = 2*0.0254; % length in meters
height = 1*0.0254; % height of a brick in meters
N = 15/height; % number of bricks in a stack
dx = 0.001; % delta x
x = 0:dx:width; % x length in meters
y = 0:dx:pcm_length; % y length in meters



%% properties
% air properties
rho_air = 0.4289; % air density at 550 C in kg/m^3
T_air = 550 + 273; % temp of air in K
mu_air = 377e-7; % air viscosity at 825 K in Pa s
cp_air = 1.1045*1000; % air specific heat in J/kgK
k_air = 58.45e-3; % air thermal conductivity in W/mK
U_air = 10; % air velocity in m/s
Re_air = rho_air * U_air * height / mu_air; % reynolds number of air
Pr_air = 0.7125; % air prandtl number

% PCM properties
rho_pcm = [2320, 2480, 2044]; % pcm density in kg/m^3
k_pcm = [2.09, 0.6, 0.5]; % pcm thermal conductivity in W/mK
cp_pcm = [2360, 1450, 1470]*1000; % pcm specific heat in J/kgK
T_melt = [496, 435, 380]+273; % pcm melting temp in K
L_pcm = [370, 350, 150]*1000; % pcm latent heat of fusion in J/kg
alpha_pcm = k_pcm./(rho_pcm.*cp_pcm); % pcm alpha in m^2/s
Ti = 500; % initial temp in K

h_pcm = k_pcm./(height.*N) .* (2 + 1.1.*(Re_air.^(0.6).*Pr_air.^(1/3))); % convective heat transfer coefficient in W/m^2K
Bi_pcm = h_pcm .* dx ./ k_pcm; % pcm biot number

% find dt
[dt, dt_all] = find_dt_please_2(dx, alpha_pcm, Bi_pcm);
%dt = 3600; % manually set dt
% Fourier number
Fo_pcm = k_pcm.*dt./(rho_pcm.*cp_pcm.*dx.^2);




%% loops!!
iter = [5e5, 1e5, 1e5];
t_melt = [0,0,0];

% outer loop for pcm
for k = 1:length(rho_pcm)
    
    % initialize the arrays for temp, total heat, and heat of fusion
    T = zeros(length(x),length(y))+ Ti;
    q = zeros(length(x), length(y));
    q_pc = zeros(length(x),length(y));
    
    % set Fo, Bi, h,
    Fo = Fo_pcm(k);
    Bi = Bi_pcm(k);
    h = h_pcm(k);
    
    
    % loop for time -- iteration limit and test for all nodes melted
    for i = 1:iter(k)
        
        
        % x loop
        for m = 1:length(x)
            
            
            
            % y loop
            for n = 1:length(y)
                
                % if temp is less than melting temp --> normal finite diff
                if T(m,n) < T_melt(k)
                    % add properties for phase change? currently using
                    % average for everything
                    
                    % boundary conditions
                    if m==1 && n==1 % bottom left corner
                        T(m,n) = (1-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m+1,n)+T(m,n+1))) + (2*Fo*Bi*T_air);
                    elseif m==1 && n==length(y) % top left corner
                        T(m,n) = (1-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m+1,n)+T(m,n-1))) + (2*Fo*Bi*T_air);
                    elseif m==length(x) && n==1 % bottom right corner
                        T(m,n) = (1-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m-1,n)+T(m,n+1))) + (2*Fo*Bi*T_air);
                    elseif m==length(x) && n==length(y) % top right corner
                        T(m,n) = (1-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m-1,n)+T(m,n-1))) + (2*Fo*Bi*T_air);
                    elseif m==1 && n<length(y) && n>1 % left side
                        T(m,n) = (1-(2*Fo))*T(m,n) + (0.5*Fo*(2*T(m+1,n) + T(m,n-1) + T(m,n+1)));
                    elseif m==length(x) && n<length(y) && n>1 % right side
                        T(m,n) = (1-(2*Fo))*T(m,n) + (0.5*Fo*(2*T(m-1,n) + T(m,n-1) + T(m,n+1)));
                    elseif m<length(x) && m>1 && n==1 % bottom
                        T(m,n) = (1- 2*Fo - Bi*Fo)*T(m,n) + (0.5*Fo*(2*T(m,n+1) + T(m+1,n) + T(m-1,n))) + (Bi*Fo*T_air);
                    elseif m<length(x) && m>1 && n==length(y) % top
                        T(m,n) = (1- 2*Fo - Bi*Fo)*T(m,n) + (0.5*Fo*(2*T(m,n-1) + T(m+1,n) + T(m-1,n))) + (Bi*Fo*T_air);
                    elseif m<length(x) && m>1 && n<length(y) && n>1 % interior nodes
                        T(m,n) = (1- 4*Fo)*T(m,n) + Fo*(T(m+1,n) + T(m-1,n) + T(m,n+1) + T(m,n-1));
                    else
                        fprintf('You messed up the indexing, dummy\n')
                    end
                end
                
                % if temp is above melting temp and the latent heat is not
                % to the total for fusion, convert additional temp to heat
                if T(m,n) > T_melt(k) && q_pc(m,n) ==0
                    q_pc(m,n) = cp_pcm(k) * (T(m,n)-T_melt(k)); % calc latent heat
                    if q_pc(m,n) <= L_pcm(k)
                        T(m,n) = T_melt(k); % set temp of node to the melting temp if not phase changed
                    elseif q_pc(m,n) > L_pcm(k)
                        T(m,n) = T_melt(k) + (q_pc(m,n) - L_pcm(k))/cp_pcm(k); % if latent heat is enough to fully phase change, convert additional heat to temp
                        q_pc(m,n) = L_pcm(k);
                    end
                elseif T(m,n) == T_melt(k) && q_pc(m,n)<=L_pcm(k) && q_pc(m,n)>=0
                    % during phase change -- calc heat instead of temp
                    % using finite diff temp eqns
                    % boundary conditions
                    if m==1 && n==1 % bottom left corner
                        q_pc(m,n) = q_pc(m,n) + cp_pcm(k)*((-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m+1,n)+T(m,n+1))) + (2*Fo*Bi*T_air));
                    elseif m==1 && n==length(y) % top left corner
                        q_pc(m,n) = q_pc(m,n) + cp_pcm(k)*((-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m+1,n)+T(m,n-1))) + (2*Fo*Bi*T_air));
                    elseif m==length(x) && n==1 % bottom right corner
                        q_pc(m,n) = q_pc(m,n) + cp_pcm(k)*((-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m-1,n)+T(m,n+1))) + (2*Fo*Bi*T_air));
                    elseif m==length(x) && n==length(y) % top right corner
                        q_pc(m,n) = q_pc(m,n) + cp_pcm(k)*((-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m-1,n)+T(m,n-1))) + (2*Fo*Bi*T_air));
                    elseif m==1 && n<length(y) && n>1 % left side
                       q_pc(m,n) = q_pc(m,n) + cp_pcm(k)*((-(2*Fo))*T(m,n) + (0.5*Fo*(2*T(m+1,n) + T(m,n-1) + T(m,n+1))));
                    elseif m==length(x) && n<length(y) && n>1 % right side
                        q_pc(m,n) = q_pc(m,n) + cp_pcm(k)*((-(2*Fo))*T(m,n) + (0.5*Fo*(2*T(m-1,n) + T(m,n-1) + T(m,n+1))));
                    elseif m<length(x) && m>1 && n==1 % bottom
                        q_pc(m,n) = q_pc(m,n) + cp_pcm(k)*((-2*Fo - Bi*Fo)*T(m,n) + (0.5*Fo*(2*T(m,n+1) + T(m+1,n) + T(m-1,n))) + (Bi*Fo*T_air));
                    elseif m<length(x) && m>1 && n==length(y) % top
                        q_pc(m,n) = q_pc(m,n) + cp_pcm(k)*((-2*Fo - Bi*Fo)*T(m,n) + (0.5*Fo*(2*T(m,n-1) + T(m+1,n) + T(m-1,n))) + (Bi*Fo*T_air));
                    elseif m<length(x) && m>1 && n<length(y) && n>1 % interior nodes
                        q_pc(m,n) = q_pc(m,n) + cp_pcm(k)*((-4*Fo)*T(m,n) + Fo*(T(m+1,n) + T(m-1,n) + T(m,n+1) + T(m,n-1)));
                    else
                        fprintf('You messed up the indexing, dummy\n')
                    end
                    
                    % double check temp and phase change
                    if q_pc(m,n) <= L_pcm(k)
                        T(m,n) = T_melt(k); % set temp of node to the melting temp if not phase changed
                    elseif q_pc(m,n) > L_pcm(k)
                        T(m,n) = T_melt(k) + (q_pc(m,n) - L_pcm(k))/cp_pcm(k); % if latent heat is enough to fully phase change, convert additional heat to temp
                        q_pc(m,n) = L_pcm(k);
                    end
                    
                elseif T(m,n) == T_melt(k) && q_pc(m,n)>L_pcm(k)
                    T(m,n) = T_melt(k) + (q_pc(m,n) - L_pcm(k))/cp_pcm(k);
                    
                elseif T(m,n)>T_melt(k) && q_pc(m,n)==L_pcm(k)
                    % normal finite diff for after phase change -- add
                    % liquid properties if desired
                    
                    % boundary conditions
                    if m==1 && n==1 % bottom left corner
                        T(m,n) = (1-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m+1,n)+T(m,n+1))) + (2*Fo*Bi*T_air);
                    elseif m==1 && n==length(y) % top left corner
                        T(m,n) = (1-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m+1,n)+T(m,n-1))) + (2*Fo*Bi*T_air);
                    elseif m==length(x) && n==1 % bottom right corner
                        T(m,n) = (1-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m-1,n)+T(m,n+1))) + (2*Fo*Bi*T_air);
                    elseif m==length(x) && n==length(y) % top right corner
                        T(m,n) = (1-(4*Fo)-(2*Fo*Bi))*T(m,n) + (2*Fo*(T(m-1,n)+T(m,n-1))) + (2*Fo*Bi*T_air);
                    elseif m==1 && n<length(y) && n>1 % left side
                        T(m,n) = (1-(2*Fo))*T(m,n) + (0.5*Fo*(2*T(m+1,n) + T(m,n-1) + T(m,n+1)));
                    elseif m==length(x) && n<length(y) && n>1 % right side
                        T(m,n) = (1-(2*Fo))*T(m,n) + (0.5*Fo*(2*T(m-1,n) + T(m,n-1) + T(m,n+1)));
                    elseif m<length(x) && m>1 && n==1 % bottom
                        T(m,n) = (1- 2*Fo - Bi*Fo)*T(m,n) + (0.5*Fo*(2*T(m,n+1) + T(m+1,n) + T(m-1,n))) + (Bi*Fo*T_air);
                    elseif m<length(x) && m>1 && n==length(y) % top
                        T(m,n) = (1- 2*Fo - Bi*Fo)*T(m,n) + (0.5*Fo*(2*T(m,n-1) + T(m+1,n) + T(m-1,n))) + (Bi*Fo*T_air);
                    elseif m<length(x) && m>1 && n<length(y) && n>1 % interior nodes
                        T(m,n) = (1- 4*Fo)*T(m,n) + Fo*(T(m+1,n) + T(m-1,n) + T(m,n+1) + T(m,n-1));
                    else
                        fprintf('You messed up the indexing, dummy\n')
                    end
                    
                elseif T(m,n)>T_melt(k) && q_pc(m,n)<L_pcm(k)
                    fprintf('T>Tmelt and q_pc<L\n')
                elseif T(m,n)>T_melt(k) && q_pc(m,n)>L_pcm(k)
                    fprintf('T>Tmelt and q_pc>L\n')
                %else
                    %fprintf('Something went wrong with phase change calcs\n')
                    
                end
                
                
                
                
                % if temp is above melting temp and the latent heat is
                % total required for melting --> node is melted: use liquid
                % phase properties for the finite diff
                        
    
    

    
    
    
    
            end
            
            
        end
    
    
    
    
        % break condition for all nodes melting
        if min(q_pc,[],'all') >= L_pcm(k)
            t_melt(k) = dt(k)*i;
            fprintf('Time to melt:%f s\n',t_melt(k))
            disp(i)
            break
        end
    
    
    end
    
    
    
    
    
    
    
    
    
    
    
end
