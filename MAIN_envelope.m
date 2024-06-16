
%% MAIN envelope runs trhough the amount of water allowed through to the base.

% - cavity size is fixed.
% - annual mean air temp is fixed.
% - annual accumulation is fixed.
% - we vary K the basal hydraulic transmissitvity. 
% - we vary the amount of water let through to the base. 

% Terleth, April 2024.

%% 
close all
clearvars
clc

% run counter for progress report
counter = 1 ; 

%% PATH TO SAVE TO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_path = 'D:\APRIL_envelope\' ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% override existing files of the same name?
overide_output = 1 ; 

%% load constants 

C = constants() ; 


%% define vectors

% proportions of water let through to the base.
p_array = sort([0:0.025:0.4,0.5,0.75,1]) ; 

% values of [times] k 
k_array =  [0.001,0.0025,0.003,0.004, 0.005,0.006,0.0075,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.75,1,5, 10] ; 


%% define other parameters
 
a = 0.6 ; 
Ta = -8; 

l = 10e3 ; 
slope = 1 ; 

seasonal_amp = 15 ; % deg C

P.u1 = 10 ;         % m a-1 
P.u2 = 100 ;        % m a-1

%% set cavity step height 
P.h_cavity = 0 ; 

%% normalise other parameters 
P.l = l / C.l_0 ; 
P.slope = slope / 1 ; 

%% ODE solver parameters  

% intial conditions (H/Ho and E/Eo and S/So)
init = [1,1,1,0] ; 

% timespan (t1/t0 tend/t0) 
time = 0:1e-4:60 ; % 

% solver options 
options = odeset('RelTol',1e-6,'Stats','off','OutputFcn', @odewbar) ; 



%% iterate over options 
for j = 1:length(p_array) 
    disp(['Glacier macro-permeability set to ' num2str(p_array(j)*100) '%.'] )

    for i = 1:length(k_array)
        
        %% K ; 
        K = C.K * k_array(i) ; 
        disp(['Basal hydraulic conductivity set to ' num2str(k_array(i)) ' times K_0.'])

        %% check if exists 
        if ~overide_output && isfile([save_path 'params_gp_' num2str(p_array(j)) '_kp_' num2str(k_array(i)) '.mat'])
            disp('File already exists, skipping.')
            counter = counter+1 ; 
            continue 
        end 

        %% Apply seasonality to Ta  
        [TaS,TaConstant] = make_seasonality(time*C.t_0,Ta,seasonal_amp) ; 
        timeS = time ; 

        %% normalise parameters 
        P.a = a / C.a_0 ; 
        TaSnorm = TaS / C.T_0 ; 
        P.l = l / C.l_0 ; 
        P.slope = slope / 1 ; 
        
        % compute melt 
        m = C.DDF * (Ta * C.T_0 - C.Toffset) / C.a_0 ; 
        P.m = max(m,0) ;


        %% compute the adjusted DDF
        T1 = max(TaConstant - C.Toffset,0) ; 
        T2 = max(TaSnorm .* C.T_0 - C.Tm,0) ;  
        DDF2 = C.DDF .* sum(T1)/sum(T2) ; 
        
        C.DDF2 = DDF2 ; 

        
        % adjust constants 
        [C_mod,P] = adjust_constants(C,P,K) ; 

        %% set the value of glacier permeability
        g_p = p_array(j) ; 

        %% solve ODE 
        disp('Solving ODE...')
        [t,y] = ode23s(@(t,y)ode_MAIN(t,y,C_mod,P,timeS,TaSnorm,g_p),time,init, options)  ;  
        
        %% save the associated workspace
        disp('Saving...')

        if exist(save_path, 'dir') ~= 7
            % Folder doesn't exist, create it
            mkdir(save_path);
            disp(['Folder "', save_path, '" created.']);
        else
            % Folder already exists
            disp(['Folder "', save_path, '" already exists.']);
        end
        
        
        save([save_path 'params_gp_' num2str(p_array(j)) '_kp_' num2str(k_array(i)) '.mat'])
        
        disp(['Done with r un nb.', num2str(counter), ' of ', num2str(length(p_array)*length(k_array))])
        counter = counter+1 ; 
    end 
end 

%% save and end 

save([save_path 'FINAL_workspace.mat'])


