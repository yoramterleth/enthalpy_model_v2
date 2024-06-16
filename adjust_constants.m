function [C,P] = adjust_constants(C,P,K)
    
    %% this function recalculates all the "constants". or non-free parameters, that change as a result of 
    %% changing K and epsilon to modify the drainage system 
    
    % Yoram Terleth - summer 2023. 

    % completed and tested on 18 Aug 2023. Reproduces the results obtained
    % by using Benn et al values directly, as well as non-rounded versions
    % of all the parameters in Table 1. Changing K affects the phase
    % diagram. 


    %%  convert the silly units 
    C.a_0_ms = C.a_0 / 3.154e7 ; % m a to m s 
    

    %% inital quantities 

    % E0 
    C.E0 = ((C.g * C.sin_theta * C.a_0_ms * (C.l_0^2)) / (C.Lf * K))^(1/C.alpha) ; 
    % C.E0 = 1.8e8 ; 
    % T0
    C.T_0 = C.E0 / (C.p_ice *  C.Cp * C.d) ; 

    % w0
    C.W_0 = C.E0 / (C.p_ice*C.Lf) ; 

    % N0 
    C.N_0 = C.C / C.E0 ; 

    % H0
    C.H0 = ((C.R * (C.C^C.q) * (C.a_0_ms^C.p) *(C.l_0^C.p)) / (C.p_ice * C.g * C.sin_theta * (C.E0^C.q)))^(1/(C.p+1)) ; 

    % t0
    % C.t_0 = C.H0 / C.a_0 ; 

    % u0
    C.u_0_ms = ((C.p_ice * C.g * C.sin_theta * (C.E0^C.q) * C.a_0_ms * C.l_0) / (C.R * (C.C^C.q)))^(1/(C.p+1)) ; 
    C.u_0 = C.u_0_ms / 3.168874e-8 ; 

    % Q0
    C.Q_0 = (C.g * C.sin_theta * C.a_0_ms * (C.l_0^2)) / C.Lf ; 

    % S0 
    C.S_0 = ((C.g * C.sin_theta * C.a_0_ms * (C.l_0^2) * C.Wc) / (C.Lf * C.Kc *  ((C.p_ice * C.g * C.sin_theta)^(1/2))))^(3/4) ; 
    

    %% other dimensionless parameters
    
    C.tau_0 = 0.0029 ; % not in paper but worked back from other constants using (C.G / C.gamma) / C.u_0 

    % gamma 
    C.gamma = C.G / (C.tau_0 * C.u_0) ; 

    % kappa 
    C.kappa = (C.k * C.T_0) / (C.tau_0 * C.u_0 * C.H0) ; 

    % delta 
    C.delta = (C.p_ice * C.Lf * C.a_0_ms) / (C.tau_0 * C.u_0) ; 

    % mu 
    C.mu = (C.E0 * C.a_0_ms) / (C.tau_0 * C.u_0 * C.H0) ; 

    % chi 
    C.chi = C.N_0 / (C.p_ice * C.g * C.H0) ;

    % lambda 
    C.lambda = ((2 * C.A) * ((C.p_ice * C.g * C.sin_theta)^C.n) * (C.H0^(C.n+1))) / ((C.n+2) * C.u_0_ms) ; 

    % v 
    C.t_0_s = C.t_0 * 3.154e+7 ; 
    C.v = 1 / (C.t_0_s * C.A_t * C.N_0^C.n) ; 

    % sigma 
    C.sigma = (C.Kc * ((C.p_ice * C.g * C.sin_theta)^(3/2)) * (C.S_0^(1/3))) / (C.p_ice * C.Lf * C.A_t * (C.N_0^C.n)) ; 

    % S0 hat 
    C.S_0_hat = C.So / (C.S_0 * C.A_t * (C.N_0^C.n)) ; 

    
    % add the epsilon value to the constants structure
    % C.epsilon = epsilon ;  

end