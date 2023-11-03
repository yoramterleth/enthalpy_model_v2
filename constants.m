function [C] = constants() 

% constants listed in table 1 in Benn et al.(2019)

C.p_ice = 916 ;     % kg m-3
C.g = 10 ;        % m s-2

C.sin_theta = 0.05 ;    % sin(bed slope)

C.Lf = 3.3e5 ;      % J kg-1

C.Cp = 2000 ;       % K kg-1 K-1
C.k =  2.1 ;        % W m-1 K-1
C.G = 0.06 ;        %  W m-2
C.d = 10 ;          % m 

C.n = 3 ;           % 
C.A = 2.4e-25 ;     % Pa-3 s-1
C.p = 1/3 ; 
C.q = 1 ; 
C.R = 15.7 ;        % m-1/3 s1/3
C.alpha = 5 ;  % 2 % 5 ;       % 

C.K = 2.3e-47 ;     % kg-5 m2 s9
C.C = 9.2e13 ;      % Pa J m-2 

C.DDF = 0.1 ;         % m a-1 K-1 pr year! 
C.Tm =  0 ; % 273.15 ;      % K 
C.Toffset =  -10 ; % 263.15 ; % K 

C.Kc = 0.04 ;       % m4/3 kg1/2
C.Wc = 1000 ;       % m 

C.A_t = 1.8e-25 ;   % Pa-3 s-1
C.So = 3e-13 ;      % m2 s-2
C.a_0 = 1 ;         % m a -1 
C.l_0 = 10000 ;     % m 
C.E0 = 1.8e8 ;     % J m-2
C.T_0 = 10 ; % - 263.15 ;  % 283.15 ;        % K 
C.W_0 = 0.6 ;       % m 
C.N_0 = 0.5e6 ;     % Pa 
C.t_0 = 5 ;       % a (years) 
C.H0 = 200 ;       % m 
C.u_0 = 50 ;        % m a -1 
C.Q_0 = 5e-6 ;      % m2 s-1 
C.S_0 = 0.02 ;      % m2 

C.gamma = 0.41 ;  %  
C.kappa = 0.7 ;     % 
C.delta = 66 ; 
C.mu = 0.2 ;
C.chi = 0.27 ; 
C.lambda = 0.009 ;
C.v = 0.007 ;
C.sigma = 16  ; %7 ; 
C.S_0_hat = 0.0007 ; 


%%% aditional constants that are needed to add equation 4 from Bartholomaus
%%% et al, describing cavity size 
C.h_cavity = 1 ; % m 
C.rothl = 0.313 ; % rothlisberger dimensionless constant 
C.cavity_sinuosity = 3 ; % following Bartholomaus
C.p_water = 1000 ; % water density in kg m3
C.cav_frac = 0.25 ; % space of bed occupied by cavities


end 
