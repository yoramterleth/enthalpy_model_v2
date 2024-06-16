function dydt = ode_MAIN(t,y,C,P,tS,TaS, g_p)


        %% Arpil 2024 // Yoram Terleth

        % this function includes the system of PDEs. 
        
        
        % output: y output values of the 4 PDEs
        % dydt(:,1)= ice thickness
        % dydt(:,2)= basal enthalpy
        % dydt(:,3)= channel size 
        % dydt(:,4)= subglacial storage space (not included in model submitted
        % for review)

        % input: 
        % - t: time vector
        % - y: ealier values
        % - C constants structure
        % - P parameters sturcture
        % - tS initial time vector (for seasonality)
        % - TaS: air temperatures with seasonality
        % - minimum fraction of surface melt that is allowed to the
        % distiribued draonage system. 
        


    %%  make array 
    dydt = zeros(4,1) ;

    
    %% define functions

    % get the average annual temperature 
    annualTa = mean(TaS) ; 

    % this time with seasonality 
    Ta = interp1(tS, TaS,t,'linear') ; 

    % melt %% !! C.DDF2 implements the concentration of all the melt to
    % only when air temperature is >0 deg C! 
    m = max(0, C.DDF2 * (Ta * C.T_0 - C.Tm) / C.a_0) ; 

    % Enthalpy plus
    Eplus = max(y(2)*C.E0,0)/C.E0 ; 


    %enthalpy minus 
    Eminus = min(y(2)*C.E0, 0)/C.E0 ; 


    % E that is in storage: any enthalpy goes into storage first (see eq
    % lower down for y(4) description. 
    Estorage = min(y(4), Eplus) ;

    % remaining E, that is not stored and can be drained
    Edrainable = max(Eplus - Estorage,0) ; 


    % N - effective pressure
    N = min(y(1)/C.chi, 1/(Edrainable)) ;

    % u - sliding velocity
    u = P.slope^(1/C.p) * y(1).^(1+(1/C.p)) .* N.^(-C.q/C.p) ; 

    % Beta - protion of surface melt that enters the distributed subglacial
    % drainage system. 
    Beta = min(max(g_p,((u*C.u_0 - P.u1)/(P.u2-P.u1))),1) ; 

    % Phi - fill fraction
    Phi = min(1,(Eplus/(y(1)/C.chi))) ;
    

    % H - describing change in ice thickness. 
    dydt(1) = P.a - m - (1/P.l)* (P.slope^(1/C.p) * y(1)^(1+(1/C.p)) * N^(-C.q/C.p) + (C.lambda * (P.slope^(C.n)))) ; 

    % E - describing change in basal enthalpy
    dydt(2) = ((P.slope^(1+(1/C.p)) * y(1)^(1+(1/C.p)) * N^(-C.q/C.p) + C.gamma - C.kappa * ((Eminus - annualTa)/y(1)) - (1/P.l)*(P.slope * Eplus^(C.alpha) + (Phi * P.slope^(1/2) *y(3)^(4/3))) + (C.delta*Beta*m))/C.mu) ;
    

    % S - describing change in subglacial channel size. 
    dydt(3) = ((C.sigma * Phi * P.slope^(3/2) * y(3)^(4/3)) - (y(3) * N^(C.n)) + C.S_0_hat)/ C.v ; 
    
   
    % Ac !Not included in the model described in Terleth et al. 24.! This variable is based on Bartholomaus et al. 2011. and
    % describes change subglacial storage space. It is disabled by setting
    % p.h_cavity to 0. It describes transient subglacial storage space.
    % That exists below the glacier but does not contribute to effective
    % pressure drops.
    dydt(4) = ((u * P.h_cavity) + ((1 - C.rothl) * ((C.p_water*C.g)/(C.p_ice* C.Lf))*(C.cav_frac * ((1/P.l)*(P.slope * Estorage^(C.alpha))))*(P.slope/C.cavity_sinuosity)) - ( y(4) * ((N/C.n)^C.n)) ) / C.v ;  
    


end 
