function dydt = ode_CASE3_seasonality(t,y,C,P,tS,TaS)

    % this is case 3 in Benn et al. 2019: it includes transfer of surface
    % melt to the glacier bed and the consideration of an adaptive drainage
    % system. 

        % this case is a a lot stiffer than the others! Use ode23s to
        % solve... 

        % now using seasonality 

    %%  make array 
    dydt = zeros(3,1) ;

    
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

    % N 
    N = min(y(1)/C.chi, 1/(Eplus)) ;

    % u 
    u = P.slope^(1/C.p) * y(1).^(1+(1/C.p)) .* N.^(-C.q/C.p) ; 

    % Beta 
    %Beta = min(max(0,((u*C.u_0 - P.u1)/(P.u2-P.u1))),1) ; 

    % changed so that it always lets through some of the melt even at zero
    % speed: 
    Beta = min(max(0.1,((u*C.u_0 - P.u1)/(P.u2-P.u1))),1) ; 

    % Phi 
    Phi = min(1,(Eplus/(y(1)/C.chi))) ; 

    % H %% added this back to P.m so that it can go back to annual melt rate! 
    dydt(1) = P.a - m - (1/P.l)* (P.slope^(1/C.p) * y(1)^(1+(1/C.p)) * N^(-C.q/C.p) + (C.lambda * (P.slope^(C.n)))) ; 

    % E 
    dydt(2) = ((P.slope^(1+(1/C.p)) * y(1)^(1+(1/C.p)) * N^(-C.q/C.p) + C.gamma - C.kappa * ((Eminus - annualTa)/y(1)) - (1/P.l)*(P.slope * Eplus^(C.alpha) + (Phi * P.slope^(1/2) *y(3)^(4/3))) + (C.delta*Beta*m))/C.mu) ;

    % S
    dydt(3) = ((C.sigma * Phi * P.slope^(3/2) * y(3)^(4/3)) - (y(3) * N^(C.n)) + C.S_0_hat)/ C.v ; 

end 
