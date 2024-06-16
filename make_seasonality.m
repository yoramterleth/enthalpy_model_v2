function [TaS,TaConstant] = make_seasonality(t,Ta,seasonal_amp) 

% this function makes an array of seasonally variable air temperatures. 

TaS = zeros(length(t),1); 
TaConstant = zeros(length(t),1); 

for i = 1:length(t)
TaS(i) = Ta +  (seasonal_amp * sin(2*pi*t(i)));
TaConstant(i) = Ta ; 
end 
end
