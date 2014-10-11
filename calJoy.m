function Flux = calJoy(omega,thickness_list, epsilon_list,kx_list,sourceLayer,targetLayer,isLossy)
%This function calculate the heat flux for source joy
%   omega:          fequency, scalar 
%   thickness_list: thickness from layer 1 to N, should start and end with
%                   Inf, vector
%   epsilon_list:   complex vector, containing epsilon for every layer
%   sourceLayer:    at which layer source is placed, should be a lossy layer
%   targetLayer:    at which layer the flux is calculated, should be a
%                   lossless layer
%   isLossy:        a vector containing 1 and 0 indicate whether a layer is
%                   lossy

numOfLayer=size(epsilon_list,1);
Flux=zeros(size(kx_list));
kz = zeros(numOfLayer,1);
f_list=zeros(numOfLayer,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kx_index=length(kx_list)
    %% initializing
    kx = kx_list(kx_index);
    for i=1:numOfLayer
        kz(i)=sqrt(omega^2/c^2*epsilon_list(i)-kx^2);
        if imag(kz(i))>=0
            kz(i)=-kz(i);
        end
        if thickness_list(i)==Inf
            f_list(i)=1;
        else
            f_list(i)=exp(-1i*kz(i)*thickness_list(i));
        end
    end
    %% if source is placed at the bottom or top layer, then put it exactly at surface
    if sourceLayer==1 || sourceLayer==numOfLayer
        [ e,h ] = pointJoy(omega,thickness_list, epsilon_list,kx,sourceLayer,0,isLossy);
        Flux(kx_index)=-e(targetLayer)*conj(h(targetLayer))*(-1/(2*imag(kz(sourceLayer))));
    %% if source is placed at other layers, then put it exactly at top and bottom
    else
        [e_up,h_up] = pointJoy(omega,thickness_list, epsilon_list,kx,sourceLayer,thickness_list(sourceLayer),isLossy);
        [e_down,h_down] = pointJoy(omega,thickness_list, epsilon_list,kx,sourceLayer,0,isLossy);
        A=[f_list(sourceLayer),1;1,f_list(sourceLayer)]\[e_up(targetLayer);e_down(targetLayer)];
        B=[f_list(sourceLayer),1;1,f_list(sourceLayer)]\[h_up(targetLayer);h_down(targetLayer)];
        Flux(kx_index)=-(A(1)*conj(B(1))+A(2)*conj(B(2)))*(exp(2*imag(kz(sourceLayer))*thickness_list(sourceLayer))-1)/(2*imag(kz(sourceLayer)))-...
            A(2)*conj(B(1))*exp(-1i*kz(sourceLayer)*thickness_list(sourceLayer))*(exp(2*1i*kz(sourceLayer)*thickness_list(sourceLayer))-1)/(2*1i*real(kz(sourceLayer)))-...
            A(1)*conj(B(2))*exp(-1i*conj(kz(sourceLayer))*thickness_list(sourceLayer))*(1-exp(-2*1i*kz(sourceLayer)*thickness_list(sourceLayer)))/(2*1i*real(kz(sourceLayer)));
    end
    %% finally obtain the flux
    Flux(kx_index)=-omega*imag(epsilon_list(sourceLayer)*ep0)/pi^2*real(Flux(kx_index));
end

end

