clear;
%% open matlabpool in Hera
%% matlabpool open 4;

%% initializing parameters
materialNames={'InAs';'Air'};
thickness_list=[Inf;Inf];   % thickness should start and end with Inf
isLossy=[1;0];
numOfLayer=size(thickness_list,1);
sourceLayer=1;   % make sure the source layer is a lossy layer
targetLayer=2;
d=-1e-7;         % make sure the source is inside the layer
kx=0;

epsilon=cell(numOfLayer,1);
for i=1:numOfLayer
    [omega_list,epsilon{i}]=addMaterial(materialNames{i});
end

exx=cell(size(omega_list,1),1);
hyx=cell(size(omega_list,1),1);
eyy=cell(size(omega_list,1),1);
hxy=cell(size(omega_list,1),1);
exz=cell(size(omega_list,1),1);
hyz=cell(size(omega_list,1),1);
fluxSpectrumTE=zeros(size(omega_list,1),1);
fluxSpectrumTM=zeros(size(omega_list,1),1);
fluxSpectrum=zeros(size(omega_list,1),1);

%% calculating the heat flux for all lossless layers
for i=1:size(omega_list,1)
%parfor i=1:size(omega,1)
    epsilon_list=size(numOfLayer,1);
    for j=1:numOfLayer
       epsilon_list(j,1)=epsilon{j}(i); 
    end   
    omega=omega_list(i);
    %% this function is aimed at calculating the flux spectrum at all lossless layers from source layer 
    %% This layer should be lossless
%    fluxSpectrumTE(i)= quadgk(@(kx) calJoy(omega,thickness_list,epsilon_list,kx,sourceLayer,targetLayer,isLossy),0,Inf);
%    fluxSpectrumTM(i)= quadgk(@(kx) (calJox(omega,thickness_list, epsilon_list,kx,sourceLayer,targetLayer,isLossy)...
%        + calJoz(omega,thickness_list, epsilon_list,kx,sourceLayer,targetLayer,isLossy)),0,Inf);
%    fluxSpectrum(i)=fluxSpectrumTE{i}+fluxSpectrumTM{i};
    %% this function is aimed at calculating the fields if souce is at z=0. Source should be in lossy media
    [exx{i},hyx{i}] = pointJox(omega,thickness_list, epsilon_list,kx,sourceLayer,d,isLossy);
    [eyy{i},hxy{i}] = pointJoy(omega,thickness_list, epsilon_list,kx,sourceLayer,d,isLossy);
    [exz{i},hyz{i}] = pointJoz(omega,thickness_list, epsilon_list,kx,sourceLayer,d,isLossy);
end

%% saving data
save output
%% close matlabpool
%% matlabpool close;
