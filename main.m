clear;
%% open matlabpool in Hera
%% matlabpool open 4;

materialNames={'InAs';'Air'};
thickness_list=[Inf;Inf];   % thickness should start and end with Inf
numOfLayer=size(thickness_list,1);


epsilon=cell(numOfLayer,1);
for i=1:numOfLayer
    [omega,epsilon{i}]=addMaterial(materialNames{i});
end

exx=zeros(size(omega,1),1);
hyx=zeros(size(omega,1),1);
eyy=zeros(size(omega,1),1);
hxy=zeros(size(omega,1),1);
exz=zeros(size(omega,1),1);
hyz=zeros(size(omega,1),1);

for i=1:size(omega,1)
%parfor i=1:size(omega,1)
    epsilon_list=size(numOfLayer,1);
    for j=1:numOfLayer
       epsilon_list(j)=epsilon{j}(i); 
    end   
    %% this function is aimed at calculating the flux spectrum at layer l=2 when source is at layer 1 
    %% This layer should be lossless
%    fluxSpectrum(i) = quadgk(@(kx) (calJox(omega(i),thickness_list, epsilon_list,kx,2,1)...
%       + calJoy(omega(i),thickness_list, epsilon_list,kx,2,1)...
%        + calJoz(omega(i),thickness_list, epsilon_list,kx,2,1)),0,Inf);
    
    %% this function is aimed at calculating the fields if souce is at z=0. Source should be in lossy media
%    [exx(i),hyx(i)] = pointJox(omega(i),thickness_list, epsilon_list,kx,0);
    [eyy(i),hxy(i)] = pointJoy(omega(i),thickness_list, epsilon_list,kx,0);
%    [exz(i),hyz(i)] = quadgk(@(kx) pointJoz(omega(i),thickness_list, epsilon_list,kx,0),0,Inf);
end

%% saving data
save output
%% close matlabpool
%% matlabpool close;
