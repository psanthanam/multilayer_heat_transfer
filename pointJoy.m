function [ e,h ] = pointJoy(omega,thickness_list, epsilon_list,kx,sourceLayer,d,isLossy)
%This function calculates the flux when source joy is placed at a given
%point d
%   omega:          frequency
%   thickness_list: thickness of every layer, should start and end with Inf
%   epsilon_list:   epsilon for every layer, complex vector
%   kx:             the parellel wave vector
%   sourceLayer:    where the source is. Origin point is put at the upper
%                   surface of layer 1
%   d:              if the source layer is layer 1, then this is the z
%                   coordinate of the source, otherwise it is the relative 
%                   position to the lower bound of source layer 
%   isLossy:        vector contain information of material. If lossy, then
%                   the value is 1, otherwise is 0
%   e,h:            calculated electric and magnetic fields at the target
%                   layer, if the target layer has infinite thickness, the 
%                   fields are defined at the surface, otherwise they are 
%                   defined at the upper bound of target layer

PhysicsConst;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initializing layers
numOfLayer=size(thickness_list,1);
if sourceLayer==1
    thickness_list=[Inf;-d;thickness_list(2:end)];
    epsilon_list=[epsilon_list(1);epsilon_list];
    isLossy=[1;isLossy];
elseif sourceLayer == numOfLayer
    thickness_list=[thickness_list(1:sourceLayer-1);d;Inf];
    epsilon_list=[epsilon_list;epsilon_list(numOfLayer)];
    isLossy=[isLossy;1];
else
    thickness_list=[thickness_list(1:sourceLayer-1);d;thickness_list(sourceLayer)-d;thickness_list(sourceLayer+1:end)];
    epsilon_list=[epsilon_list(1:sourceLayer-1);epsilon_list(sourceLayer);epsilon_list(sourceLayer);epsilon_list(sourceLayer+1:end)];
    isLossy=[isLossy(1:sourceLayer-1);1;1;isLossy(sourceLayer+1:end)];
end
numOfLayer=numOfLayer+1;
f_list=zeros(numOfLayer,1);
I=cell(numOfLayer-1,1);
M=cell(numOfLayer,1);
kz=zeros(numOfLayer,1);
S=cell(numOfLayer,1);
e=zeros(numOfLayer-1,1);
h=zeros(numOfLayer-1,1);

%epsilon_list,thickness_list,sourceLayer,numOfLayer,isLossy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% defining f_list, M matrix, I matrix, S matrix
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
        M{i}=[1,1;-kz(i)/(omega*mu0),kz(i)/(omega*mu0)];
        S{i}=eye(2,2);
    end
    for i=1:numOfLayer-1
        I{i}=M{i}\M{i+1};
    end
    %% updating S from layer 2 to source layer
    if sourceLayer~=1
        for i=1:sourceLayer-1
            S{i+1}(1,1)=1/(I{i}(1,1)-f_list(i)*S{i}(1,2)*I{i}(2,1))*f_list(i)*S{i}(1,1);
            S{i+1}(1,2)=1/(I{i}(1,1)-f_list(i)*S{i}(1,2)*S{i}(2,1))*...
                (I{i}(2,2)*f_list(i)*S{i}(1,2)-I{i}(1,2))*f_list(i+1);
            S{i+1}(2,1)=S{i}(2,1)+I{i}(2,1)*S{i}(2,2)*S{i+1}(1,1);
            S{i+1}(2,2)=S{i}(2,2)*(I{i}(2,1)*S{i+1}(1,2)+I{i}(2,2)*f_list(i+1));
        end
    end
    %% updating S from layer sourcelayer+1 to top layer
    if sourceLayer~=numOfLayer-1
        for i=sourceLayer+1:numOfLayer-1
            S{i+1}(1,1)=1/(I{i}(1,1)-f_list(i)*S{i}(1,2)*I{i}(2,1))*f_list(i)*S{i}(1,1);
            S{i+1}(1,2)=1/(I{i}(1,1)-f_list(i)*S{i}(1,2)*I{i}(2,1))*...
                (I{i}(2,2)*f_list(i)*S{i}(1,2)-I{i}(1,2))*f_list(i+1);
            S{i+1}(2,1)=S{i}(2,1)+I{i}(2,1)*S{i}(2,2)*S{i+1}(1,1);
            S{i+1}(2,2)=S{i}(2,2)*(I{i}(2,1)*S{i+1}(1,2)+I{i}(2,2)*f_list(i+1));
        end
    end
    %% solving for the a_l and b_l
    temp=zeros(2,2);
    temp(1,1)=M{sourceLayer+1}(1,1)+M{sourceLayer+1}(1,2)*S{numOfLayer}(2,1)*f_list(sourceLayer+1);
    temp(1,2)=-(M{sourceLayer}(1,1)*S{sourceLayer}(1,2)*f_list(sourceLayer)+M{sourceLayer}(1,2));
    temp(2,1)=M{sourceLayer+1}(2,1)+M{sourceLayer+1}(2,2)*S{numOfLayer}(2,1)*f_list(sourceLayer+1);
    temp(2,2)=-(M{sourceLayer}(2,1)*S{sourceLayer}(1,2)*f_list(sourceLayer)+M{sourceLayer}(2,2));
    sourceLayerFields=temp\[0;1];
    bSource=sourceLayerFields(2);
    aSourceUp=sourceLayerFields(1);
    bSourceUp=S{numOfLayer}(2,1)*aSourceUp;
    b1=S{sourceLayer}(2,2)*bSource;
    a1=0;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% if the first layer is not lossy
    if isLossy(1)==0
        targetFields=M{1}*[a1;b1];
        e(1)=targetFields(1);  
        h(1)=targetFields(2);
    end
    
    %% if the layer above the source is not lossy, actually this is impossible in the setup
    if isLossy(sourceLayer+1)==0
        targetFields=M{sourceLayer+1}*[aSourceUp*f_list(sourceLayer+1);bSourceUp];
        e(sourceLayer)=targetFields(1);  
        h(sourceLayer)=targetFields(2);
    end
 
    %% from layer 2 to sourcelayer, calculate the fields for the lossless layer
    for i=2:sourceLayer
        S_Target=S{i};
        if isLossy(i)==0
            b_Target=b1/S_Target(2,2);
            a_Target=S_Target(1,2)*b_Target;
            targetFields=M{i}*[a_Target*f_list(i);b_Target];
            e(i)=targetFields(1);
            h(i)=targetFields(2);
         end
    end
     
    %% from sourcelayer+1 to the top layer, remember the index of e and h is i-1 because we artifically add a layer.
    for i=sourceLayer+2:numOfLayer
        S_Target=S{i};
        if isLossy(i)==0
            b_Target=(bSourceUp-S_Target(2,1)*aSourceUp)/S_Target(2,2);
            a_Target=S_Target(1,1)*aSourceUp+S_Target(1,2)*b_Target;
            targetFields=M{i}*[a_Target*f_list(i);b_Target];
            e(i-1)=targetFields(1);
            h(i-1)=targetFields(2);
        end
    end

 %   e,h
end
