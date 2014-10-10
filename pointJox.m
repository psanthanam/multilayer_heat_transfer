function [ e,h ] = pointJox(omega,thickness_list, epsilon_list,kx,sourceLayer,d,targetLayer)
%This function calculates the flux when source joy is placed at a given point
%   omega:          frequency
%   thickness_list: thickness of every layer, should start and end with Inf
%   epsilon_list:   epsilon for every layer, complex vector
%   kx:             the parellel wave vector
%   sourceLayer:    where the source is. Origin point is put at the upper
%                   surface of layer 1
%   d:              if the source layer is layer 1, then this is the z
%                   coordinate of the source, otherwise it is the relative 
%                   position to the lower bound of source layer 
%   targetLayer:    the layer where the probe is
%   e,h:            calculated electric and magnetic fields at the target
%                   layer, if the target layer has infinite thickness, the 
%                   fields are defined at the surface, otherwise they are 
%                   defined at the upper bound of target layer

PhysicsConst;

%% initializing layers
numOfLayer=size(thickness_list,1);
if sourceLayer==1
    thickness_list=[Inf;-d;thickness_list(2:end)];
    epsilon_list=[epsilon_list(1);epsilon_list];
elseif sourceLayer == numOfLayer
    thickness_list=[thickness_list(1:sourceLayer-1);d;Inf];
    epsilon_list=[epsilon_list;epsilon_list(numOfLayer)];
else
    thickness_list=[thickness_list(1:sourceLayer-1);d;thickness_list(sourceLayer)-d;thickness_list(sourceLayer+1:end)];
    epsilon_list=[epsilon_list(1:sourceLayer-1);epsilon_list(sourceLayer);epsilon_list(sourceLayer);epsilon_list(sourceLayer+1:end)];
end
numOfLayer=numOfLayer+1;
if targetLayer>sourceLayer
    targetLayer=targetLayer+1;
end
f_list=zeros(numOfLayer,1);
I=cell(numOfLayer-1,1);
M=cell(numOfLayer,1);
kz=zeros(numOfLayer,1);
S_Up=eye(2,2);
S_Down=eye(2,2);
S_Target=eye(2,2);

% epsilon_list,thickness_list,sourceLayer,targetLayer,numOfLayer

    %% defining f_list
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
        M{i}=[kz(i)/(omega*epsilon_list(i)*ep0),-kz(i)/(omega*epsilon_list(i)*ep0);1,1];
    end
    for i=1:numOfLayer-1
        I{i}=M{i}\M{i+1};
    end
    %% updating S_Down
    if sourceLayer~=1
        for i=1:sourceLayer-1
            S_Down(1,1)=1/(I{i}(1,1)-f_list(i)*S_Down(1,2)*I{i}(2,1))*f_list(i)*S_Down(1,1);
            S_Down(1,2)=1/(I{i}(1,1)-f_list(i)*S_Down(1,2)*I{i}(2,1))*...
                (I{i}(2,2)*f_list(i)*S_Down(1,2)-I{i}(1,2))*f_list(i+1);
            S_Down(2,1)=S_Down(2,1)+I{i}(2,1)*S_Down(2,2)*S_Down(1,1);
            S_Down(2,2)=S_Down(2,2)*(I{i}(2,1)*S_Down(1,2)+I{i}(2,2)*f_list(i+1));
        end
    end
    %% updating S_up
    if sourceLayer~=numOfLayer-1
        for i=sourceLayer+1:numOfLayer-1
            S_Up(1,1)=1/(I{i}(1,1)-f_list(i)*S_Up(1,2)*I{i}(2,1))*f_list(i)*S_Up(1,1);
            S_Up(1,2)=1/(I{i}(1,1)-f_list(i)*S_Up(1,2)*I{i}(2,1))*...
                (I{i}(2,2)*f_list(i)*S_Up(1,2)-I{i}(1,2))*f_list(i+1);
            S_Up(2,1)=S_Up(2,1)+I{i}(2,1)*S_Up(2,2)*S_Up(1,1);
            S_Up(2,2)=S_Up(2,2)*(I{i}(2,1)*S_Up(1,2)+I{i}(2,2)*f_list(i+1));
        end
    end
    %% solving for the a_l and b_l
    temp=zeros(2,2);
    temp(1,1)=M{sourceLayer+1}(1,1)+M{sourceLayer+1}(1,2)*S_Up(2,1)*f_list(sourceLayer+1);
    temp(1,2)=-(M{sourceLayer}(1,1)*S_Down(1,2)*f_list(sourceLayer)+M{sourceLayer}(1,2));
    temp(2,1)=M{sourceLayer+1}(2,1)+M{sourceLayer+1}(2,2)*S_Up(2,1)*f_list(sourceLayer+1);
    temp(2,2)=-(M{sourceLayer}(2,1)*S_Down(1,2)*f_list(sourceLayer)+M{sourceLayer}(2,2));
    sourceLayerFields=temp\[0;-1];
    bSource=sourceLayerFields(2);
    aSourceUp=sourceLayerFields(1);
    bSourceUp=S_Up(2,1)*aSourceUp;
    b1=S_Down(2,2)*bSource;
    a1=0;
    if targetLayer==1
        targetFields=M{1}*[a1;b1];
    elseif targetLayer == sourceLayer+1
        targetFields=M{targetLayer}*[aSourceUp*f_list(targetLayer);bSourceUp];
    %% case when target layer is below sourcelayer, from layer 1
    elseif targetLayer<sourceLayer
        for i=1:targetLayer-1
            S_Target(1,1)=1/(I{i}(1,1)-f_list(i)*S_Target(1,2)*I{i}(2,1))*f_list(i)*S_Target(1,1);
            S_Target(1,2)=1/(I{i}(1,1)-f_list(i)*S_Target(1,2)*I{i}(2,1))*...
                (I{i}(2,2)*f_list(i)*S_Target(1,2)-I{i}(1,2))*f_list(i+1);      
            S_Target(2,1)=S_Target(2,1)+I{i}(2,1)*S_Target(2,2)*S_Target(1,1);      
            S_Target(2,2)=S_Target(2,2)*(I{i}(2,1)*S_Target(1,2)+I{i}(2,2)*f_list(i+1));
        end
        b_Target=b1/S_Target(2,2);
        a_Target=S_Target(1,2)*b_Target;
        
        targetFields=M{targetLayer}*[a_Target*f_list(targetLayer);b_Target];
   %% case when target layer is above sourcelayer, from source layer + 1
    else
        for i=sourceLayer+1:targetLayer-1
            S_Target(1,1)=1/(I{i}(1,1)-f_list(i)*S_Target(1,2)*I{i}(2,1))*f_list(i)*S_Target(1,1);
            S_Target(1,2)=1/(I{i}(1,1)-f_list(i)*S_Target(1,2)*I{i}(2,1))*...
                (I{i}(2,2)*f_list(i)*S_Target(1,2)-I{i}(1,2))*f_list(i+1);      
            S_Target(2,1)=S_Target(2,1)+I{i}(2,1)*S_Target(2,2)*S_Target(1,1);      
            S_Target(2,2)=S_Target(2,2)*(I{i}(2,1)*S_Target(1,2)+I{i}(2,2)*f_list(i+1));
        end
        b_Target=(bSourceUp-S_Target(2,1)*aSourceUp)/S_Target(2,2);
        a_Target=S_Target(1,1)*aSourceUp+S_Target(1,2)*b_Target;
        targetFields=M{targetLayer}*[a_Target*f_list(targetLayer);b_Target];
    end
    %% solving for the E and H fields.
    e=targetFields(1);  
    h=targetFields(2);
    %e,h
end
