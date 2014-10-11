function [ omega,epsilon ] = addMaterial( materialName )
%This function reads in the material name and output the material
%dielectric funcion
  
materialFileName=sprintf('%s%s',materialName,'.mat');
materialData=load(materialFileName);
omega=materialData.omega;
epsilon=materialData.epsilon;
end
