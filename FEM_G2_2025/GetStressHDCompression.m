function [Stress]=GetStressHDCompression(E0,SigmaMax_C,Strain)
%Caculate stress and stiffness from the Hardin Drenvich Stress Strain relationship
Stress= E0*Strain/(1+E0*Strain/SigmaMax_C);
end
