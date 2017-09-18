%GOPH 547
%Ye Sun
%Safian Omar Qureshi
%ID: 10086638
%Collaberated with Tian Yu, Younghwan Ryan Ok, Ahtisham Mohummad 

%creating function that will be used later in script. inputs from course notes formula,
%output gravity effect

function [gz]=grav_eff_point(x,xm,m)
 G=6.674*10^(-11);
 r=norm(x-xm);%r is the scalar distance between 2 points
 r3=r^3;
 dz=(x(3)-xm(3)); 
 gz=(G*m*dz)/r3;
end

