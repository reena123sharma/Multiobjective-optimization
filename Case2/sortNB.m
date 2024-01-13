function [ sol ] = sortNB( Y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
   
     

        Y = sortrows(Y,9); 
        sol(1,:) = Y(1,:);
        Y = sortrows(Y,10);
        sol(2,:) = Y(1,:);
        Y = sortrows(Y,11);
        sol(3,:) = Y(1,:);
   
end