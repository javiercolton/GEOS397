clear all
clc

n1=9;
n2=5;

if n1 > n2
    fprintf 'n1 is greater'
else
    if n1 == n2
        fprintf 'equal'
    else
        fprintf 'n2 is greater'
    end
end