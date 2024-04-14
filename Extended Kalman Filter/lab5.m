clc;
clear all;

data = readtable("data.txt");
ground_t = data(:,1);
measurements = data(:,2);

t = 1;
q = 0.001;
r = 0.1;

dfda = [0 0 0; 0 1 0; 0 0 0];
dgdn = 1;
dgdx = [0 0 1];

Q = [0 0 0; 0 q 0; 0 0 0];
R = r;

I = eye(3);
S = eye(3);
X = [0;0;measurements(1,"Var2")];