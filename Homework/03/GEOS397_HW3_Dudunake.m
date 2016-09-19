%% Homework 3: Presented by Taylor Dudunake and Javier Colton
%
%
%% Part 1: Finite Difference Approximation
% 
%% _Step 1: Taylor Series Expansion_
% This is the Taylor Series expansion of the fininite difference 
% approximation with respect to time _t_ using the variable _k_:
%
% $$ f(x,t+k) = f(x,t) + k\frac{\partial f(x,t)}{\partial t} 
% + \frac{k^2}{2} \frac{\partial^2 f(x,t)}{\partial t^2} + \frac{k^3}{6} 
% \frac{\partial^3 f(x,t)}{\partial t^3} + .... $$
%
%% _Step 2: The forward difference operator for first derivatives_
% The term _O(h)_ represents the difference between the result of the
% equation and what the real answer is expected to be. Also, it's the
% difference between the terms of order h and higher
%
% The first order _time_ derivative using the forward difference operator
% is:
%
% $$ \frac{\partial f(x,t)}{\partial t} = \frac{f(x, t+k) - f(x)}{k} + 
% \frac{k}{2} \frac{\partial^2 f(x,t)}{\partial t^2} = \frac{f(x, t+k)
% - f(x,t)}{k} + O(k) $$
%
%% _Step 3: The centered difference operator for second derivatives_
% This is centered difference approximation coming from the second
% derivative in the lecture 05 notes: 
%
% $$ \frac{\partial^2 f(x)}{\partial x^2} \approx \frac{f(x+h) -
% 2f(x)+f(x-h)}{h^2} $$
%
%% _Step 4: Approximate the derivative_
% Next, we calculate the approximation of the first order time derivative
% which is as follows:
%
% $$ \frac{f(x, t+k)-f(x,t)}{k} = \frac{\partial f(x,t)}{\partial t} $$
%
%% Part 2: The finite difference solution to the diffusion equation
%
%% _Step 1: Approximate the partial differential equation_
% Rewritten in terms of the finite difference approximations:
%
% $$ \frac{\partial z}{\partial t} = \kappa \frac{\partial^2 z}{\partial
% x^2} $$
% 
% $$ \frac{z(x, t+k) - z(x,t)}{k} = \kappa \frac{z(x+h, t)-2z(x, t) +
% z(x-h,t)}{h^2} $$
%
%% _Step 2: Solve for the values of the function at time t+k_
% In code the following code, we use k = dt and h = dx. In the following
% equations, we have solved the equation for z(x, t+k):
%
% $$ z(x, t+k) = k(\kappa \frac{z(x+h,t) - 2z(x,t) + z(x-h,t)}{h^2} +
% z(x,t) $$
%
% $$ z(x, t+dt) = \kappa dt \frac{z(x+dx,t) - 2z(x,t) + z(x-dx,t)}{\partial
% x^2} + z(x,t)$$
%
%% Part 3: Implementing the numerical solution
%
%% _Step 1: Define parameters and constants_
%
clc 
clear all
%
%
dt = 1; %years (yr)
dx = 1; %meters (m)
k = 2 * exp(-3); %square meters per year (m^2/yr)
%
%% _Step 2: Make the initial model_
%
% The following code generates our triangle moraine at t = 1 year
z = [0 0 0 0 0 1 2 3 4 5 6 7 8 9 10 9 8 7 6 5 4 3 2 1 0 0 0 0 0];
nNode = numel(z);
xArray = (0 : nNode - 1) * dx;
%
% The fifteenth element of the array (xArray=14) is the maximum elevation
% in this model. This is the same thing as z(15) = 10.
%
z(15);
%
%%
% *Figure and figure properties:*
H = plot(xArray, z);
hold on
title ('Topography of Moraine');
xlabel ('Meters');
ylabel ('Height of Moraine (m)');
axis([0, 28, 0, 12]);
legend toggle;
legend ('Initial Topography');
set(gcf,'PaperUnits', 'inches','PaperPosition', [0 0 4 4]);
saveas (H, 'initialtopography', 'png');

%% _Step 3: Loop through time to compute the topography at t+dt_
tMax = 100;
t0 = dt;
%
for t = (t0+dt:dt:tMax);
    
    for x = (2:nNode-1);
        z(1) = 0;
        z(nNode) = 0;
        
        z(x) = ((z(x-dx)+z(x+dx)-2*z(x))*dt*k)/(dx^2)+z(x);
    end
end

%% _Step 4: Plot results_
%
G = plot(xArray, z);
legend('Initial Elevation', 'Final Elevation');
hold off
%% Part 4: Discussion 
% When we increase the time of which the moraine experiences erosion to
% tMax = 1e2, 1e3, 1e4, 1e5, and 1e6 years, we see varying levels of 
% erosion and transport. At 100 years,  the height of the moraine is 6.3440 meters. At 
% 1000 years, the height of the moraine is 1.603 meters. At 10,000 years,
% the height of the moraine is 5.88e-06] meters. At 100,000 years, the height
% of the moraine is 2.57e-60 meters. Finally, at 1,000,000 years, the
% height of the moraine is 2.39e-321 meters. We notice that at some point
% between 1000 years and 10,000 years, the moraine is completely eroded
% away and the sediments are transported outwards from the base of the
% moraine. Logically, this makes sense from my knowlde of soil erosion and
% transport.
%
% When we vary the diffusivity constant but hold the maximum time constant
% (100 years), we see different results. When the value is significantly 
% smaller than 2 * exp(-3), such as 2 * exp(-15),the moraine experience 
% hardly any erosion with a maximum height of 9.99 meters. However, if we 
% increase the diffusivity constant to a value such as 2 * exp(-1), we see 
% a lot of erosion resulting in a maximum height of 0.19 meters.