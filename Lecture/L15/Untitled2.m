%1
y=linspace(0,3,6);
x=sqrt(y+1);

Qt=trapz(y,x)
%% 2
f=@(r) (5*r)./(4+r.^2).^2;

Qi=integral(f,-1,1,'AbsTol',1e.-8)

%% 3
f3=@(theta)(cos(2* theta).*sin(2* theta)).^-3
upperLimit= pi/6
Qi=integral(f3,0,upperLimit, 'AbsTol', 1.e-8)

%%
x=linspace(0,5,101);
y=x.*exp(-x);

figure
plot(x,y,'g-')
title('Function Plot')
xlabel('x')
ylabel('y')

I=@(a) 1.-exp(-a)-a*exp(-a);
i5=I(5)

Qt=trapz(x,y)

f=@(q) q.*exp(-q);
Qi=integral(f,0,5,'AbsTol',1.e-4)

%Exact value is 0.9596. The trapz method yielded 0.9594. The integral
%method produced a value of 0.9596.

