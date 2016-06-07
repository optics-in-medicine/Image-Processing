function y = logistic_func(x,xdata)
K = x(3);
r = x(2);
t0 = x(1);


y=K./(1+exp(-r*(xdata-t0)));