function J=logn_fit_zero(x,Data)

J=log(lognpdf(Data,x(1),x(2)));
J(Data==0)=log(logncdf(1,x(1),x(2)));

J=-sum(J);
end