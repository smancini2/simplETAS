function [A]=nm2N(B)


a=size(B);

nlength=a(1); mlength=a(2);

for n=1:nlength,
    for m=1:mlength,
    A((n-1)*mlength + m ,1)=B(n,m);
    end
end




