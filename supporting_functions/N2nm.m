function [A]=N2nm(B,nlength,mlength)


for n=1:nlength
A(:,n)=B((n-1)*mlength +1:n*mlength,1);
end
    