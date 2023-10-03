function F = root2D(x,mu,C,t2,t1,H,Delta,N_eqk_inside)

FF1 = 0;
FF2 = 0;
for j = 1:N_eqk_inside
    FF1 = FF1+mu(j)/(x(1)*mu(j)+x(2)*C(j));
    FF2 = FF2+(C(j)/(x(1)*mu(j)+x(2)*C(j)));
end
F(1) = FF1-(t2-t1)*Delta;
F(2) = FF2-H;

end