INV := proc(f::polynom, p::integer, k::integer)
	local inv;
	Gcdex(f,Randprime(k,x) mod p, x, 'inv') mod p;
	return inv;
end proc:

p := 5: k := 3: minpol := r^3+r+1:
alias(x=RootOf(minpol));
f := 4*x^2+ 3;
inv := INV(f,p,k);
Expand(f*inv) mod p;
