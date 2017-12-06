INV := proc(f::polynom, p::integer, k::integer)
	local inv;
	Gcdex(f,Randprime(k,r) mod p, r, 'inv') mod p;
	return inv;
end proc:

p := 5: k := 3: minpol := x^3+x+1:
alias(r=RootOf(minpol));
f := 4*r^2+ 3;
inv := INV(f,p,k);
Expand(f*inv) mod p;
