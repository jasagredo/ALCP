irrtest := proc(f::polynom, q::integer)
	local p, n, rep, a, factores, b, i;

	p := NumberTheory[PrimeFactors](q)[1];
	n := degree(f);
	factores := NumberTheory[PrimeFactors](n);
	a := Rem(x^(q^n)-x, f, x) mod p;
	if a <> 0 then
		return `reducible`
	end if;

	for i to numelems(factores) do
		#b := Powmod(x, q^(n/factores[i]),f , x) mod p;
		if Gcd(x^(q^(n/factores[i]))-x,f) mod p <> 1 then
			return `reducible`
		end if;
	end do;
	return `irreducible`
end proc:

p := 2;
k := 4;
minpol := Randprime(k, x) mod p;
alias(r=RootOf(minpol));
f := Randpoly(4, x, r) mod p;
g := Factor(f) mod p;
l := r^2+r+x;
Irreduc(f) mod p;
irrtest(f, 16);
