sfd := proc(poly::polynom, p::integer)
	local L, s, j, g, h, m, f;
	f := poly;
	L := [];
	s := 1;
	do
		j := 1;
		g := Quo(f, Gcd(f, diff(f,x)) mod p, x) mod p;
		while g <> 1 do
			f := Quo(f,g,x) mod p;
			h := Gcd(f,g) mod p;
			m := Quo(g,h,x) mod p;
			if m <> 1 then
				L := [op(L),[m, j*s]];
			end if;
			g := h;
			j := j+1;
		end do;
		if f <> 1 then
			f := f^(1/p) mod p;
			s := p*s;
		end if;
		if f = 1 then
			break;
		end if;
	end do;
	return L;
end proc:

sfd(x^3+3*x^2-4, 5);

ddf := proc(poly::polynom, p::integer, ex::integer)
	local L, h, k, g, f, q;
	q := p^ex;
	f := poly;
	L := [];
	h := Rem(x, f, x) mod p;
	k := 0;
	while f <> 1 do
		h := Powmod(h, q, f, x) mod p;
		k := k + 1;
		g := Gcd(h-x, f) mod p;
		if g <> 1 then
			L := [op(L), [g,k]];
			f := Quo(f,g,x) mod p;
			h := Rem(h,f,x) mod p;
		end if;
	end do;
	return L;
end proc:

ddf(x^4+3*x^3+4*x^2+3, 5, 4);

edf := proc(f::polynom, k::posint, p::posint, q::posint)
	local l, Mk, r, H, Hp, randdeg, a, d, i, deg, h;
	l := degree(f);
	Mk := sum(x^(2^j), j=0..q*k-1);
	r := l/k;
	H := {f};
	while numelems(H) < r do
		Hp := {};
		for i to numelems(H) do
			h := H[i];
			randdeg := rand(q);
			deg := randdeg();
			while deg = 0 do
				deg := randdeg();
			end do;
			aleat := Randpoly(deg, x, alpha) mod p;
			a := Rem(aleat, h, x) mod p;
			Mks := Expand(subs(x=alpha, Mk)) mod p;
			d := Gcd(Mks, h) mod p;
			if d = 1 or d = h then
				Hp := Hp union {h};
			else
				Hp := Hp union {d, Quo(h, d) mod p};
			end if;
		end do;
		H := Hp;
	end do;
	return H;
end proc:
#stopat(edf);
#edf(x^2+5*x+4, 1, 7, 2 );
#unstopat(edf);
