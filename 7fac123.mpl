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

eds := proc(f::polynom, p::integer, k::integer, d::integer)
local a, g1, b, g2;
    a := Randpoly(rand(degree(f))(), x) mod p;
    if degree(a) < 1 then
        return -1
    end if;

    g1 := Gcd(a, f) mod p;
    if g1 <> 1 then
        return g1
    end if;

    b := Powmod(a, ((p^k)^d-1)/2, f, x)  mod p;
    g2 := Gcd(b-1, f) mod p;
    if g2 <> 1 and g2 <> f then
        return g2
    else
        return -1
    end if;
end proc:

edf := proc(f::polynom, p::integer, k::integer, d::integer)
local n, g;
    n := degree(f);
    if n <= d then
        return [f]
    end if;
    do
        g := eds(f, p, k, d);
        if g <> -1 then
            break
        end if;
    end do;

    return [op(edf(g, p, k, d)), op(edf(Quo(f,g, x) mod p, p, k, d))]
end proc:

#a := x^4+3*x^3+4*x^2+3;
#edf(a, 5, 1, 2);
