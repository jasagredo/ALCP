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


facpol := proc(f::polynom, p::integer, k::integer) #p=2??
local h, v, i, U, g, facg, j, e;
    #1
    h := x;
    v := f/lcoeff(f);
    i := 0;
    U := [];

    do
        #2
        i := i + 1;
        h := Powmod(h, p^k, f, x) mod p;
        g := Gcd(h-x, v) mod p;

        #3
        if g <> 1 then
            facg := edf(g, p, k, i);

            #4
            for j to nops(facg) do
                e := 0;
                while Rem(v, facg[j], x) mod p = 0 do
                    v := Quo(v,facg[j], x) mod p;
                    e := e + 1
                end do;
                U := [op(U), [facg[j],e]]
            end do;
        end if;
        if v = 1 then
            break
        end if;
    end do;
    return U
end proc:

pa := 5;
a := x^4+3*x^3+4*x^2+3;
facpol(a, pa, 1);
(`mod`(Factors(a), pa))[2];

pb := 3;
b := x^11+2*x^9+2*x^8+x^6+x^5+2*x^3+2*x^2+1;
facpol(b, pb, 1);
(`mod`(Factors(b), pb))[2];

pc := 13;
c := Randpoly(42, x) mod pc:
c := Quo(c, lcoeff(c), x) mod pc;
facpol(c, pc, 1);
(Factors(c) mod pc)[2];
