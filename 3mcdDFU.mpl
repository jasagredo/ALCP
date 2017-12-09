MCDDFU := proc(a::polynom, b::polynom, p::integer)
local r, n, i, ar, an;
    if nargs = 2 then
        r := [a,b]
    elif nargs = 2 then
        r := [a mod p, b mod p]
    end if;
    n := map(t->degree(t),r);
    i := 2;
    while r[i] <> 0 do
        if nargs = 2 then
            ar := primpart(rem(lcoeff(r[i], x)^(n[i-1]-n[i]+1)*r[i-1], r[i], x), x)
        else
            ar := Primpart(Rem(lcoeff(r[i], x)^(n[i-1]-n[i]+1)*r[i-1] mod p, r[i], x) mod p, x) mod p
        end if;
        r := [op(r), ar];
        an := degree(r[i+1]);
        n := [op(n), an];
        i:= i+1
    end do;
    return r[i-1]
end proc:

a := x^3 + 2*x^2 + x + 2;
b := x^3 - 3*x^2 + x - 3;
mcd := MCDDFU(a, b);

a := x^8+2*x^7+4*x^6-2*x^5-3*x^4+4*x^3+9*x^2-2*x-5;
b := 5*x^5-2*x^4-13*x^3-11*x^2+7*x+6;
mcd := MCDDFU(a, b);

a := 824*x^5 - 65*x^4 - 814*x^3 -741*x^2 -979*x - 764;
b := 216*x^4 + 663*x^3 +880*x^2 +916*x+617;
mcd := MCDDFU(a,b);
