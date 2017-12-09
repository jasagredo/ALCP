with(GaussInt, GIgcd, GIrem):

COPRIMOS := overload([
    proc(a::integer, b::integer)
        return evalb(igcd(a,b) = 1)
    end proc,

    proc(a::complex(integer), b::complex(integer))
        return evalb(GIgcd(a,b) = 1)
    end proc,

    proc(a::polynom, b::polynom, p::integer)
        if nargs = 2 then
            return evalb(gcd(a,b) = 1)
        elif nargs = 3 then
            return evalb(Gcd(a,b) mod p = 1)
        end if;
    end proc
                     ]):

RESTOSCHINO := overload([
    proc(n::list(posint), a::list(integer))
    option overload;
    local i, j, nt, res, coc, mcd, u, v, s, t, c;
        for i from 1 to nops(n) do
            for j from 1 to nops(n) do
                if j <> i then
                    if not COPRIMOS(n[i], n[j]) then
                        error("Los valores introducidos no son coprimos entre si")
                    end if
                end if
            end do
        end do;
        if nops(n) <> nops(a) then
            error("Las listas introducidas tienen dimensiones distintas")
        end if;
        nt := mul(n[i], i = 1..nops(n));
        res := 0;
        for i from 1 to nops(n) do
            coc := nt/n[i];
            mcd := igcdex(coc, n[i], 'u', 'v');
            s[i] := u;
            t[i] := v;
            c[i] := irem(a[i]*s[i],n[i]);
            res := res + (c[i]*nt)/n[i]
        end do;
        return res mod nt
    end proc,

    proc(n::list(complex(integer)), a::list(complex(integer)))
    option overload;
    local i, j, nt, res, coc, mcd, u, v, s, t, c;
        for i from 1 to nops(n) do
            for j from 1 to nops(n) do
                if j <> i then
                    if not COPRIMOS(n[i], n[j]) then
                        error("Los valores introducidos no son coprimos entre si")
                    end if
                end if
            end do
        end do;
        if nops(n) <> nops(a) then
            error("Las listas introducidas tienen dimensiones distintas")
        end if;
        nt := mul(n[i], i = 1..nops(n));
        res := 0;
        for i from 1 to nops(n) do
            coc := nt/n[i];
            mcd := GIgcdex(coc, n[i], 'u', 'v');
            s[i] := u;
            t[i] := v;
            c[i] := GIrem(a[i]*s[i],n[i]);
            res := res + (c[i]*nt)/n[i]
        end do;
        return GImod(res,nt)
    end proc,

    proc(n::list(polynom), a::list(polynom), p::integer)
    option overload;
    local i, j, nt, res, coc, mcd, u, v, s, t, c;
        for i from 1 to nops(n) do
            for j from 1 to nops(n) do
                if j <> i then
                    if nargs = 2 and not COPRIMOS(n[i], n[j]) then
                        error("Los valores introducidos no son coprimos entre si")
                    elif not COPRIMOS(n[i], n[j], p) then
                        error("Los valores introducidos no son coprimos entre si")
                    end if;
                end if
            end do
        end do;
        if nops(n) <> nops(a) then
            error("Las listas introducidas tienen dimensiones distintas")
        end if;
        nt := mul(n[i], i=1..nops(n));
        if nargs = 3 then
            nt := Expand(nt) mod p
        end if;
        res := 0;
        for i from 1 to nops(n) do
            if nargs = 2 then
                coc := nt/n[i];
                mcd := gcdex(coc, n[i], x, 'u', 'v')
            else
                coc := Expand(nt/n[i]) mod p;
                mcd := Gcdex(coc, n[i], x, 'u', 'v') mod p
            end if;
            s[i] := u;
            t[i] := v;
            if nargs = 2 then
                c[i] := rem(a[i]*s[i],n[i], x);
                res := res + (c[i]*nt)/n[i]
            else
                c[i] := Rem(a[i]*s[i] mod p, n[i], x) mod p;
                res := res + ((c[i]*nt mod p)/n[i] mod p) mod p
            end if;
        end do;
        if nargs = 2 then
            return rem(res, nt, x)
        end if;
        return Rem(res, nt, x) mod p
    end proc]);

n := [3, 4, 5];
a := [2, 3, 1];
RESTOSCHINO(n, a);

n := [5, 7, 11];
a := [2, 4, 5];
RESTOSCHINO(n, a);

n := [11, 13];
a := [2,7];
RESTOSCHINO(n, a);
