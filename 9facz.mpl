with(LinearAlgebra):
with(PolynomialTools):

henselstep := proc(m, f, g, h, s, t)
local e, q, r, gs, hs, b, c, d, ss, ts;
    e := expand(f - g*h)  mod m^2;
    q := Quo(s*e, h, x) mod m^2;
    r := Rem(s*e, h, x) mod m^2;
    gs := expand(g + t*e + q*g) mod m^2;
    hs := expand(h + r) mod m^2;

    b := expand(s*gs+t*hs-1) mod m^2;
    c := Quo(s*b, hs, x) mod m^2;
    d := Rem(s*b, hs, x) mod m^2;
    ss := expand(s - d) mod m^2;
    ts := expand(t - t*b - c*gs) mod m^2;
    return [gs,hs,ss,ts]
end proc:


hensel := proc(f::polynom, g, p, l)
local r, k, d, g0, h0, i, s0, t0, henst, llam1, llam2;
    r := numelems(g);
    if r = 1 then
        return [f/lcoeff(f) mod p^l]
    end if;
    k := floor(r/2);
    d := ceil(log[2](l));

    g0 := lcoeff(f);
    h0 := 1;
    for i to k do
        g0 := expand(g0 * g[i]) mod p
    end do;
    for i from k+1 to r do
        h0 := expand(h0*g[i]) mod p
    end do;

    Gcdex(g0, h0, x, 's0', 't0') mod p;

    for i to d do
        henst := henselstep(p^(2^(i-1)), f, g0, h0, s0, t0);
        g0, h0, s0, t0 := henst[1], henst[2], henst[3], henst[4]
    end do;

    llam1 := hensel(g0, g[1..k], p, l);
    llam2 := hensel(h0, g[k+1..r], p, l);
    return [op(llam1), op(llam2)]
end proc:

facz := proc(f::polynom)
local n, A, b, B, C, gamma, p, cotaPrimo, fb, fbd, l, fmon, rf, T, s, G, fs, S, gs, hs, i, j, TS;
    n := degree(f);
    A := VectorNorm(<CoefficientList(f, x)>);
    if n = 1 then
        return {f}
    end if;
    b := lcoeff(f);
    B := sqrt(n+1)*(2^n)*A*b;
    C := (n+1)^(2*n)*A^(2*n-1);
    gamma := ceil(2*log[2](C));

    p := 4;
    cotaPrimo := 2*gamma*ln(gamma);

    do
        p := nextprime(p);
        fb := f mod p;
        fbd := diff(fb, x) mod p;
        if Gcd(fb, fbd) mod p = 1 and irem(p,b) <> 0 then
            break
        end if;
    end do;
    l := ceil(log[p](2*B+1));

    fmon := f/b mod p;
    rf := map(x -> x[1], (Factors(fmon) mod p)[2]);

    rf := hensel(f/b mod p^l, rf, p, l);
    T := {seq(1..numelems(rf))};
    s := 1;
    G := {};
    fs := f;

    while 2*s <= numelems(T) do
        S := combinat[choose](T,s);
        for i to numelems(S) do
            gs := b;
            hs := b;
            for j to numelems(S[i]) do
                gs := Expand(gs * rf[S[i][j]]) mod p^l
            end do;
            TS := T minus S[i];
            for j to numelems(TS) do
                hs := Expand(hs * rf[TS[j]]) mod p^l
            end do;

            if Norm(CoefficientVector(gs, x), 1) * Norm(CoefficientVector(hs,x), 1) <= evalf(B) then

                T := T minus S;
                G := G union {primpart(gs)};
                fs := primpart(hs);
                b := lcoeff(fs);
                s := s - 1;
                break
            end if;
        end do;
        s := s + 1
    end do;
    return G union {fs}
end proc:

f := 6*x^4+5*x^3+15*x^2+5*x+4;
facz(f);
evalb(expand(mul(%)) = f);
