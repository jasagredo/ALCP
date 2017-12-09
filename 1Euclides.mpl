with(GaussInt, GIrem, GInormal):

EUCLIDES := overload([
    proc (a::integer, b::integer, p::integer) #Z, Zp
    local r0, r1, r2;
    option overload;
        if nargs > 2 then
            r0, r1 := abs(a) mod p, abs(b) mod p;
            while r1 <> 0 do
                r2 := irem(r0, r1) mod p;
                r0 := r1;
                r1 := r2
            end do;
            return r0 mod p
        end if;
        r0, r1 := abs(a), abs(b);
        while r1 <> 0 do
            r2 := irem(r0, r1);
            r0 := r1;
            r1 := r2
        end do;
        return r0
    end proc,

    proc (a::complex(integer), b::complex(integer)) #GaussIntegers
    local r0, r1, r2;
    option overload;

        r0, r1 := a, b;
        while r1 <> 0 do
            r2 := GIrem(r0, r1);
            r0 := r1;
            r1 := r2
        end do;
        return GInormal(r0)
    end proc,

    proc(a::polynom, b::polynom, p::integer) #Q[x], Fp[x], Fq[x]
    local r0, r1, r2;
    option overload;
        if nargs = 3 then
            r0 := a mod p;
            r1 := b mod p
        elif nargs = 2 then
            r0, r1 := a, b
        else
            error(`Número incorrecto de parámetros`)
        end if;
        while degree(r1) > 0 do
            if nargs = 2 then
                r2 := rem(r0, r1, x)
            else
                r2 := Rem(r0, r1, x) mod p
            end if;
            r0:=r1;
            r1 := r2
        end do;
        if nargs = 2 then
            return r0 / lcoeff(r0)
        end if;
        return r0
    end proc
                     ]):


#Ejemplo con Z
print(`Ejemplo con Z`);
a := 2;
b := 140;
mcd := EUCLIDES(a,b);

#Ejemplo con Zp
print(`Ejemplo con Zp`);
p := 3;
a := 8;
b := 140;
mcd := EUCLIDES(a,b,p);

#Ejemplo con GaussInt
print(`Ejemplo con GaussInt`);
p1 := Complex(2,-3):
p2 := Complex(5,2):
p3 := Complex(5,1):
a := p1*p2;
b := p1*p3;
mcd := EUCLIDES(a,b);

#Ejemplo con Q[x]
print(`Ejemplo con Q[x]`);
a := (x^2 + 2*x - 15)/7;
b := (x^2 + 13*x + 40)/5;
mcd := EUCLIDES(a,b);

#Ejemplo con Fp
print(`Ejemplo con Fp`);
p:=7:
a:= 6;
b := 3*r^6+2*r^2+r+5;
mcd:=EUCLIDES(a, b, p);


#Ejemplo con Fp[x]
print(`Ejemplo con Fp[x]`);
a := 6^2*x^2 + 2*6*x + 1;
b := x+6;
mcd := EUCLIDES(a, b, p);

#Ejemplo con Fq
print(`Ejemplo con Fq`);
p:=2: k:=4: minpol:=x^4+x+1;
alias(r=RootOf(minpol)):
a := r^3+r^2;
b := a*a mod p;
mcd := EUCLIDES(a, b, p);

#Ejemplo con Fq[x]
print(`Ejemplo con Fq[x]`);
a := Expand(r^2*x) mod p;
b := Expand(r^3*x) mod p;
mcd := EUCLIDES(b, a, p);
