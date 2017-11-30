with(GaussInt, GIrem, GIquo, GInormal):

EUCLIDESEXT := overload([
	proc (a::integer, b::integer, p::integer)
		local r0, r1, r2, s0, s1, s2, t0, t1, t2, q;
		option overload;
		if nargs = 2 then
			r0, r1, s0, s1, t0, t1 := a, b, 1, 0, 0, 1;
			while r1 <> 0 do
				q := iquo(r0, r1);
				r2 := irem(r0, r1);
				s2 := s0 - s1*q;
				t2 := t0 - t1*q;
				r0, s0, t0 := r1, s1, t1;
				r1, s1, t1 := r2, s2, t2
			end do;
			return r0, s0, t0
		end if;
		r0, r1, s0, s1, t0, t1 := a mod p, b mod p, 1, 0, 0, 1;
		while r1 <> 0 do
			q := iquo(r0, r1) mod p;
			r2 := irem(r0, r1) mod p;
			s2 := s0 - s1*q;
			t2 := t0 - t1*q;
			r0, s0, t0 := r1, s1, t1;
			r1, s1, t1 := r2, s2, t2
		end do;
		return r0, s0, t0
	end proc,

	proc (a::(complex(integer)), b::(complex(integer)))
		local r0, r1, r2, s0, s1, s2, t0, t1, t2, q;
		option overload;
		r0, r1, s0, s1, t0, t1 := a, b, Complex(1,0), Complex(0,0), Complex(0,0), Complex(1,0);
		while r1 <> Complex(0,0) do
			q := GIquo(r0, r1);
			r2 := GIrem(r0, r1);
			s2 := s0 - s1*q;
			t2 := t0 - t1*q;
			r0, s0, t0 := r1, s1, t1;
			r1, s1, t1 := r2, s2, t2
		end do;
		return r0, s0, t0;
	end proc,

	proc(a::polynom, b::polynom, p::integer)
		local r0, r1, r2, s0, s1, s2, t0, t1, t2, q, FX;
		option overload;

		if nargs = 3 then
			r0 := a mod p;
			r1 := b mod p
		elif nargs = 2 then
			r0, r1 := a, b
		end if;

		s0, s1, t0, t1 := 1, 0, 0, 1;
		while degree(r1) <> 0 and degree(r1) <> -infinity do
			if nargs = 3 then
				q := Quo(r0, r1, x) mod p;
				r2 := Rem(r0, r1, x) mod p
			else
				q := quo(r0, r1, x);
				r2 := rem(r0, r1, x)
			end if;
			s2 := s0 - s1*q;
			t2 := t0 - t1*q;
			r0, s0, t0 := r1, s1, t1;
			r1, s1, t1 := r2, s2, t2
		end do;
		if nargs = 2 then
      return r0, s0, t0;
    else
      return Expand(r0) mod p, Expand(s0) mod p, Expand(t0) mod p
    end if;
	end proc
	]):


#Ejemplo con Z
print(`Ejemplo con Z`);
a := 2;
b := 140;
mcd, u, v := EUCLIDESEXT(a,b);
evalb(mcd -u*a -v*b=0);

#Ejemplo con Zp
print(`Ejemplo con Zp`);
a := 8;
b := 140;
p := 3;
mcd, u, v := EUCLIDESEXT(a,b,p);

#Ejemplo con GaussInt
print(`Ejemplo con GaussInt`);
p1 := Complex(2,-3):
p2 := Complex(5,2):
p3 := Complex(5,1):
a := p1*p2;
b := p1*p3;
mcd, u, v := EUCLIDESEXT(a,b);

#Ejemplo con Q[x]
print(`Ejemplo con Q[x]`);
a := (x^2 + 2*x - 15)/7;
b := (x^2 + 13*x + 40)/5;
mcd, u, v := EUCLIDESEXT(a,b);

#Euclides con Fp
print(`Ejemplo con Fp`);
p:=7:
a:=6;
b:=3*r^6+2*r^2+r+5;
mcd, u, v :=EUCLIDESEXT(a,b, p);

#Euclides con Fp[x]
print(`Ejemplo con Fp[x]`);
a := 3*x^6+2*x^2+x+5;
b := 6*x^4+x^3+2*x+4;
mcd, u, v:=EUCLIDESEXT(a,b,p);

#Euclides con Fq
print(`Ejemplo con Fq`);
p:=2: k:=4: minpol:=x^4+x+1:
alias(r=RootOf(minpol)):
a:=r^3+r^2;
b:=r*r mod p;
mcd, u, v := EUCLIDESEXT(a, b, p);


print(`Ejemplo con Fq[x]`);
a := r^2*x^2+2*x*r+1;
b := r*x+1;
EUCLIDESEXT(a,b, p);
