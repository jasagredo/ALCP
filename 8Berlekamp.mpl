with(PolynomialTools):
with(LinearAlgebra):
berlek := proc(f::polynom, p::integer, k::integer)
    local n, q, i, Q, r, j, QI, b, bbase, a, g1, g2;
	n := degree(f);
	q := p^k;
	i := 0;
    Q := Matrix(n,n);
	while i < n do
		r[i+1] := Powmod(x^q, i, f, x) mod p;
		r[i+1] := CoefficientList(r[i+1], x);
		while numelems(r[i+1]) < n do
			r[i+1] := [op(r[i+1]), 0];
		end do;
		i := i + 1;
	end do;
	for i to n do
		for j to n do
			Q[i, j] := r[i][j];
		end do;
	end do;
	QI := Q-IdentityMatrix(n);
    QI := LUDecomposition(QI, output='U');
    b := NullSpace(QI);
    r := numelems(b);

    for i to r do
        bbase[i] := FromCoefficientVector(b[i], x);
    end do;

    a := 0;
    for i to r do
        a := a + bbase[i]*rand(p)();
    end do;
    g1 := Gcd(a,f) mod p;
    if g1 <> 1 and g1 <> f then
        return g1;
    end if;

    if a=0 and iquo(q-1,2)=0 then
        b := 1;
    else
        b := Powmod(a, iquo(q-1,2), f, x) mod p;
    end if;

    g2 := Gcd(b-1, f) mod p;
    if g2 <> 1 and g2 <> f then
        return g2;
    end if;
    return "failure";
end proc:

berlek(x^3+x+2, 5, 1);
