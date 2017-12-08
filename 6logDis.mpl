babystep := proc(a, b, p::posint)
	local n, m, alphaA, i, inv, beta, idx2, idx, claseA, claseB;
	if nargs = 4 then
		f := args[4];
		n := p^degree(f);
	else
		n := p;
		f := Randprime(1,x) mod p;
	end if;
	m := ceil(sqrt(n-1));
	claseA := Rem(a, f, x) mod p;
	claseB := Rem(b, f, x) mod p;
	alphaA := table();
	alphaA[0] := 1;
	for i from 1 to m-1 do
		alphaA[i] := alphaA[i-1]*(claseA^(m) mod p) mod p;
	end do;
	Gcdex(claseA, f, x, 'inv') mod p;
	beta := claseB;
	idx2 := 0;
	do
		if member(beta, alphaA, idx) then
			break;
		end if;
		beta := Rem(beta*inv, f, x) mod p;
		idx2 := idx2+1;
	end do;

	return idx*m+idx2;
end proc:

7^6 mod 37;

babystep(7, %, 37);
