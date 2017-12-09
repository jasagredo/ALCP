testaks := proc(n::integer)
local a, r, j;
    a := iperfpow(n, 'b');
    if not a then
        if a > 1 and b > 1 then
            return false
        end if;
    end if;

    r := 2;
    do
        if igcd(n, r) > 1 or (igcd(n,r) = 1 and NumberTheory[MultiplicativeOrder](n, r) > 4*(degree(n)+1)^2) then
            break
        end if;
        r := r + 1
    end do;

    if r=n then
        return true
    end if;

    if igcd(n,r) > 1 then
        return false
    end if;

    for j to 2*(degree(n)+1)*floor(sqrt(r))+1 do
        if Rem((x+j)^n - (x^n + j), x^r-1, x) mod n <> 0 then
            return false
        end if;
    end do;

    return true
end proc:

testaks(41);

testaks(1223);

testaks(1219);
