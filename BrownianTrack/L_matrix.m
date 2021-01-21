function L = L_matrix(Fk,Dk,ck)
N = 2*length(Fk)-2;
offset = length(Fk);

L = zeros(N+1,N+1);
for r = -(N/2):(N/2)
    for k = -(N/2):(N/2)
            if (r-k + offset > 0) && (r-k +offset <= N+1)
                if r-k < 0
                    L(r + offset,k+offset) = conj(Fk(abs(r-k)+1))*(2*pi*1j*k) ...
                        - conj(Dk(abs(r-k)+1))*(2*pi*k)^2-conj(ck(abs(r-k) +1));
                else
                    L(r + offset,k+offset) = Fk(r-k+1)*(2*pi*1j*k) - Dk(r-k+1)*(2*pi*k)^2-ck(r-k +1);
                end
            end
    end
end

end