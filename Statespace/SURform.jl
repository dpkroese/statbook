function SURform(X)
   r, c = size(X);
   idi = vec(Int64.(kronecker(1:r,ones(c))))
   idj = 1:r*c
   return sparse(idi,idj,X'[:])
end
