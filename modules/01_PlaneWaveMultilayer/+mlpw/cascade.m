function S=cascade(A,B)
%CASCADE Redheffer product: A precedes B, incoming [left;right].
d=1-A(2,2)*B(1,1);
S=[A(1,1)+A(1,2)*B(1,1)*A(2,1)/d, A(1,2)*B(1,2)/d; ...
   B(2,1)*A(2,1)/d, B(2,2)+B(2,1)*A(2,2)*B(1,2)/d];
end
