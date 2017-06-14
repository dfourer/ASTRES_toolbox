function ort = io(x,imf)
% ort = IO(x,imf) computes the index of orthogonality
%
% inputs : - x    : analyzed signal
%          - imf  : empirical mode decomposition

n = size(imf,1);

s = 0;

for i = 1:n
  for j =1:n
    if i~=j
      s = s + abs(sum(imf(i,:).*conj(imf(j,:)))/sum(x.^2));
    end
  end
end

ort = 0.5*s;
end