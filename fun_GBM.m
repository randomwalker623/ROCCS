function SS = fun_GBM( dt,mu1,sigma1,N,M,S0 )
ST=ones(N,M);
for i=1:N
  for j=1:M
    ST(:,j)=randn(N,1);
  end
end
STT= ST';

for i=1:M
  SS(i,1)= S0;
    for j=2:N
      SS(i,j)= SS(i,j-1)*exp(mu1*dt+(sigma1)*sqrt(dt)*STT(i,j));
    end
end

end

