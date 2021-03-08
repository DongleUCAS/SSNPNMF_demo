function [E_rmse,E_sad,ang_theta]=evaluation(U,V,A,abf,p,M,N,D)

Aest = U;
sest = V;
% permute results
CRD = corrcoef([A Aest]);
DD = abs(CRD(p+1:2*p,1:p));
perm_mtx = zeros(p,p);
aux=zeros(p,1);
for i=1:p
    [ld cd]=find(max(DD(:))==DD);
    ld=ld(1);cd=cd(1); % in the case of more than one maximum
    perm_mtx(ld,cd)=1;
    DD(:,cd)=aux; DD(ld,:)=aux';
end

perm_mtx=abs(perm_mtx);
Aest = Aest*perm_mtx;
sest = sest'*perm_mtx;
Sest = reshape(sest,[M,N,p]);
sest = sest';

coef_A = sqrt(diag(Aest'*Aest)./diag(A'*A));
Aest  = Aest ./ repmat(coef_A',[D 1]);
coef_s = sqrt(diag(sest'*sest)./diag(abf'*abf));
sest  = sest ./ repmat(coef_s',[p 1]);



E_rmse = sqrt(sum(sum((abf-sest).*(abf-sest)))/(M*N*p));

nA = diag(A'*A);
nAest = diag(Aest'*Aest);
ang_theta = acos( diag(A'*Aest)./sqrt(nA.*nAest) );
E_sad = mean(ang_theta.^2)^.5;


% show the estimations

figure(1),
for i=1:p
    subplot(2,2,i),
    plot(A(:,i),'r'); axis([0 250 0 1])
    hold on
    plot(Aest(:,i),'g'); axis([0 250 0 1])
end
figure(2)
for i=1:p
    subplot(2,p,i),
    imagesc(reshape(abf(i,:),M,N));
    if i==1 title('True abundance'); end
    hold on
    subplot(2,p,p+i),
    imagesc(Sest(:,:,i));
    if i==1 title('Estimated abundance'); end
end

end
