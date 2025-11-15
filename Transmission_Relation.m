function beta_j = Transmission_Relation(ii,X)

if(isempty(X))
    X=[-505.266972443922	0.0925585588236106	0.194908155731669	1.14018129218886	17.8796481564311];
end

beta_j=X(4)./(1+exp(X(1).*(ii-X(2)).*(ii-X(3))));
end