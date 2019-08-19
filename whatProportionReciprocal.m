function pRecip = whatProportionReciprocal(adjMatrix);

% There is a connection in either direction:
adjMatrixConnEither = (adjMatrix==1 | adjMatrix'==1);
% There is a connection in both directions:
adjMatrixConnBoth = (adjMatrix==1 & adjMatrix'==1);
adjMatrixConnEitherUpper = adjMatrixConnEither(triu(true(size(adjMatrix))));
adjMatrixConnBothUpper = adjMatrixConnBoth(triu(true(size(adjMatrix))));
pRecip = mean(adjMatrixConnBothUpper(adjMatrixConnEitherUpper==1));
fprintf(1,'%.1f%% reciprocal connectivity rate\n',pRecip*100);

end
