function B=Norm01(A)

B=(A-min(A(:)))/(max(A(:))-min(A(:)));