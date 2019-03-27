% test what fraction of neurons end up positive in 1 distribution

iu1 = -0.01;
isig1 = 0.008;
n1 = 100000000;


I = [ iu1+isig1*randn(1,n1) ];


disp(['Intrinsic freq that are positive = ' mat2str(length(find(I>0))*100/length(I)) '%'])

max(I)