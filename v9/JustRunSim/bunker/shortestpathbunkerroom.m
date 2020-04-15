% testing bunker, comparing my shortest mean path function with matlabs
main;

o1 = zeros(100,1);

for i = 1:100
    o1(i) = mean(graphshortestpath(sparse(O.A),i));
end

O1 = mean(o1)


O2 = graphshortD(O.A)
