function out = uniRand(a,b)
%%% 22-11-20        first revision
%%% Benoit Duchet, University of Oxford

%%% out is a uniformly distributed random number between a and b

range = b - a;
out = rand(1,1)*range + a;

end
