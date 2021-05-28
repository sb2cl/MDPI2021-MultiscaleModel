% First we configure the simulation options.
simOptions AbsTol = 1e-3;
simOptions RelTol = 1e-9;
simOptions TimeSpan  = [0 10];

% Then we define the model.
variable x(start = 0);
parameter k(value = 1);

equation der_x == k - x;
