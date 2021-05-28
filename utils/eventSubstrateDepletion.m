function [position,isterminal,direction] = eventSubstrateDepletion(t,y,p,m)

  % Look for the index of the substrate.
  ind = find(contains(m.varsName,'bio__s'));

  position = y(ind) - m.varsStart(ind)*0.02; % The value that we want to be zero

  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction

end
