function [position,isterminal,direction] = eventFeedDepletion(t,y,p,m)

  % Look for the index of the substrate.
  V_final     = m.getSymbolByName('bio__V_final').value;
  ind_V_feed  = find(contains(m.varsName,'bio__V_feed'));
  ind_V       = find(contains(m.varsName,'bio__V'));

  position = y(ind_V_feed) + y(ind_V(1)) - V_final; % The value that we want to be zero

  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction

end

