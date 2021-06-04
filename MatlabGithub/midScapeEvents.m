function [value,isterminal,direction] = midScapeEvents(t,y)
sizy=size(y);
numFields=sizy(1)/20;
value = y(7*numFields+(numFields/2)-20)-10;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 1;   % Negative direction only