function [value,isterminal,direction] = midScapeEvents(t,y)
sizy=size(y);
numFields=sizy(1)/20;
value = y(7*numFields+(numFields/2)-20)-10;     % Detecting when epidemic reaches 10% incidence in 20th last field
isterminal = 1;   
direction = 1;   