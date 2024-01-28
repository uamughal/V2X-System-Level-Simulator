% N=StationNumElements(S)
% return a vector with the number of array elements for each array/station
% in S
%
% Author: Martin K�ske (TUI)
function N=StationNumElements(S)
N=zeros(1,length(S));
for i=1:length(S)
    N(i)=length(S(i).Elements);
end;