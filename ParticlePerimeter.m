% Get the single particle perimeter
%
%   Inputs:
%      image_buffer - n x photodiodes/8 raw image buffer without timestamps
%   Outputs:
%      Perimeter 
%
%   *  Created by Wei Wu, July 4th, 2014 
function [pperimeter] = ParticlePerimeter(image_buffer)

[m, n] = size(image_buffer);
pperimeter = 0;
c1=[49*ones(1,n+2);49*ones(m,1),image_buffer,49*ones(m,1);49*ones(1,n+2)];
for i=2:m+1
    for j=2:n+1
        if (48==c1(i,j) && ( 48~=c1(i+1,j) || 48~=c1(i-1,j) || 48~=c1(i,j+1) || 48~=c1(i,j-1) ) )
            pperimeter = pperimeter+1;
        end
    end
end

end
