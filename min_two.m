function [ m, I ] = min_two( A )
%MIN_TWO Summary of this function goes here
%   Detailed explanation goes here
[rows,cols] = size(A);

m1 = [];
m2 = [];
x1 = [];
x2 = [];
y1 = [];
y2 = [];

[m1, I_t] = min([ min(A(1,:)) min(A(2,:)) ]);

if I_t == 1
    [m2, y2] = min(A(2,:));
    [~,  y1] = min(A(1,:));
    x1 = 1; x2 = 2;
else
    [m2, y2] = min(A(1,:));
    [~,  y1] = min(A(2,:));
    x1 = 2; x2 = 1;
end

for i = 1:rows
    [min_i, I_i] = min( A(i,:) );
    
    if min_i < m1
        m2 = m1;
        x2 = x1; y2 = y1;
        m1 = min_i;
        x1 = i; y1 = I_i;
        
    elseif min_i >= m1 & min_i < m2
            m2 = min_i;
            x2 = i; y2 = I_i;
    end
end

m = [m1; m2];
I = [x1 y1; x2 y2];

end

