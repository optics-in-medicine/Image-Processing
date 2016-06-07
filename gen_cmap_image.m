function [ imc ] = gen_cmap_image( qF, cmap )
%GEN_CMAP_IMAGE Summary of this function goes here
%   Detailed explanation goes here

max_qF = max( qF(:) );                                                         % find max
min_qF = min( qF(:) );
[len_cm, ch] = size( cmap );

qF_n = uint16( floor( (len_cm) .* ( qF - min_qF) ./ (max_qF - min_qF)) + 1);
qF_n(qF_n > len_cm) = len_cm;

% define how each colorchannel maps to values of I_par_n (1-256).
imc = [];
for i = 1:ch
    col_map = cmap(:,i);
    imc = cat(3, imc, col_map( qF_n ));
end

end

