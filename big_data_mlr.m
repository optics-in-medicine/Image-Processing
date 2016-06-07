
X = []; Y = [];
mf = [5 5];
aby_t = []; aby_b = [];
ppix_t = []; ppix_b = [];
ird_t = []; ird_b = [];


fnames = {'405B_2','405C_2','407B_1','410B_2', '408B_1','406A_2'};

for i = 1:length(fnames)
% pick control points, number them after each click.
load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
    fnames{i} 'wrp.mat'])

    P1 = medfilt2(nm(im1a.I_a), mf);
    P2 = medfilt2(nm(im2a.I_a), mf);
    P3 = medfilt2(nm(im3a.I_a), mf);


    X1 = [P1(mask.tumor) P2(mask.tumor) P3(mask.tumor)];
    Y1 = ones(sum(mask.tumor(:)),1);

    X_add = [X1; [P1(mask.brain) P2(mask.brain)  P3(mask.brain)]];
    Y_add = [Y1; zeros(sum(mask.brain(:)),1)];

    X = [X; X_add];
    Y = [Y; Y_add];
    
    aby_t = [aby_t; P1(mask.tumor)];
    aby_b = [aby_b; P1(mask.brain)];
    
    ppix_t = [ppix_t; P2(mask.tumor)];
    ppix_b = [ppix_b; P2(mask.brain)];
    
    ird_t = [ird_t; P3(mask.tumor)];
    ird_b = [ird_b; P3(mask.brain)];
    

end

[mp, prob1, r1, t01] = get_logistic( aby_t, aby_b, 100 );
[mp, prob2, r2, t02] = get_logistic( ppix_t, ppix_b, 100 );
[mp, prob3, r3, t03] = get_logistic( ird_t, ird_b, 100 );

figure;
plot(mp,prob1,'.',mp,logistic_func(mp,1,r1,t01));


prob_t = []; prob_b = [];
fnames = {'405B_2','405C_2','407A_2','407B_1','410B_2', '408B_1','406A_2'};
for i = 1:length(fnames)
    
    load(['X:\#5 - Data\# 2015 RAT F98 ZEISS IMAGING\# F98 Rats Batch 4 (F401-F415)\processed\',...
    fnames{i} 'wrp.mat'])

    mprob = (logistic_func(nm(im1a.I_a),1,r1,t01) + ...
             logistic_func(nm(im2a.I_a),1,r2,t02)) ./ 2;

    prob_t = [prob_t; mprob(mask.tumor)];
    prob_b = [prob_b; mprob(mask.brain)];
     
    
    subplot(2,3,i)
    imagesc(mprob); axis image; axis off; colormap(goodmap('comet'));
end


[X,Y,T,AUC] = perfcurve([ones(length(prob_t),1);zeros(length(prob_b),1)]+1,[aby_t; aby_b],2);
plot(X,Y);
hold on;
[X,Y,T,AUC] = perfcurve([ones(length(prob_t),1);zeros(length(prob_b),1)]+1,[ppix_t; ppix_b],2);
plot(X,Y);
[X,Y,T,AUC] = perfcurve([ones(length(prob_t),1);zeros(length(prob_b),1)]+1,[ird_t; ird_b],2);
plot(X,Y);
[X,Y,T,AUC] = perfcurve([ones(length(prob_t),1);zeros(length(prob_b),1)]+1,[prob_t; prob_b],2);
plot(X,Y)
legend('ABY029','PpIX','IRD680','MTCI')
ylabel('Sensitivity')
xlabel('1-Specificity');
