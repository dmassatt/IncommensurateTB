function value = intra_graphene(x,y,o,m)
delta = .01;
rs = x.^2+y.^2;

% d_1 = 1.43, d_2 = 2.4768, d_3 = 3.7834, d_4 = 4.29
% d_1^2 = 2.0449, d_2^2 = 6.1345, d_3^2 = 14.3141, d_4^2 = 18.4041
% t_1 = -2.892, t_2 = .243, t_3 = -.266,t_4 = .024


value = (-2.892)*(rs > 2.0449-delta & rs < 2.0449+delta)...
    + .243*(rs > 6.1345-delta & rs < 6.1345+delta)...
    - .266*(rs > 14.3141 - delta & rs < 14.3141 + delta)...
    +.024*(rs > 18.4041 - delta & rs < 18.4041 + delta);




% value = 0;
% if rs > 2.0449-delta && rs < 2.0449+delta
%     value = -2.892;
% elseif rs > 6.1345-delta && rs < 6.1345+delta
%     value = .243;
% elseif rs > 14.3141 - delta && rs < 14.3141 + delta
%     value = -.266;
% elseif rs > 18.4041 - delta && rs < 18.4041 + delta
%     value = .024;
% end

end