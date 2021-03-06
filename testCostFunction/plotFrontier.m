function pFrontier = plotFrontier(ah, frontier, x)

pFrontier = plot3(ah, x(1,[frontier.inds frontier.inds(1)]), ...
    x(2,[frontier.inds frontier.inds(1)]), ...
    x(3,[frontier.inds frontier.inds(1)]), 'r-',...
    'linewidth',4,'color',rand(1,3));
