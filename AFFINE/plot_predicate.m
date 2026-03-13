function pl_pred = plot_predicate(zls,t0,t1,Color)
    pl_pred = plot(polyshape([t0,t1, t1, t0,t0],[zls(1),zls(1), zls(2),zls(2),zls(1)]),'FaceAlpha',0.25,'FaceColor',Color,'EdgeColor','none');
end