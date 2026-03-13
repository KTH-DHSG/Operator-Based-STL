function pl_pred_bou = plot_predicate_bounds(zls,t0_min,t0_max,t1_min,t1_max,Color)
    dx = 0.1; 
    pl_pred_bou{1} = plot(polyshape([t0_min,t0_max, t0_max, t0_min,t0_min],[dx + zls(1),dx + zls(1), -dx + zls(2),-dx + zls(2),dx + zls(1)]),'FaceAlpha',0.1,'FaceColor',Color,'EdgeColor',Color,'LineStyle','--');
    pl_pred_bou{2} = plot(polyshape([t1_min,t1_max, t1_max, t1_min,t1_min],[dx + zls(1),dx + zls(1), -dx + zls(2),-dx + zls(2),dx + zls(1)]),'FaceAlpha',0.1,'FaceColor',Color,'EdgeColor',Color,'LineStyle','-.');

end