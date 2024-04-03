function [] = labelTiledLayout(f, sz, shift)
% [] = labelTiledLayout(figure f, size sz, [shift = 0])
% lables figure f with ABCD.. of font size sz. Optional shift to starting letter
if(nargin == 2)
    shift = 0;
end
lables = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890';
t = get(f,'children');
numTiles = t.GridSize(1)*t.GridSize(2);
for i = 1:numTiles
    nexttile(i)
    ax = gca;
    if(strcmp(ax.YDir,'reverse'))
        text(min(xlim), min(ylim),lables(i + shift), 'Horiz','left', 'Vert','bottom','FontSize',sz)
    else
        text(min(xlim), max(ylim),lables(i + shift), 'Horiz','left', 'Vert','bottom','FontSize',sz)
    end
end