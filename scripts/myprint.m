fig = gcf;
fig.PaperPositionMode = 'auto';
print(fig, filename, '-depsc', '-painters')
% print(fig, filename, '-depsc', '-opengl')