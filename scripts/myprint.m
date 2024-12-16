fig = gcf;
fig.PaperPositionMode = 'auto';
print(fig, filename, '-depsc', '-vector')
% print(fig, filename, '-depsc', '-opengl')