function coords = wingbox(geometry)
% Takes geometry (struct with lifting surface geometric parameters)
% Outputs x and y values of corners

% Fit curves to upper and lower surfaces
index   = find(~geometry.polar(:,1));
n       = 10;
disc = linspace(0,1,100);
polyU   = polyfit(geometry.polar(1:index,1), geometry.polar(1:index,2), n);
polyL   = polyfit(geometry.polar(index:end,1), geometry.polar(index:end,2), n);
rsparU  = polyval(polyU, geometry.rspar);
rsparL  = polyval(polyL, geometry.rspar);
fsparU  = polyval(polyU, geometry.fspar);
fsparL  = polyval(polyL, geometry.fspar);

% Find coordinates of corners
ftop    = [geometry.fspar,  fsparU];
fbottom = [geometry.fspar,  fsparL];
rbottom = [geometry.rspar,  rsparL];
rtop    = [geometry.rspar,  rsparU];

% Combine into useful array
coords  = [ftop;
            fbottom;
            rbottom;
            rtop;
            ftop];

%% Uncomment to view plots of polyfit to assess validity
plot(disc,polyval(polyU,disc), "-k", "LineWidth", 2)
hold on
grid on
axis equal
plot(disc,polyval(polyL,disc), "-k", "LineWidth", 2)
plot(geometry.polar(:,1),geometry.polar(:,2), "-k", "LineWidth", 2)
plot(coords(:,1), coords(:,2), "-r", "LineWidth", 2);
hold off

xlabel("Normalised chord, $x/c$", "Interpreter","latex", "FontSize", 16)
ylabel("Normalised height, $y/c$", "Interpreter","latex", "FontSize", 16)
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
% axis equal

end