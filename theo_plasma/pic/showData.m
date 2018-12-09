function showData
global t ekin epot etherm
t = []; ekin = []; epot = []; etherm = [];
files = dir('data/partData.*');
for k = 1:numel(files)
    filenameP = [ 'data/', files(k).name ];
    showParticles( filenameP );
    filenameF = strrep(filenameP, 'part', 'field');
    showFields( filenameF );
    pause(.2)
end
end

function showParticles( fileName )
fid = fopen(fileName, 'r');
np = fscanf(fid, '# np = %d', [1 1]);
time = fscanf(fid, ' time = %g *', [1 1]);
pData = fread(fid, [2, np], 'double');
fclose(fid);
xp = pData(1, :);
vp = pData(2, :);

subplot(3, 2, 1);
plot( xp, vp, '.' );
title( [ 'time = ', num2str(time) ] );
xlabel('x');
ylabel('v');
xlim([0, 3]);
ylim([-3, 3]); % osz.
%ylim([.15, .25]); % freiflug
%ylim([.0, 1.5]); % zweistr
end

function showFields( fileName )
global t ekin epot etherm

fid = fopen(fileName, 'r');

nx = fscanf(fid, '# nx = %d', [1 1]);
time = fscanf(fid, ' time = %g *', [1 1]);
f = fread(fid, [nx, 5], 'double'); % fields x, n, u, T, Ex
fclose(fid);
x = f(:, 1);
n = f(:, 2);
u = f(:, 3);
T = f(:, 4);
Ex = f(:, 5);

dx = x(2) - x(1);   

subplot(3, 2, 2);
plot( x, n, '-*' );
title('n');

subplot(3, 2, 3);
plot( x, u, '-*' );
title('u');

subplot(3, 2, 4);
plot( x, T, '-*' );
title('T');

subplot(3, 2, 5);
plot( x - 0.5*dx, Ex, '-*' );
title('Ex');

eth = dx*sum( 0.5 * n .* T );
ek = dx*sum( 0.5 * n .* u.^2 );
ep = dx*sum( 0.5 * Ex.^2 );
t = [ t ; time ];
ekin = [ ekin ; ek ];
epot = [ epot ; ep ];
etherm = [ etherm ; eth ];
subplot(3, 2, 6);
plot(t, ekin, t, epot, t, ekin+epot+etherm, t, etherm);
legend('E_{kin}', 'E_{pot}', 'E_{tot}', 'E_{therm}', 'Location','best');

end
