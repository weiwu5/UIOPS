function lwc = lwc_calc(conc,bins)

% calculates lwc in g/m3, needs conc in cm^-3 and bins in microns

%rho_i = 0.91;
lwc = sum(conc(:,:) .* pi/6 * (bins(:)*(1e-4)) .^3,2)*1e6;

