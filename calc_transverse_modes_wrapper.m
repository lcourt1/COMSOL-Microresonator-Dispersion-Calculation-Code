

%save file inputs - need to input current save number and save number for
%centered frequency data
%increments automatically so each iteration is saved with a different
%number
reset_save_num = 0; %set to 0 - dont reset, not set to 0 - reset (use 1?)
str_currsave = strcat(string(datetime('today')),"_savenum.mat");
if isfile(str_currsave) && reset_save_num == 0
    load(str_currsave);
    save_num = save_num_check + 1;
    save_num_check = save_num;
    save(str_currsave,"save_num_check")
else
    save_num = 1;
    save_num_check = save_num;
    save(str_currsave,"save_num_check");
end

%string name of centered frequency data
% str_cenfreq = "transverse_mode_data_011624_4.mat";


%ring resonator parameters
ring_thickness = 687e-9;
ring_width = 802e-9;
% wg_width = 6e-6;
ring_radius = 23e-6;
domain_dis = 10e-6; %size of one side of entire model domain (square)
pml_width = 1e-6; %on each side of domain
m_x = 0.2; %AlGaAs mole fraction
index_m_wg = 4.5e2; %starting mode number guess for notched design (1.1e3 for bowers algaas design at 1550)
% index_m_wg = 1.2e3;
%need to update the index_m guess at higher frqeuencies - should atleast be
%1500 at 2.5e14 Hz

% index_m_wg_s = 1.1e3;
% index_m_wg_mid = 1.3e3;
% index_m_wg_e = 1.5e3;
index_m_wg_s = 1.3e2;
% index_m_wg_mid = 1.1e3;
index_m_wg_e = 4.5e2;

n_period = 0.5e-6; %general width of notches
nw1 = 1*n_period;
ns1 = 8*n_period;
nw2 = 1*n_period;
ns2 = 9*n_period;
nw3 = 0.5*n_period;
ns3 = 10*n_period;
nw4 = 0.1*n_period;
ns4 = 11*n_period;
%notch width and starting position vectors for passing into function
n_w_arr = [nw1;nw2;nw3;nw4];
n_s_arr = [ns1;ns2;ns3;ns4];
% n_w_arr = [0 0 0 0]; %starting with no notches (4 spots)
% n_s_arr = [0 0 0 0]; %notch positions
% n_height = 5e-8; %notch height
n_height = 0;

% index_m_wg = 1e3; %starting mode number guess
% n_w_arr = [0 0 0 0]; %starting with no notches (4 spots)
% n_s_arr = [0 0 0 0]; %notch positions
% n_height = 5e-8; %notch height


%frequencies considered - starting with just 1550 nm
lmd_cen = 1.55e-6;
c_const = 2.9979e8;
freq_cen = c_const / lmd_cen; %1.55 um
freq_cen = 2.82e14;
% dfreq = 5e12;
% freq_arr = f_sweep;
% freq_l = 1.85e14; %roughly 1.62 um (+7)
freq_l = 1.6e14;
freq_midtolow = freq_cen - freq_l;
freq_h = freq_cen + freq_midtolow;
num_freq = 121;
dfreq = (freq_h - freq_l) / (num_freq-1);
freq_arr = freq_l:dfreq:freq_h;
%frequency offset for calculating group index
freq_diff = 1e9;
freq_offset = freq_arr - freq_diff;

%array for index_m starting guesses
dind = (index_m_wg_e - index_m_wg_s) / (num_freq-1);
index_m_wg_arr = index_m_wg_s:dind:index_m_wg_e;

h_const = 6.6261e-34;
c_const = 2.9979e8;
x = 0.2;
ED0 = (1.765 + 1.115*x + 0.37*(x^2))*1.60218e-19;
E0 = (1.425 + 1.155*x + 0.37*(x^2))*1.60218e-19;

for a  = 1:num_freq
    lmd = c_const / freq_arr(a);
    %index of refraction data
    % n_sell_AlGaAs = 3.80448 - 403388*lmd; %fitted from ARL data
    %using Sellmeier equation to see if results better match bowers paper
    %(for AlGaAs)
    kai0 = h_const*c_const / (lmd*ED0);
    kai = h_const*c_const / (lmd*E0);
    f_kai0 = (2 - sqrt(1+kai0) - sqrt(1-kai0))/(kai0^2);
    f_kai = (2 - sqrt(1+kai) - sqrt(1-kai))/(kai^2);
    n_sell_AlGaAs = sqrt((6.3 + 19*x)*(f_kai+(f_kai0/2)*(E0/ED0)^1.5) + (9.4 - 10.2*x));

    %Sellmeier equations for SiO2
    n_sell_SiO2 = sqrt(((0.6961663*((lmd*1e6)^2))/(((lmd*1e6)^2)-(0.0684043^2)))+...
        ((0.4079426*((lmd*1e6)^2))/(((lmd*1e6)^2)-(0.1162414^2)))+...
        ((0.8974794*((lmd*1e6)^2))/(((lmd*1e6)^2)-(9.896161^2)))+1);

    %Sellmeier equation for SiN from Greg's COMSOL code: CouplingCMT.m
    % n_SiN_71_Boulder = @(x) sqrt(1 + 2.906 * x^2/(x^2 - 0.13782^2)  - 0.01003*x^2);
    ln = lmd*1e6; %normalized lambda
    n_sell_SiN_71_Boulder(a) = sqrt(1 + 2.906 * ln^2/(ln^2 - 0.13782^2)  - 0.01003*ln^2);

    n_sell_core(a) = n_sell_SiN_71_Boulder(a);
    n_sell_clad(a) = n_sell_SiO2;
    % n_sell_core = 1.976; %from gregs code N_sin_bolder
    % n_sell_clad = 1;
end

%straight waveguide parameters
swg_thickness = wg_thickness; %same thickness as resonator
% swg_width = wg_width; %same width as resonator
% swg_thickness = 400e-9; %nm
% swg_width = 700e-9; %nm
swg_width = 1e-6;
% index_m_swg = 1e7; %starting mode number guess
if swg_width > swg_thickness
    m = 1;
    n = 0;
    mu = 4*pi*1e-7; %mu = mu0
    epsilon = n_sell_core^2 * 8.854e-12;
    index_m_swg = sqrt((2*pi*f_sweep)^2*mu*epsilon - (m*pi/swg_width)^2 - (n*pi/swg_thickness)^2);
else
    m = 0;
    n = 1;
    mu = 4*pi*1e-7; %mu = mu0
    epsilon = n_sell_core^2 * 8.854e-12;
    index_m_swg = sqrt((2*pi*f_sweep)^2*mu*epsilon - (m*pi/swg_width)^2 - (n*pi/swg_thickness)^2);
end


domain_dis = 6e-6;
cmt_mesh_height = 1.5*ring_thickness;
%other system parameters
num_sol = 50; %number of mode solutions to search for

%Get comsol data
% tic
% %resonator model function call
% index_m_wg = index_m_wg_s;
for a = 1:num_freq
    % if a > num_freq / 3
    %     if a > num_freq * 2 / 3
    %         index_m_wg = index_m_wg_e;
    %     else
    %         index_m_wg = index_m_wg_mid;
    %     end
    % end
    index_m_wg = index_m_wg_arr(a);
    index_m_wg
    freq_arr(a)
    if a == 10
        bbb = 3;
    end
    tic
    [dat_ring_res(a),dat_ring_res_off(a),ring_ng(a,:),ring_neff(a,:),mode_pol(a,:),mode_ind(a,:)] = cmt_ring_model_rawfields_modepol_greg(ring_thickness,ring_width,ring_radius,domain_dis,...
        pml_width,index_m_wg,freq_arr(a),freq_offset(a),n_sell_core(a),n_sell_clad(a),num_sol,n_w_arr,n_s_arr,n_height);
    time_ring = toc
end

str_datsave = strcat(string(datetime('today')),"_transverse_mode_data_iter",num2str(save_num),".mat");
% save("transverse_mode_data_011624_5.mat","dat_ring_res","dat_ring_res_off","ring_ng","ring_neff","mode_pol","mode_ind","ring_thickness","ring_width","ring_radius","freq_arr","n_sell_clad","n_sell_core","freq_diff");
save(str_datsave,"dat_ring_res","dat_ring_res_off","ring_ng","ring_neff","mode_pol","mode_ind","ring_thickness","ring_width","ring_radius","freq_arr","n_sell_clad","n_sell_core","freq_diff");



%% Disperison calculations

% load('transverse_mode_data_011624_5.mat')


%primarily use this
str_currsave = strcat(string(datetime('today')),"_savenum.mat");
load(str_currsave);
save_num = save_num_check;
str_datload = strcat(string(datetime('today')),"_transverse_mode_data_iter",num2str(save_num),".mat");
load(str_datload);

%use this for yesterdays data
% str_currsave = strcat(string(datetime('yesterday')),"_savenum.mat");
% load(str_currsave);
% save_num = save_num_check;
% str_datload = strcat(string(datetime('yesterday')),"_transverse_mode_data_iter",num2str(save_num),".mat");
% load(str_datload);


% str_datload = str_datsave;
% str_datload = strcat(string(datetime('today')),"_transverse_mode_data_iter",num2str(save_num),".mat");
% str_datload = "transverse_mode_data_011624_5.mat";
% load(str_datload);


mode_check = mode_ind(1,:); %all the same for now
num_freq = length(freq_arr);
cen = ceil(num_freq / 2);
%currently TE mode
% mu_r = freq_mu_TE;
%     %currently TM mode
%     mu_r = freq_mu_TM;

%calculate FSR for last four modes
% w_r = 2*pi * freq_arr;
% w0 = w_r(cen);
% w_rel = w_r - w0;
% for b = 1:length(mode_check)
%     m = mode_check(b);
%     for a = 1:num_freq
%         mu_t_on = dat_ring_res(a).d2(m,1);
%         % mu_t_off = dat_ring_res_off(a).d2(m,1);
%         mu_r(a) = mu_t_on;
%     end
%     mu0 = mu_r(cen);
%     mu_rel = mu_r - mu0;
%     %fitting using rise over run between main and offset frequency mu's
%     mu_c_on = dat_ring_res(cen).d2(m,1);
%     mu_c_off = dat_ring_res_off(cen).d2(m,1);
%     slope_c = abs(mu_c_on - mu_c_off) / freq_diff;
%     D1_slope(b) = (2*pi) / slope_c;
%     %fitting first derivative using polyfitn
%     % D_poly = polyfit(mu_rel,w_rel,2);
%     D1_poly2 = polyfitn(mu_rel,w_rel,'x x^2');
%     % D1_fit2(b,:) =  1 / D1_poly2.Coefficients;
% end

w_r = 2*pi * freq_arr;
w0 = w_r(cen);
w_rel = w_r - w0;
% chosen_ind = 4;
num_modes = length(mode_check);




%calculate indices for each mode
%TE1 - first 1
%TM1 - first 0
%TE0 - second 1
%TM0 - second zero
%index order in matrix - TM1, TE1, TM0, TE0
%this is very limited - only works if both TM1 and TE1 are present

%add if statement based on summation - if sums to more than 2, than assume
%there is a spurious mode, which means I should only look for the first two
%modes, which I assume will always be there

sum_ind = sum(mode_pol(1,:));
if sum_ind > 2
    %looking for two modes only
    num_modes = 2;
    mode_ind_sort = zeros(num_freq,num_modes);
    mode_check = mode_check(end-1:end);
    for b = 1:num_modes
        for a = 1:num_freq
            if b == 1 %TM0
                for c = 1:num_modes
                    if mode_pol(a,end-c+1) == 0
                        mode_ind_sort(a,1) = mode_check(end-c+1);
                        break
                    end
                end
            else %TE0
                for c = 1:num_modes
                    if mode_pol(a,end-c+1) == 1
                        mode_ind_sort(a,2) = mode_check(end-c+1);
                        break
                    end
                end
            end
        end
    end
else
    %looking for 4 modes
    mode_ind_sort = zeros(num_freq,num_modes);
    for b = 1:num_modes
        for a = 1:num_freq
            if b == 1 %TM1
                for c = 1:num_modes
                    if mode_pol(a,c) == 0
                        mode_ind_sort(a,1) = mode_check(c);
                        break
                    end
                end
            elseif b == 2 %TE1
                for c = 1:num_modes
                    if mode_pol(a,c) == 1
                        mode_ind_sort(a,2) = mode_check(c);
                        break
                    end
                end
            elseif b == 3 %TM0
                found_zero = 0;
                for c = 1:num_modes
                    if mode_pol(a,c) == 0 && ~found_zero
                        found_zero = 1;
                    elseif mode_pol(a,c) == 0 && found_zero
                        mode_ind_sort(a,3) = mode_check(c);
                        break
                    end
                end
            else %b == 4, TE0
                found_zero = 0;
                for c = 1:num_modes
                    if mode_pol(a,c) == 1 && ~found_zero
                        found_zero = 1;
                    elseif mode_pol(a,c) == 1 && found_zero
                        mode_ind_sort(a,4) = mode_check(c);
                        break
                    end
                end
            end
        end
    end
end

%following same pattern - TM1, TE1, TM0, TE0
for b = 1:length(mode_check)
    % m = mode_check(b);
    for a = 1:num_freq
        m = mode_ind_sort(a,b);
        mu_t_on = dat_ring_res(a).d2(m,1);
        % mu_t_off = dat_ring_res_off(a).d2(m,1);
        mu_r(a) = mu_t_on;
    end
    mu0 = mu_r(cen);
    mu_r_save(b,:) = mu_r;
    mu_rel = mu_r - mu0;
    mu_rel_save(b,:) = mu_rel;
    %fitting using rise over run between main and offset frequency mu's
    mu_c_on = dat_ring_res(cen).d2(m,1);
    mu_c_off = dat_ring_res_off(cen).d2(m,1);
    slope_c = abs(mu_c_on - mu_c_off) / freq_diff;
    D1_slope(b) = (2*pi) / slope_c;
    %fitting first derivative using polyfitn
    D1_poly2 = polyfitn(w_rel,mu_rel,'x');
    D1_fit2(b,:) =  1 / D1_poly2.Coefficients;
    % D1 = D1_fit2;
    D1 = D1_slope;

    %Calculate D2
    D_int(b,:) = (w_rel - mu_rel * D1(b)) / (2*pi);
    % D_int(b,:) = (w_rel - mu_rel * D1_TE0) / (2*pi);

    % figure
    % plot(mu_rel,D_int(b,:))

    %calculate D2 / 2pi
    D2_poly = polyfitn(mu_rel,D_int(b,:),'x^2');
    D2(b) = 2 * D2_poly.Coefficients;
end
D1 = D1_fit2;
D1 / (2*pi)

% figure
% plot(freq_arr,mu_r_save)
% 
% figure
% plot(mu_r_save,freq_arr)

aa = 1;


%calculate chosen mode D1 from frequency data located around central freq
%{
str_censave = "transverse_mode_data_011624_4.mat";
str_D1save = strcat(string(datetime('today')),"calcD1_for_lambda1550nm_centered_iter",num2str(save_num),".mat");
if ~isfile(str_D1save)
    calc_D1_chosenmode(str_censave,str_D1save,chosen_ind);
end
load(str_D1save);
%}
% [freq_arr_cen,mu_chosen_cen,D1_chosen] = calc_D1_chosenmode(str_censave,str_D1save,chosen_ind);

% load(str_censave);
% freq_arr_cen = freq_arr;
% mu_chosen_cen = mu_r_save(chosen_ind,:);
% D1_chosen = D1(chosen_ind);
% str_D1save = strcat(string(datetime('today')),"calcD1_for_lambda1550nm_centered_iter",num2str(save_num),".mat");
% save(str_D1save,"freq_arr_cen","D1_chosen","mu_chosen_cen");

% str_datload = str_datsave;
% str_datload = strcat(string(datetime('today')),"_transverse_mode_data_iter",num2str(save_num),".mat");
% str_datload = "transverse_mode_data_011624_5.mat";
% load(str_datload);



% D1_TE0 = D1(chosen_ind);
% mu_TE0 = mu_rel_save;

% %calculate Dint for last four modes
% for b = 1:length(mode_check)
%     % m = mode_check(b);
%     for a = 1:num_freq
%         m = mode_ind_sort(a,b);
%         mu_t_on = dat_ring_res(a).d2(m,1);
%         % mu_t_off = dat_ring_res_off(a).d2(m,1);
%         mu_r(a) = mu_t_on;
%     end
%     % mu_r_save2(b,:) = mu_r;
%     mu0 = mu_r(cen);
%     if b == chosen_ind
%         mu0_save = mu0;
%     end
%     mu_rel = mu_r - mu0;
%     % mu_rel_save2(b,:) = mu_rel;
%     D_int(b,:) = (w_rel - mu_rel * D1(b)) / (2*pi);
%     % D_int(b,:) = (w_rel - mu_rel * D1_TE0) / (2*pi);
% 
%     % figure
%     % plot(mu_rel,D_int(b,:))
% 
%     %calculate D2 / 2pi
%     D2_poly = polyfitn(mu_rel,D_int(b,:),'x^2');
%     D2(b) = 2 * D2_poly.Coefficients;
% end

%plot for checking D2 fit
% for b = 1:num_modes
%     figure
%     plot(mu_rel_save(b,:),D_int(b,:),'-b','Linewidth',3)
%     hold on
%     plot(mu_rel_save(b,:),0.5*D2(b)*mu_rel_save(b,:).^2,'--r','Linewidth',2)
%     hold off
% end

aaa =3;

%{
for b = 1:length(mode_check)
    m = mode_check(b);
    for a = 1:num_freq
        mu_t_on = dat_ring_res(a).d2(m,1);
        % mu_t_off = dat_ring_res_off(a).d2(m,1);
        mu_r(a) = mu_t_on;
    end
    mu0 = mu_r(cen);
    mu_rel = mu_r - mu0;
    %fitting using rise over run between main and offset frequency mu's
    mu_c_on = dat_ring_res(cen).d2(m,1);
    mu_c_off = dat_ring_res_off(cen).d2(m,1);
    slope_c = abs(mu_c_on - mu_c_off) / freq_diff;
    D1_slope(b) = (2*pi) / slope_c;
    %fitting first derivative using polyfitn
    D1_poly2 = polyfitn(w_rel,mu_rel,'x');
    D1_fit2(b,:) =  1 / D1_poly2.Coefficients;
end
D1 = D1_fit2;
D1 / (2*pi)
D1_TE0 = D1(chosen_ind);

%calculate Dint for last four modes
for b = 1:length(mode_check)
    m = mode_check(b);
    for a = 1:num_freq
        mu_t_on = dat_ring_res(a).d2(m,1);
        % mu_t_off = dat_ring_res_off(a).d2(m,1);
        mu_r(a) = mu_t_on;
    end
    mu_r_save(b,:) = mu_r;
    mu0 = mu_r(cen)
    if b == chosen_ind
        mu0_save = mu0;
    end
    mu_rel = mu_r - mu0;
    mu_rel_save(b,:) = mu_rel;
    % D_int(b,:) = (w_rel - mu_rel * D1(b)) / (2*pi);
    D_int(b,:) = (w_rel - mu_rel * D1_TE0) / (2*pi);

    % figure
    % plot(mu_rel,D_int(b,:))

    %calculate D2 / 2pi
    D2_poly = polyfitn(mu_rel,D_int(b,:),'x^2');
    D2(b) = 2 * D2_poly.Coefficients;
end
%}


%need to fix this - due to the inversion, the mu's dont match - mu here
%refers to the mode that I'm primarily looking at, which means I need to
%calculate the range of frequencies that allows the other modes to match
%the absolute mode value (m) of the primary mode (which is a lot in the
%case of TE0)

%alternatively, can try to set up an inversion where there is a cross on
%the x axis beforehand, which would then translate to a proper avoided
%crossing

chosen_ind = num_modes;

%calculating at chosen frequency - need to calculate D1 for this frequency
%before running this part, then need to load it
% freq_arr_cen = freq_arr;
% mu_TE0_cen = mu_r_save(chosen_ind,:);
% save("calcD1_for_lambda1550nm_centered_011624_1.mat","freq_arr_cen","D1_TE0","mu_TE0_cen")
% load calcD1_for_lambda1550nm_centered_011624_1.mat;
% D1_chosen = 
% cen_chosen = ceil(length(freq_arr_cen) / 2);
% w_cen = 2*pi * freq_arr_cen(cen_chosen);
% w_rel_cen = 2*pi * freq_arr - w_cen;
% mu0 = mu_chosen_cen(1,cen_chosen);
w_rel_cen = w_rel;
mu0 = mu_r_save(chosen_ind,cen);
D1_use = D1(chosen_ind);
% D1_use = 7.9567e11;
% mu0 = 1.156e3;
for b = 1:length(mode_check)
    % mu_r_other0(b,:) = mu_r_save(b,:) - mu_r_save(b,cen);
    % D_mulint(b,:) = (w_rel_cen - mu_r_other0(b,:) * D1(b)) / (2*pi);
    mu_r_other0(b,:) = mu_r_save(b,:) - mu0;
    D_mulint(b,:) = (w_rel_cen - mu_r_other0(b,:) * D1_use) / (2*pi);
end
figure
plot(mu_r_other0(1,:),D_mulint(1,:))
hold on
for a = 2:length(mode_check)
    plot(mu_r_other0(a,:),D_mulint(a,:))
end
hold off
% ylim([0 max(D_mulint(chosen_ind,:))])
% xlim([-50 50])
% ylim([0 1e10])



% %already calculated without 2pi
% D2
% 
% figure
% plot(mu_rel,D_int(1,:))
% hold on
% for a = 2:length(mode_check)
%     plot(mu_rel,D_int(a,:))
% end
% hold off
% ylim([0 max(D_int(end,:))])
% 
% D2_t = D2(end-3);
% D2_bowers = 10.8e6;
% mu = -20:20;
% disp_t = D2_t * 0.5 * mu.^2;
% disp_bowers = D2_bowers * 0.5 * mu.^2;

% figure
% plot(mu,disp_bowers,'Linewidth',2)
% hold on
% plot(mu,disp_t,'--','Linewidth',3)
% hold off
% xlabel("Mode number $\mu$",'interpreter','latex')
% ylabel("$D_\textrm{int}$",'interpreter','latex')
% legend('Bowers paper result','Calculated result','interpreter','latex')



% save("dispcalc_from_transverse_mode_data_011624_5.mat","D1","mu_r_save","D_int","D2")