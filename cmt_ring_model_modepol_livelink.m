%Created 5/29/22 for second and third set of calculations for James

% function [dat_beta,dat_neff,dat_Er,dat_Ez,dat_Ephi] = cmt_ring_model(wg_thickness,wg_width,ring_radius,domain_dis,swg_edge_outer,cmt_mesh_height,...
%             pml_width,index_m,f_sweep, n_sell_core,n_sell_clad,num_sol,n_w_arr,n_s_arr,n_height)
% function [dat_beta,dat_neff,dat_Er,dat_Ez,dat_Ephi,dat_Hr,dat_Hz,dat_Hphi] = cmt_ring_model(wg_thickness,wg_width,ring_radius,domain_dis,swg_edge_outer,cmt_mesh_height,...
%             pml_width,index_m,f_sweep, n_sell_core,n_sell_clad,num_sol,n_w_arr,n_s_arr,n_height,interp_coords)
function [dat_res,dat_res_off,ng,neff,mode_pol,mode_ind] = cmt_ring_model_rawfields_modepol_greg(wg_thickness,wg_width,ring_radius,domain_dis,...
            pml_width,index_m,f_sweep,freq_offset,n_sell_core,n_sell_clad,num_sol,n_w_arr,n_s_arr,n_height)

%Model creation - created by COMSOL

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

% model.modelPath('C:\Users\comsol\Desktop\Logan\AlGaAs Work\CLEO21 Figures and COMSOL text files');
model.modelPath('C:\Users\comsol\Desktop\Logan\AlGaAs Work\Sang-Yeon Ring Notches');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);
model.component('comp1').geom('geom1').axisymmetric(true);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');



%Studies
model.study.create('std1');
model.study('std1').create('mode', 'ModeAnalysis');
model.study('std1').feature('mode').set('shiftactive', false);
model.study('std1').feature('mode').set('conrad', '1');
model.study('std1').feature('mode').set('ngen', '1');
model.study('std1').feature('mode').set('ngenactive', false);
model.study('std1').feature('mode').activate('ewfd', true);
model.study('std1').feature('mode').set('transform', 'outofplanewavenumber');
model.study('std1').feature('mode').set('modeFreq', 'f_offset');
model.study('std1').feature('mode').set('neigsactive', true);
model.study('std1').feature('mode').set('neigs', num_sol);
model.study('std1').feature('mode').set('shiftactive', true);
model.study('std1').feature('mode').set('shift', 'index_m');

%Paramters
model.param.set('wg_thickness', wg_thickness);
model.param.set('wg_width', wg_width);
model.param.set('ring_radius', ring_radius);
model.param.set('domain_dis', domain_dis);
model.param.set('pml_width', pml_width);
model.param.set('index_m', index_m);
model.param.set('f_sweep', f_sweep);
model.param.set('f_offset', freq_offset);
model.param.set('n_sell_core', n_sell_core);
model.param.set('n_sell_clad', n_sell_clad);
% model.param.set('swg_edge_outer', swg_edge_outer);
% model.param.set('cmt_mesh_height', cmt_mesh_height);
%notch parameters
model.param.set('notch_thickness',n_height);

%if a structure would be created that is less than 5 nm in width (like a bump), make it go away
%basically, if n_w >= 495 nm, n_w = 500 nm
% for aa = 1:4
%     if n_w_arr(aa) >= 495e-9
%         n_w_arr(aa) = 500e-9;
%     end
% end

model.param.set('nw1',n_w_arr(1));
model.param.set('ns1',n_s_arr(1));
model.param.set('nw2',n_w_arr(2));
model.param.set('ns2',n_s_arr(2));
model.param.set('nw3',n_w_arr(3));
model.param.set('ns3',n_s_arr(3));
model.param.set('nw4',n_w_arr(4));
model.param.set('ns4',n_s_arr(4));


%Geometry
%%{
%Configuration with large square cladding around core, nothing else
%Core
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('base', 'corner');
model.component('comp1').geom('geom1').feature('r1').set('size', {'wg_width' 'wg_thickness'});
model.component('comp1').geom('geom1').feature('r1').set('pos', {'ring_radius - wg_width/2' '-wg_thickness/2'});
% model.component('comp1').geom('geom1').run('r1');
%Cladding
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('base', 'corner');
model.component('comp1').geom('geom1').feature('r2').set('size', {'2*domain_dis/2' 'domain_dis/2-wg_thickness/2'});
model.component('comp1').geom('geom1').feature('r2').set('pos', {'ring_radius - domain_dis/2' '-domain_dis/2'});
model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('r3').set('base', 'corner');
model.component('comp1').geom('geom1').feature('r3').set('size', {'2*domain_dis/2' 'domain_dis/2+wg_thickness/2'});
model.component('comp1').geom('geom1').feature('r3').set('pos', {'ring_radius - domain_dis/2' '-wg_thickness/2'});
%PML
%left bottom-half vert
model.component('comp1').geom('geom1').create('r4', 'Rectangle');
model.component('comp1').geom('geom1').feature('r4').set('base', 'corner');
model.component('comp1').geom('geom1').feature('r4').set('size', {'pml_width' '(domain_dis/2+pml_width-wg_thickness/2)'});
model.component('comp1').geom('geom1').feature('r4').set('pos', {'ring_radius-domain_dis/2-pml_width' '-(domain_dis/2+pml_width)'});
%right bottom-half vert
model.component('comp1').geom('geom1').create('r5', 'Rectangle');
model.component('comp1').geom('geom1').feature('r5').set('base', 'corner');
model.component('comp1').geom('geom1').feature('r5').set('size', {'pml_width' '(domain_dis/2+pml_width-wg_thickness/2)'});
model.component('comp1').geom('geom1').feature('r5').set('pos', {'ring_radius+domain_dis/2' '-(domain_dis/2+pml_width)'});
%left top-half vert
model.component('comp1').geom('geom1').create('r6', 'Rectangle');
model.component('comp1').geom('geom1').feature('r6').set('base', 'corner');
model.component('comp1').geom('geom1').feature('r6').set('size', {'pml_width' '(domain_dis/2+pml_width+wg_thickness/2)'});
model.component('comp1').geom('geom1').feature('r6').set('pos', {'ring_radius-domain_dis/2-pml_width' '-wg_thickness/2'});
%right top-half vert
model.component('comp1').geom('geom1').create('r7', 'Rectangle');
model.component('comp1').geom('geom1').feature('r7').set('base', 'corner');
model.component('comp1').geom('geom1').feature('r7').set('size', {'pml_width' '(domain_dis/2+pml_width+wg_thickness/2)'});
model.component('comp1').geom('geom1').feature('r7').set('pos', {'ring_radius+domain_dis/2' '-wg_thickness/2'});
%bottom horizontal
model.component('comp1').geom('geom1').create('r8', 'Rectangle');
model.component('comp1').geom('geom1').feature('r8').set('base', 'corner');
model.component('comp1').geom('geom1').feature('r8').set('size', {'2*(domain_dis/2+pml_width)' 'pml_width'});
model.component('comp1').geom('geom1').feature('r8').set('pos', {'ring_radius-domain_dis/2-pml_width' '-(domain_dis/2+pml_width)'});
%top horizontal
model.component('comp1').geom('geom1').create('r9', 'Rectangle');
model.component('comp1').geom('geom1').feature('r9').set('base', 'corner');
model.component('comp1').geom('geom1').feature('r9').set('size', {'2*(domain_dis/2+pml_width)' 'pml_width'});
model.component('comp1').geom('geom1').feature('r9').set('pos', {'ring_radius-domain_dis/2-pml_width' 'domain_dis/2'});


%create geometry
model.component('comp1').geom('geom1').run;
%%}

% figure
% mphgeom(model)


% figure
% mphgeom(model)
% for aa = 1:13
%     figure
%     mphgeom(model,'geom1','entity','domain','selection',[aa]);
%     title(num2str(aa))
% end


%Materials
% model.component('comp1').physics('ewfd').feature('wee1').set('DisplacementFieldModel', 'RefractiveIndex'); %not sure if needed
%{
%for case with smaller mesh box region
% model.component('comp1').material.create('mat1', 'Common');
% model.component('comp1').material('mat1').selection.set([7]);
% model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
% model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'n_sell_core'});
% model.component('comp1').material.create('mat2', 'Common');
% model.component('comp1').material('mat2').selection.set([1 2 3 4 5 6 8 9 10 11]);
% model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
% model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'n_sell_clad'});

%for larger mesh box region - core area is now 8
model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').selection.set([7]);
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'n_sell_core'});
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat2').selection.set([1 2 3 4 5 6 8 9 10]);
model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'n_sell_clad'});
%%}

%Boundary Conditions
model.component('comp1').coordSystem.create('pml1', 'PML');
model.component('comp1').coordSystem('pml1').selection.set([1 2 3 4 6 8 9 10]);
model.component('comp1').coordSystem('pml1').set('ScalingType', 'Cylindrical');
% model.component('comp1').coordSystem('pml1').selection.set([1 2 3 4 5 8 10 11 12 13]);


%Mesh
% model.component('comp1').mesh('mesh1').automatic(false);
% model.component('comp1').mesh('mesh1').feature('size').set('hmax', '1e-6');
%%{
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([7]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmax', '2e-8');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmin', '5e-9');
% model.component('comp1').mesh('mesh1').run('ftri1');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([5]);
% model.component('comp1').mesh('mesh1').run('ftri2');
model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('map1').selection.set([1 2 3 4 6 8 9 10]);
model.component('comp1').mesh('mesh1').feature('map1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('map1').feature('size1').set('hmax', 'pml_width / 5');
% model.component('comp1').mesh('mesh1').run('map1');
model.component('comp1').mesh('mesh1').run;
%}

%%{
%for case with smaller mesh box region
% model.component('comp1').material.create('mat1', 'Common');
% model.component('comp1').material('mat1').selection.set([7]);
% model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
% model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'n_sell_core'});
% model.component('comp1').material.create('mat2', 'Common');
% model.component('comp1').material('mat2').selection.set([1 2 3 4 5 6 8 9 10 11]);
% model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
% model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'n_sell_clad'});

%for larger mesh box region - core area is now 8
model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').selection.set([9]);
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'n_sell_core'});
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat2').selection.set([1 2 3 4 5 6 7 8 10 11 12 13]);
model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'n_sell_clad'});
%%}

%Boundary Conditions
model.component('comp1').coordSystem.create('pml1', 'PML');
model.component('comp1').coordSystem('pml1').selection.set([1 2 3 4 5 8 10 11 12 13]);
model.component('comp1').coordSystem('pml1').set('ScalingType', 'Cylindrical');
% model.component('comp1').coordSystem('pml1').selection.set([1 2 3 4 5 8 10 11 12 13]);


%Mesh
% model.component('comp1').mesh('mesh1').automatic(false);
% model.component('comp1').mesh('mesh1').feature('size').set('hmax', '1e-6');
%%{
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([9]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmax', '2e-8');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmin', '5e-9');
% model.component('comp1').mesh('mesh1').run('ftri1');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([6 7]);
% model.component('comp1').mesh('mesh1').run('ftri2');
model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('map1').selection.set([1 2 3 4 5 8 10 11 12 13]);
model.component('comp1').mesh('mesh1').feature('map1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('map1').feature('size1').set('hmax', 'pml_width / 5');
% model.component('comp1').mesh('mesh1').run('map1');
model.component('comp1').mesh('mesh1').run;

% figure
% mphmesh(model)



%Study - any changes/parameters that need to be set
model.sol.create('sol1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'mode');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'mode');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').set('neigs', num_sol);
model.sol('sol1').feature('e1').set('shift', '1');
model.sol('sol1').feature('e1').set('control', 'mode');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol('sol1').feature('e1').feature('d1').set('linsolver', 'mumps');
model.sol('sol1').feature('e1').feature('d1').label('Suggested Direct Solver (ewfd)');
model.sol('sol1').feature('e1').set('maxeigit', 10000);
model.sol('sol1').attach('std1');

% Save model
curdir=pwd;
% n = 1;
% FileName = ['\NotchtWG_LiveLink_',num2str(n),'.mph'];
% str_n = '011624_1';
str_currsave = strcat(string(datetime('today')),"_savenum.mat");
load(str_currsave,"save_num_check");
save_num = save_num_check;
str_n = strcat(string(datetime('today')),"_",num2str(save_num));
FileName2 = strcat('\NotchtWG_LiveLink_',str_n,'.mph');
model.save(strcat(curdir,FileName2));

RunFlag = 1;

if RunFlag

    model.sol('sol1').runAll;
    % model.save(strcat(curdir,FileName2));

%     curdir=pwd;
%     n = 5;
%     FileName = ['\NotchtWG_LiveLink_',num2str(n),'.mph'];
%     model.save(strcat(curdir,FileName));
    
    m = num_sol;

    eval_exp = {'ewfd.neff','ewfd.beta','ewfd.Er','ewfd.Ez','ewfd.Ephi',...
        'ewfd.Dr','ewfd.Dz','ewfd.Dphi','ewfd.Hr','ewfd.Hz','ewfd.Hphi'};
%     eval_exp = {'ewfd.neff','ewfd.beta','ewfd.Er','ewfd.Ez','ewfd.Ephi',...
%         'ewfd.Hr','ewfd.Hz','ewfd.Hphi'};
    dat_res_off = mpheval(model,eval_exp);

    if  abs(imag(dat_res_off.d2(m,1)*1e4)) > real(dat_res_off.d2(m,1))
        disp("spurious mode found in offset check")
    end

    

    %change frequency and run study again
    model.study('std1').feature('mode').set('modeFreq', 'f_sweep');

    model.sol('sol1').runAll;

    eval_exp = {'ewfd.neff','ewfd.beta','ewfd.Er','ewfd.Ez','ewfd.Ephi',...
        'ewfd.Dr','ewfd.Dz','ewfd.Dphi','ewfd.Hr','ewfd.Hz','ewfd.Hphi'};
%     eval_exp = {'ewfd.neff','ewfd.beta','ewfd.Er','ewfd.Ez','ewfd.Ephi',...
%         'ewfd.Hr','ewfd.Hz','ewfd.Hphi'};
    dat_res = mpheval(model,eval_exp);

    for b = 1:num_sol
        Er_sum(b) = sum(abs(dat_res.d3(b,:)));
        Ez_sum(b) = sum(abs(dat_res.d4(b,:)));

        r_rat(b) = Er_sum(b) / (Er_sum(b) + Ez_sum(b));
        z_rat(b) = Ez_sum(b) / (Er_sum(b) + Ez_sum(b));
    end
    
    %need to look at all modes that arent spurious - we can just look at
    %the first four for now (laziness)

    %find if modes is TE or TM
    TE_arr = zeros(1,num_sol);
    for b = 1:num_sol
        if r_rat(b) > 0.5
            TE_arr(b) = 1;
        end
    end
    
    mode_ind = num_sol-3:num_sol;

    mode_pol = TE_arr(mode_ind);

    for a = 1:length(mode_ind)


        mo = mode_ind(a);

        if abs(imag(dat_res.d2(mo,1)*1e4)) > real(dat_res.d2(mo,1))
            disp("spurious mode found")
        end



        %mode number to look at - primary mode = last mode number
        % m = curr_mode;

        neff(a) = abs(dat_res.d1(mo,1)) / ring_radius;
        neff_off = abs(dat_res_off.d1(mo,1)) / ring_radius;
        dn = (neff(a) - neff_off);
        dw = (f_sweep - freq_offset);
        group_ind_off = f_sweep * dn/dw;
        ng(a) = neff(a) + f_sweep * dn/dw;
        c_const = 2.9979e8;
        % Pow(a) = mphglobal(model,'ewfd.intWe+ewfd.intWm') * c_const/(2*pi*ring_radius)/ng(a);

    end



end %end RunFlag loop


end
% out = model;

