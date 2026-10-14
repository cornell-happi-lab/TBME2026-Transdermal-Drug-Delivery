function BSA_Drug_Delivery_Results_Cam_Ready
% Release analysis (Franz cell) for BSA from per-sample BSA concentrations.
% This script loads per-sample BSA concentrations (180 samples) from
% NNLS_Cam_Ready_results.csv, which must sit in the same folder as this script.
% Produces the following metrics to console:
% Q(afer 180 minutes), steady-state Jss / lag time (tau) / Papp per geometry x condition (PASS vs CIP).
% CIP denotes iontophoretic delivery; PASS denotes passive diffusion.

%% SETTINGS 
% Franz 
Vs = 1.0; % aliquot mL

% Diffusion area. Change to your setting.
A_cm2 = 175.49275 / 100;

% Donor concentration (hydrogel loading)
C0_ugmL = 20000;   % 20 mg/mL = 20000 µg/mL

%%  PER-SAMPLE NNLS RESULTS
% Columns: file, geom, t_min, isPass, Vr_mL, C_BSA_ugmL
% Vr_mL is the receptor volume, set by the adapter shape used:
% 10.255033 mL (Rec adapter) or 9.850898 mL (Cir adapter).
nnlsCsv = fullfile(fileparts(mfilename('fullpath')), 'NNLS_Cam_Ready_results.csv');
results = readtable(nnlsCsv, 'TextType','string', 'Delimiter',',');
results.isPass = logical(results.isPass);
results.geom(results.geom == "R") = "Plain";   % display name for geometry R
results = sortrows(results, {'geom','t_min'});

results3 = results;

if ~isempty(results3)
    %%  Q at 180 min
    target_t = 180;
    geoms_u = unique(results3.geom,'stable');
    cond_vals = [true false];   % PASS, CIP

    fprintf('\nQ (cumulative amount) at %d min (study endpoint): \n', target_t);
    for gi = 1:numel(geoms_u)
        geom_g = geoms_u(gi);
        for ci = 1:numel(cond_vals)
            isPass = cond_vals(ci);

            % take only this geom x condition
            Rsub = results3(results3.geom == geom_g & results3.isPass == isPass, :);
            if isempty(Rsub)
                continue;
            end

            % cumulative curve for this geom x condition
            Vref_g = median(Rsub.Vr_mL);
            [t_c, C_cum, E_cum] = per_time_cum_conc_sampling( ...
                Rsub.t_min, Rsub.C_BSA_ugmL, Rsub.Vr_mL, Vref_g, Vs);

            % convert to Q (µg/cm^2)
            Q_cum   = (Vref_g .* C_cum) ./ A_cm2;
            Q_err   = (Vref_g .* E_cum) ./ A_cm2;

            % find the point at 180 min
            idx = find(abs(t_c - target_t) < 1e-6, 1, 'first');

            if isempty(idx)
                % if for some reason this geom/cond has no 180 min, report closest
                [~, idx] = min(abs(t_c - target_t));
                note = ' (closest timepoint)';
            else
                note = '';
            end

            Q180     = Q_cum(idx);
            Q180_err = Q_err(idx);
            cond_str = ternary(isPass, "PASS", "CIP");

            fprintf('  %s | %s : Q(%d) = %.4g ± %.2g µg/cm^2%s\n', ...
                geom_g, cond_str, target_t, Q180, Q180_err, note);
        end
    end


    %% Steady-state, SD across replicate Franz cells 
    % Each cell gets its own Q(t) and its own steady state window 
    % Jss, tau, Papp and R^2 are reported as mean ± SD across the n = 3 cells  
    r2_thresh = 0.957;
    results3.cell = double(regexprep(results3.file, '^T\d+_(v\d+_)?B(\d+).*$', '$2', 'ignorecase'));

    fprintf('\nSteady-state parameters across replicate Franz cells (mean ± SD, each cell fitted on its own steady-state window):\n');
    for gi = 1:numel(geoms_u)
        geom_g = geoms_u(gi);
        for ci = 1:numel(cond_vals)
            isPass = cond_vals(ci);
            Rsub = results3(results3.geom==geom_g & results3.isPass==isPass,:);
            if isempty(Rsub), continue; end
            cond_str = ternary(isPass,'PASS','CIP');

            cells = unique(Rsub.cell);
            Jss_k = nan(numel(cells),1);
            tau_k = nan(numel(cells),1);
            r2_k  = nan(numel(cells),1);
            win_k = strings(numel(cells),1);
            for k = 1:numel(cells)
                Rk = Rsub(Rsub.cell == cells(k), :);
                Vk = median(Rk.Vr_mL);
                [t_k, C_k] = per_time_cum_conc_sampling(Rk.t_min, Rk.C_BSA_ugmL, Rk.Vr_mL, Vk, Vs);
                Q_k = (Vk .* C_k) ./ A_cm2;
                bk = pick_ss_window_longest(t_k, Q_k, r2_thresh, 2);
                if isnan(bk.slope), win_k(k) = "none"; continue; end
                Jss_k(k) = bk.slope * 60;   % µg cm^-2 h^-1
                tau_k(k) = max(-bk.int / bk.slope, 0);   % lag time can't be negative -> clip at 0
                r2_k(k)  = bk.r2;
                win_k(k) = sprintf('%g–%g', t_k(bk.i1), t_k(bk.i2));
            end
            Papp_k = (Jss_k / C0_ugmL) / 3600; % cm/s

            fprintf('[Results] %s | %s : n=%d cells (windows %s), Jss=%.4g ± %.4g (µg/cm²/h), τ=%.3g ± %.3g min, Papp=%.3g ± %.3g (cm/s), R²=%.4f ± %.4f\n', ...
                geom_g, cond_str, nnz(~isnan(Jss_k)), strjoin(win_k, ', '), ...
                mean(Jss_k,'omitnan'),  std(Jss_k,'omitnan'), ...
                mean(tau_k,'omitnan'),  std(tau_k,'omitnan'), ...
                mean(Papp_k,'omitnan'), std(Papp_k,'omitnan'), ...
                mean(r2_k,'omitnan'),   std(r2_k,'omitnan'));
        end
    end

    %% Steady-state analysis and SE based on cumulative Q(t) 
    fprintf('\nDerived steady-state parameters from cumulative Q(t). (± on slope, Jss, τ and Papp are standard errors; see the Supplementary Information file.):\n');
    
    for gi = 1:numel(geoms_u)
        geom_g = geoms_u(gi);
        for ci = 1:numel(cond_vals)
            isPass = cond_vals(ci);
            Rsub = results3(results3.geom==geom_g & results3.isPass==isPass,:);
            if isempty(Rsub), continue; end
    
            % Build cumulative Q curve
            Vref_g = median(Rsub.Vr_mL);
            [t_c, C_cum] = per_time_cum_conc_sampling(Rsub.t_min, Rsub.C_BSA_ugmL, Rsub.Vr_mL, Vref_g, Vs);
            Q_cum = (Vref_g .* C_cum) ./ A_cm2;
    
            % Best linear region (>=3 pts, R2 > r2_thresh)
            best = pick_ss_window(t_c, Q_cum, r2_thresh);
    
            if isnan(best.slope)
                fprintf('[Results] %s | %s : no region with R² > %.3f\n', ...
                        geom_g, ternary(isPass,'PASS','CIP'), r2_thresh);
                continue;
            end
    
            % Compute derived quantities
            seg_t = t_c(best.i1:best.i2);
            dof = numel(seg_t)-2;
            R=best.S.R; Rinverse=R\eye(size(R)); df=dof;
            covp=(Rinverse*Rinverse.')*(best.S.normr^2/df);
            slope_se=sqrt(covp(1,1)); int_se=sqrt(covp(2,2));
            cov_si=covp(1,2);
    
            Jss = best.slope*60; % µg cm^-2 h^-1
            Jss_se = slope_se*60;
            tau = max(-best.int/best.slope, 0);  % lag time can't be negative -> clip at 0
            dtaudb=-1/best.slope; dtaudm=best.int/(best.slope^2);
            tau_se=sqrt(max((dtaudb^2)*int_se^2+(dtaudm^2)*slope_se^2+ ...
                            2*dtaudb*dtaudm*cov_si,0));
            Papp = (Jss / C0_ugmL) / 3600;
            Papp_se = (Jss_se / C0_ugmL) / 3600;
    
            fprintf('[Results] %s | %s : %.0f–%.0f min, slope=%.3g ± %.1g, Jss=%.3g ± %.1g (µg/cm²/h), τ=%.2f ± %.2f min, Papp=%.2e ± %.1e (cm/s), R²=%.3f\n', ...
                geom_g, ternary(isPass,'PASS','CIP'), ...
                t_c(best.i1), t_c(best.i2), ...
                best.slope, slope_se, ...
                Jss, Jss_se, tau, tau_se, ...
                Papp, Papp_se, best.r2);
        end
    end


end

end % end main

%% HELPERS 
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function best = pick_ss_window_longest(t_c, Q_cum, r2_thresh, minPts)
% Longest window of >= minPts consecutive points with R^2 >= r2_thresh
% (ties -> higher R2). Used per replicate cell. best.slope is NaN if none qualifies.
    n = numel(t_c);
    best = struct('i1',NaN,'i2',NaN,'r2',-Inf,'slope',NaN,'int',NaN);
    for i1 = 1:n-minPts+1
        for i2 = i1+minPts-1:n
            tt = t_c(i1:i2); yy = Q_cum(i1:i2);
            p  = polyfit(tt, yy, 1);
            yhat = polyval(p, tt);
            r2 = 1 - sum((yy-yhat).^2) / (sum((yy-mean(yy)).^2) + eps);
            L  = i2 - i1;
            Lb = best.i2 - best.i1; if isnan(Lb), Lb = -1; end
            if r2 >= r2_thresh && (L > Lb || (L == Lb && r2 > best.r2))
                best = struct('i1',i1,'i2',i2,'r2',r2,'slope',p(1),'int',p(2));
            end
        end
    end
end

function best = pick_ss_window(t_c, Q_cum, r2_thresh)
% Scan all windows of >= 3 consecutive points; keep the highest R² above
% r2_thresh (ties -> longer window). best.slope is NaN if none qualifies.
    n = numel(t_c);
    best = struct('i1',NaN,'i2',NaN,'r2',-Inf,'slope',NaN,'int',NaN,'S',[]);
    for i1=1:n-2
        for i2=i1+2:n
            tt=t_c(i1:i2); yy=Q_cum(i1:i2);
            [p,S]=polyfit(tt,yy,1);
            yhat=polyval(p,tt);
            r2=1-sum((yy-yhat).^2)/(sum((yy-mean(yy)).^2)+eps);
            if r2>r2_thresh && (r2>best.r2 || (r2==best.r2 && (i2-i1)>best.i2-best.i1))
                best=struct('i1',i1,'i2',i2,'r2',r2,'slope',p(1),'int',p(2),'S',S);
            end
        end
    end
end

function [t_sorted, Ccum_mean, Ccum_std] = per_time_cum_conc_sampling(t, C, Vr_mL, Vref, Vs)
    if nargin < 5 || isempty(Vs), Vs = 0; end
    if nargin < 4 || isempty(Vref)
        Vref = median(Vr_mL(~isnan(Vr_mL)));
    end

    T = table(t(:), C(:), Vr_mL(:), 'VariableNames', {'t','C','Vr'});
    [t_u,~,idx] = unique(T.t);
    C_mean = accumarray(idx, T.C, [], @(v) mean(v, 'omitnan'));
    C_std  = accumarray(idx, T.C, [], @(v) std(v,  'omitnan'));

    [t_sorted, order] = sort(t_u, 'ascend');
    C_mean = C_mean(order);
    C_std  = C_std(order);

    n = numel(t_sorted);
    C_var = C_std.^2;
    Mcum_mean = zeros(n,1);
    Mcum_var  = zeros(n,1);

    for j = 1:n
        Mcum_mean(j) = Vref * C_mean(j) + Vs * sum(C_mean(1:j-1));
        Mcum_var(j)  = (Vref^2) * C_var(j) + (Vs^2) * sum(C_var(1:j-1));
    end

    Ccum_mean = Mcum_mean ./ Vref;
    Ccum_std  = sqrt(Mcum_var) ./ Vref;
end
