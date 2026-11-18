clc; clear; close all;
% =====================================================================================================
% VERSION 3 - changes with respect to version 2
%  [FIX 1] SMC sliding surface uses a NORMALISED incremental-conductance error e = (dI/dV + I/V)/(I/V).
%          The old error scaled with irradiance, so one fixed boundary layer phi was too slow at low G
%          (stations 00, 01, 08, 09) and too aggressive at high G. The normalised error is irradiance independent.
%          The same error is used by "SMC only" and "GA-SMC", so the ablation stays fair (only the GA tuning differs).
%  [FIX 2] GA is trained on a MULTI-CONDITION set (steps, rapid cloud transition, temperature ramp and real days of
%          training stations mp.trStations) instead of one artificial profile. Fitness = mean objective + worst-profile
%          tracking error (robustness) + ripple penalty. Search bounds are widened and printed with at-bound flags.
%  [FIX 3] Honest evaluation: stations in mp.trStations are used for GA training only; held-out mean efficiency is printed.
%          Scenario D-F stations (3, 7, 10) are never used for training.
%  [FIX 4] Statistical interpretation distinguishes "lower/higher mean, not significant" from "no difference".
%  [FIX 5] Notes sheet / console claims are generated from the actual numbers of the run (no pre-written claims).
%  Earlier fixes kept: bound-safe PV calibration (gT, beta), horizontal tick labels, scaling-chain notes.
% =====================================================================================================

%% SETTINGS
dataFolder = 'E:\Satheesh\january\November\2026-11-KIT-COC-ST-091\data';
displayStation = 1;            % plots are shown only for this station (1 = station00 ... 10 = station09)
scenStations = [3 7 10];       % real profiles for Scenario D, E, F (1-based)
daylightThreshold = 20;
fontName = 'Times New Roman'; fontSize = 22; fontWt = 'bold';

% ---- STEP 5 SETTINGS ----
capacityMW    = NaN(1,10);     % installed AC capacity per station (MW); NaN = auto-estimate (99.9 percentile, rounded up to 0.5 MW)
dcAcRatio     = 1.10;          % only used for the NOMINAL (uncalibrated) array size
etaSys        = 0.97;          % FIXED system efficiency (inverter + wiring + soiling + mismatch)
testDayMod    = [3 6 9];       % days with mod(dayIndex,10) in this set form the TEST set (30 %), the rest TRAIN (70 %)
cc = struct();                 % calibration settings
cc.gT   = [0.5 1.5];  cc.gT0 = 0.75; cc.sgT = 0.25;     % thermal factor: bounds, physical prior value, prior std
cc.beta = [0.8 1.2];  cc.b0  = 1.00; cc.sb  = 0.10;     % irradiance-response exponent: bounds, prior value, prior std
cc.w    = 0.01;                                          % prior weight relative to the normalised squared misfit

% ---- STEP 9-18 SETTINGS ----
mp = struct();
mp.nSer = 2; mp.nPar = 2;
mp.L = 2e-3; mp.rL = 0.05; mp.Cin = 220e-6; mp.fsw = 20e3;
mp.Vbat = 120;
mp.Dmin = 0.05; mp.Dmax = 0.90; mp.D0 = 0.80;
mp.dt = 1e-4;
mp.TsSMC = 1e-3; mp.TsPO = 5e-3;
mp.stepD = 0.005; mp.incEps = 0.005;
mp.smcDefault = [0.002 5 0.05];                          % untuned reference SMC [lambda K phi] (normalised error)
mp.gaPop = 24; mp.gaGen = 18; mp.gaW = [0.4 0.2 0.4]; mp.gaSeed = 1; mp.gaTsim = 1.2;
mp.rippleCap = 1.15;           % GA-SMC ripple may exceed the lowest reference ripple (P&O / INC / default SMC) by at most 15 %
mp.ripPen = 5;                 % penalty weight for violating the ripple constraint
mp.lamMax = 0.02; mp.Krange = [0.5 100]; mp.phiRange = [1e-3 0.5];
mp.trStations = [1 2 9];       % real stations used ONLY for GA training (never 3, 7, 10 = scenarios D-F)
mp.dayDur = 2.0;
mp.Ncell = 30; mp.Cb = 25; mp.Rb = 0.15; mp.SOC0 = 0.8; mp.SOCmin = 0.2;
mp.socAccel = 600;
mp.etaConv = 0.96;
mp.mission = [3 52 5];
mp.startH = 7; mp.maxH = 72;
mp.tempLevels = [15 25 35 45 55 65];   % cell temperature levels for plot 9 (degC)
mp.tempDur = 2.4;                      % simulation time per temperature level (s)
mp.socChg0 = 0.30; mp.pChg = 50; mp.chgDur = 3.0;   % battery charging test
uav = struct('mass', 8, 'nRot', 4, 'dia', 0.4, 'rho', 1.225, 'FoM', 0.7, 'etaMotor', 0.85, 'etaESC', 0.95, 'tauTO', 1.4, 'tauHov', 1.0, 'tauCruise', 0.7);

%% STEP 1 - DATASET LOADING
if ~isfolder(dataFolder), error('Dataset folder does not exist: %s', dataFolder); end
csvFiles = dir(fullfile(dataFolder, 'station*.csv'));
if isempty(csvFiles), error('No station*.csv files found in: %s', dataFolder); end
[~, sortIndex] = sort(lower({csvFiles.name})); csvFiles = csvFiles(sortIndex);
fprintf('Number of CSV files found: %d\n', numel(csvFiles));
Stations = struct('Name', {}, 'Data', {});
for i = 1:numel(csvFiles)
    [~, stationName, ~] = fileparts(csvFiles(i).name);
    try, T = readStation(fullfile(dataFolder, csvFiles(i).name)); fprintf('Loaded: %s | Rows: %d | Columns: %d\n', csvFiles(i).name, height(T), width(T));
    catch ME, T = table(); fprintf('FAILED: %s | Error: %s\n', csvFiles(i).name, ME.message); end
    Stations(i).Name = stationName; Stations(i).Data = T;
end
if max(scenStations) > numel(Stations), error('scenStations exceeds the number of stations loaded (%d).', numel(Stations)); end

%% STEP 2 - DATA CLEANING
for i = 1:numel(Stations)
    T = Stations(i).Data;
    if isempty(T), continue; end
    n0 = height(T); T = unique(T, 'rows', 'stable');
    vn = string(T.Properties.VariableNames); dIdx = find(contains(lower(vn), ["date","time","时间","日期"]), 1, 'first');
    if ~isempty(dIdx)
        d = T.(vn(dIdx)); [tConv, okT] = parseTime(d);
        if okT, T.(vn(dIdx)) = tConv;
        else, fprintf('%s date conversion failed - sample index will be used. First raw values:\n', Stations(i).Name); disp(d(1:min(3,numel(d)))); end
    end
    Stations(i).Data = T;
    fprintf('%s cleaned | Original: %d | Final: %d\n', Stations(i).Name, n0, height(T));
end
fprintf('STEP 1 + STEP 2 COMPLETED\n');

%% OUTPUT FOLDERS
set(groot, 'DefaultAxesFontName', fontName, 'DefaultAxesFontSize', fontSize, 'DefaultAxesFontWeight', fontWt, 'DefaultTextFontName', fontName, 'DefaultTextFontSize', fontSize, 'DefaultTextFontWeight', fontWt);
outPlots = fullfile(fileparts(dataFolder), 'Final_Plots');  if ~isfolder(outPlots), mkdir(outPlots); end
outTab   = fullfile(fileparts(dataFolder), 'Final_Tables'); if ~isfolder(outTab),   mkdir(outTab);   end
xlFinal  = fullfile(outTab, 'Final_Result_Tables.xlsx');    if isfile(xlFinal), delete(xlFinal); end
nSt = numel(Stations); stNames = cellstr(string({Stations.Name})); letters = 'ABCDEF';

%% STEP 3 - ENVIRONMENTAL DATA EXTRACTION (column detection, no extra plots)
short = {'Irradiance','Temperature','PVPower','WindSpeed','Humidity'};
ylab  = {'Measured Irradiance (W/m^2)','Measured Temperature (°C)','PV Power (MW)','Wind Speed (m/s)','Humidity (%)'};
pats  = {["实测总辐射","measured total","total solar","total irradiance","totalirrad","global horizontal","irradiance","irrad","辐射","辐照","radiation","ghi"], ...
         ["实测温度","air temperature","temperature","temp","温度"], ...
         ["实际发电功率","power","功率","pv"], ...
         ["风速","wind speed","windspeed","wind_speed","wind"], ...
         ["相对湿度","relative humidity","humidity","humid","湿度"]};
colF = [0.64 0.08 0.18];
scCol = {[0.85 0.33 0.10],[0.00 0.45 0.74],[0.47 0.67 0.19],[0.49 0.18 0.56],[0.93 0.69 0.13],[0.64 0.08 0.18]};

Store = repmat(struct('t', [], 'Y', {cell(1,5)}, 'dayMask', [], 'hasTime', false, 'dayId', []), 1, nSt);
ColMap = table();
for i = 1:nSt
    T = Stations(i).Data; sName = Stations(i).Name;
    fprintf('\n============================================================\nDATA EXTRACTION - %s\n============================================================\n', sName);
    if isempty(T), fprintf('%s skipped: no data.\n', sName); continue; end
    vn = string(T.Properties.VariableNames); ln = lower(vn);
    dIdx = find(contains(ln, ["date","时间","日期"]) | (contains(ln, "time") & ~contains(ln, "temp")), 1, 'first');
    if ~isempty(dIdx) && isdatetime(T.(vn(dIdx))), tAll = T.(vn(dIdx)); else, tAll = (1:height(T))'; end
    [tAll, so] = sort(tAll); T = T(so, :);
    Num = cell(1, numel(vn)); for q = 1:numel(vn), Num{q} = toNum(T.(vn(q))); end
    [colIdx, how] = detectCols(Num, vn, ln, pats);
    colName = strings(1,5); for k = 1:5, if colIdx(k) > 0, colName(k) = vn(colIdx(k)); end, end
    for k = 1:5, fprintf('%-12s column : %-40s (%s)\n', short{k}, colName(k), how(k)); end
    if any(colIdx == 0) || any(how == "data")
        fprintf('Column names of %s:\n', sName);
        for q = 1:numel(vn), fprintf('  %d. %s\n', q, vn(q)); end
    end
    ColMap = [ColMap; table(string(sName), colName(1), colName(2), colName(3), colName(4), colName(5), 'VariableNames', {'Station','Irradiance','Temperature','PVPower','WindSpeed','Humidity'})]; %#ok<AGROW>
    Y = cell(1,5); for k = 1:5, if colIdx(k) > 0, Y{k} = Num{colIdx(k)}; end, end
    if ~isempty(Y{1}), dayMask = Y{1} > daylightThreshold; else, dayMask = true(height(T),1); end
    dayMask(isnan(dayMask)) = false;
    Store(i).t = tAll; Store(i).Y = Y; Store(i).dayMask = dayMask; Store(i).hasTime = isdatetime(tAll);
    if isdatetime(tAll), [~, ~, dId] = unique(dateshift(tAll, 'start', 'day')); Store(i).dayId = dId(:); else, Store(i).dayId = ceil((1:height(T))'/96); end
end
fprintf('\nCOLUMNS USED PER STATION\n'); disp(ColMap);

%% PLOT 1 - SOLAR IRRADIANCE DISTRIBUTION ANALYSIS
d = displayStation; stCol = jet(nSt);
if ~isempty(Store(d).Y{1})
    y = Store(d).Y{1}(Store(d).dayMask & ~isnan(Store(d).Y{1}));
    [cnt, edg] = histcounts(y, 60, 'Normalization', 'pdf'); ctr = (edg(1:end-1) + edg(2:end)) / 2;
    fig = newFig(); ax = axes(fig);
    plot(ax, ctr, cnt, '-o', 'Color', colF, 'LineWidth', 2.5, 'MarkerSize', 6, 'MarkerFaceColor', colF);
    xlabel(ax, ylab{1}, 'Interpreter', 'tex'); ylabel(ax, 'Probability Density'); title(ax, sprintf('%s - Solar Irradiance Distribution Analysis', Stations(d).Name), 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
    savePlot(fig, outPlots, '01a_Solar_Irradiance_Distribution_Display_Station');
end
fig = newFig(); ax = axes(fig); hold(ax, 'on'); lgS = {};
for i = 1:nSt
    if isempty(Store(i).Y{1}), continue; end
    y = Store(i).Y{1}(Store(i).dayMask & ~isnan(Store(i).Y{1})); if numel(y) < 2, continue; end
    [cnt, edg] = histcounts(y, 60, 'Normalization', 'pdf'); ctr = (edg(1:end-1) + edg(2:end)) / 2;
    plot(ax, ctr, cnt, 'Color', stCol(i,:), 'LineWidth', 2); lgS{end+1} = stNames{i}; %#ok<SAGROW>
end
xlabel(ax, ylab{1}, 'Interpreter', 'tex'); ylabel(ax, 'Probability Density'); title(ax, 'Solar Irradiance Distribution Analysis - All Stations', 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
if ~isempty(lgS), lgd(ax, lgS, 'eastoutside', fontName, fontSize - 6, fontWt); end
savePlot(fig, outPlots, '01b_Solar_Irradiance_Distribution_All_Stations');

%% PLOT 2 - ENVIRONMENTAL PARAMETER CORRELATION HEATMAP
Yd = Store(d).Y; dm = Store(d).dayMask;
avail = find(~cellfun(@isempty, Yd)); avail = avail(cellfun(@(z) any(~isnan(z)), Yd(avail)));
if numel(avail) >= 2
    D = cell2mat(cellfun(@(z) z(:), Yd(avail), 'UniformOutput', false)); D = D(dm, :);
    C = corr(D, 'Rows', 'pairwise');
    fig = newFig(); ax = axes(fig); imagesc(ax, C); colormap(ax, jet(256)); caxis(ax, [-1 1]);
    cb = colorbar(ax); cb.FontName = fontName; cb.FontSize = fontSize; cb.FontWeight = fontWt;
    set(ax, 'XTick', 1:numel(avail), 'XTickLabel', short(avail), 'YTick', 1:numel(avail), 'YTickLabel', short(avail), 'XTickLabelRotation', 0); axis(ax, 'square');
    for a = 1:numel(avail), for b = 1:numel(avail)
        if abs(C(a,b)) > 0.5, tc = 'w'; else, tc = 'k'; end
        text(ax, b, a, sprintf('%.2f', C(a,b)), 'HorizontalAlignment', 'center', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', fontWt, 'Color', tc);
    end, end
    xlabel(ax, 'Environmental Parameter'); ylabel(ax, 'Environmental Parameter');
    title(ax, sprintf('%s - Environmental Parameter Correlation Heatmap', Stations(d).Name), 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
    savePlot(fig, outPlots, '02_Environmental_Parameter_Correlation_Heatmap');
end

%% STEP 4 - REPRESENTATIVE DAYS AND SCENARIOS A-F
fprintf('\n============================================================\nSTEP 4 - REPRESENTATIVE DAYS AND SCENARIOS\n============================================================\n');
RepBest = NaN(1, nSt); RepTbl = table();
for i = 1:nSt
    Sx = Store(i); lblDay = "N/A"; ns = 0; f = NaN(1,7);
    if ~isempty(Sx.Y{1})
        [best, F, ~] = pickRepDay(Sx.Y, Sx.dayId, daylightThreshold); RepBest(i) = best;
        if ~isnan(best)
            f = F(best,:); ns = sum(Sx.dayId == best);
            if Sx.hasTime, lblDay = string(dateshift(Sx.t(find(Sx.dayId == best, 1)), 'start', 'day'), 'yyyy-MM-dd'); else, lblDay = "Day " + string(best); end
        end
    end
    RepTbl = [RepTbl; table(string(Stations(i).Name), lblDay, ns, f(1), f(2), f(4), f(5), f(6), f(7), 'VariableNames', {'Station','RepresentativeDay','Samples','MeanDaylightIrr','MaxIrr','MeanTemp','MeanPVPower','MeanWind','MeanHumidity'})]; %#ok<AGROW>
end
fprintf('\nREPRESENTATIVE DAY PER STATION\n'); disp(RepTbl);

Sc = cell(1,6);
cand = [scenStations, 1:nSt]; ref = 0;
for c = cand
    if ~isempty(Store(c).Y{1}) && ~isempty(Store(c).Y{3}), ref = c; break; end
end
if ref > 0
    Yr = Store(ref).Y; dm = Store(ref).dayMask; Tr = Yr{2}; if isempty(Tr), Tr = 25*ones(size(Yr{1})); end
    ok = dm & ~isnan(Yr{1}) & ~isnan(Yr{3}) & ~isnan(Tr);
    coef = [Yr{1}(ok), Yr{1}(ok).*(Tr(ok)-25)] \ Yr{3}(ok);
    pvModel = @(G, Tc) max(coef(1)*G + coef(2)*G.*(Tc-25), 0);
    if ~isempty(Yr{4}), medW = median(Yr{4}(dm), 'omitnan'); else, medW = 2; end
    if ~isempty(Yr{5}), medH = median(Yr{5}(dm), 'omitnan'); else, medH = 50; end
else
    pvModel = @(G, Tc) NaN(size(G)); medW = 2; medH = 50;
end
dt = 0.1;
tA = (0:dt:60)'; GA = interp1(linspace(0,60,7), [200 400 700 900 1000 700 300], tA, 'pchip'); TA = 25*ones(size(tA));
Sc{1} = makeSc('Scenario A - Normal day', 'Controlled', tA, GA, TA, pvModel(GA,TA), medW*ones(size(tA)), medH*ones(size(tA)), 1, 'Time (s)');
tB = (0:dt:16)'; GB = interp1([0 4 6 10 12 16], [900 850 300 200 800 950], tB, 'linear'); TB = 25*ones(size(tB));
Sc{2} = makeSc('Scenario B - Rapid cloud transition', 'Controlled', tB, GB, TB, pvModel(GB,TB), medW*ones(size(tB)), medH*ones(size(tB)), 1, 'Time (s)');
tC = (0:dt:30)'; TC = interp1([0 10 20 30], [25 30 35 40], tC, 'linear'); GC = 1000*ones(size(tC));
Sc{3} = makeSc('Scenario C - Temperature variation', 'Controlled', tC, GC, TC, pvModel(GC,TC), medW*ones(size(tC)), medH*ones(size(tC)), 1, 'Time (s)');
for n = 1:3
    s = scenStations(n); Sx = Store(s); sc = 3 + n; best = RepBest(s);
    if isnan(best), fprintf('Scenario %s skipped: no representative day found for %s.\n', letters(sc), Stations(s).Name); continue; end
    m = Sx.dayId == best; tt = Sx.t(m);
    if Sx.hasTime, tsec = seconds(tt - tt(1)); else, tsec = (0:sum(m)-1)' * 900; end
    vals = NaN(sum(m), 5);
    for k = 1:5, if ~isempty(Sx.Y{k}), vals(:,k) = fillmissing(Sx.Y{k}(m), 'linear', 'EndValues', 'nearest'); end, end
    vals(:,[1 3]) = max(vals(:,[1 3]), 0);
    nm = sprintf('Scenario %s - Real profile, Station %d (%s)', letters(sc), s, Stations(s).Name);
    Sc{sc} = makeSc(nm, 'Measured', tsec, vals(:,1), vals(:,2), vals(:,3), vals(:,4), vals(:,5), 1/3600, 'Time (h)');
end
fprintf('STEP 4 COMPLETED.\n');

%% STEP 5 - PV ELECTRICAL MODELLING (SINGLE-DIODE MODEL)
fprintf('\n============================================================\nSTEP 5 - PV ELECTRICAL MODELLING (SINGLE-DIODE)\n============================================================\n');
PV = struct('Name', 'Generic 72-cell poly-Si module 330 W', ...
            'Pmax', 330, 'Vmp', 37.2, 'Imp', 8.87, 'Voc', 45.9, 'Isc', 9.45, 'Ns', 72, ...
            'n', 1.3, 'Ki', 0.0005*9.45, 'Kv', -0.0029*45.9, 'Eg', 1.121, 'U0', 25, 'U1', 6.84, 'NserString', 24);
PV = fitSD(PV);
[PmSTC, VmSTC, ImSTC] = mppSD(1000, 25, PV);
fprintf('Module : %s\n', PV.Name);
fprintf('Single-diode parameters at STC : Iph = %.4f A | I0 = %.4e A | Rs = %.4f ohm | Rsh = %.1f ohm | n = %.2f | Ns = %d cells\n', PV.IphRef, PV.I0Ref, PV.Rs, PV.RshRef, PV.n, PV.Ns);
fprintf('Model check at STC  : Pmax = %.2f W (datasheet %.1f W) | Vmp = %.2f V (%.1f V) | Imp = %.3f A (%.2f A)\n', PmSTC, PV.Pmax, VmSTC, PV.Vmp, ImSTC, PV.Imp);
fprintf('Calibration scheme  : fixed eta = %.2f | fitted: gT in [%.2f %.2f] (prior %.2f), irradiance exponent beta in [%.2f %.2f] (prior %.2f), effective array size (closed form)\n', etaSys, cc.gT(1), cc.gT(2), cc.gT0, cc.beta(1), cc.beta(2), cc.b0);
fprintf('                      parameters that end at a bound are FIXED to their physical value and the rest is refitted (see CalStatus)\n');
fprintf('                      interleaved day split (test days: mod(day,10) in %s)\n', mat2str(testDayMod));
ParamTbl = table(["Pmax_W";"Vmp_V";"Imp_A";"Voc_V";"Isc_A";"Cells_Ns";"Ideality_n";"Ki_A_per_C";"Kv_V_per_C";"Iph_ref_A";"I0_ref_A";"Rs_ohm";"Rsh_ref_ohm";"ModulesPerString"], ...
                 [PV.Pmax; PV.Vmp; PV.Imp; PV.Voc; PV.Isc; PV.Ns; PV.n; PV.Ki; PV.Kv; PV.IphRef; PV.I0Ref; PV.Rs; PV.RshRef; PV.NserString], 'VariableNames', {'Parameter','Value'});
STCTbl = table(["Pmax_W";"Vmp_V";"Imp_A"], [PV.Pmax; PV.Vmp; PV.Imp], [PmSTC; VmSTC; ImSTC], 'VariableNames', {'Parameter','Datasheet','SingleDiodeModel'});

capIn = capacityMW(:)'; if numel(capIn) < nSt, capIn(end+1:nSt) = NaN; end
R5 = cell(1, nSt); CalTbl = table();
for i = 1:nSt
    Sx = Store(i); Y = Sx.Y; sName = Stations(i).Name;
    fprintf('\n------------------------------------------------------------\nPV MODEL - %s\n------------------------------------------------------------\n', sName);
    if isempty(Y{1}) || isempty(Y{3}), fprintf('%s skipped: irradiance or PV power column missing.\n', sName); continue; end
    n = numel(Y{1});
    G = Y{1}; G(isnan(G)) = 0; G = max(G, 0);
    Pm = Y{3}; Pm(Pm < 0) = 0;
    Ta = 20*ones(n,1); if ~isempty(Y{2}) && any(~isnan(Y{2})), Ta = fillmissing(Y{2}, 'linear', 'EndValues', 'nearest'); end
    ws = ones(n,1);    if ~isempty(Y{4}) && any(~isnan(Y{4})), ws = max(fillmissing(Y{4}, 'linear', 'EndValues', 'nearest'), 0); end
    cap = capIn(i); if isnan(cap), cap = ceil(2*prctile(Pm(~isnan(Pm)), 99.9))/2; end
    Np0 = max(round(cap*1e6*dcAcRatio / (PV.NserString*PV.Pmax)), 1);            % nominal (uncalibrated) array
    dm = Sx.dayMask; testDay = ismember(mod(Sx.dayId, 10), testDayMod); trainDay = ~testDay;
    valid = dm & ~isnan(Pm);
    calIdx = find(valid & trainDay & G > 100 & Pm > 0.01*cap & Pm < 0.97*cap);
    if numel(calIdx) < 20, fprintf('%s skipped: not enough calibration samples.\n', sName); continue; end
    if numel(calIdx) > 1500, calIdx = calIdx(round(linspace(1, numel(calIdx), 1500))); end
    idx = calIdx; nCal0 = numel(calIdx);
    for pass = 1:2
        [pf, kfit, calStatus] = calibrate(G(idx), Ta(idx), ws(idx), Pm(idx), PV, cc);
        if pass == 1                                                   % robust trimming of curtailment / outlier points
            Pu = unitPower(pf, G(idx), Ta(idx), ws(idx), PV);
            res = Pm(idx) - kfit*Pu; rc = res - median(res); sg = 1.4826*median(abs(rc));
            keep = abs(rc) < 3*max(sg, eps);
            if sum(keep) >= 20, idx = idx(keep); end
        end
    end
    gT = pf(1); beta = pf(2);
    effStr = kfit/etaSys; Np = max(round(effStr), 1); eta = etaSys;
    scale = PV.NserString * Np; kWp = scale*PV.Pmax/1e3;
    fprintf('Installed AC capacity = %.2f MW | nominal DC array = %.1f kWp (%d strings)\n', cap, Np0*PV.NserString*PV.Pmax/1e3, Np0);
    fprintf('Calibration : %d of %d samples used after trimming | effective DC array = %.1f kWp (%d strings, continuous %.1f) | eta fixed = %.2f\n', numel(idx), nCal0, kWp, Np, effStr, eta);
    fprintf('Calibrated  : gT = %.3f | beta = %.3f | status = %s\n', gT, beta, calStatus);
    if calStatus ~= "Fitted", fprintf('NOTE: %s -> the data could not identify that parameter inside its physical range, so it was set to its physical value (no bound-limited result is used).\n', calStatus); end
    [~, ~, Pd0] = simArray(G, Ta, ws, PV, 1, Np0);
    Psim0 = min(Pd0, cap);
    Ge = effG(G, beta);
    Tc = cellTemp(Ge, Ta, ws, PV, gT);
    [Vpv, Ipv, Pdc] = simArray(Ge, Ta, ws, PV, gT, Np);
    Psim = min(eta*Pdc, cap);
    trM = valid & trainDay; teM = valid & testDay; alM = valid;
    m0 = metrics(Pm(alM), Psim0(alM), cap);
    mTr = metrics(Pm(trM), Psim(trM), cap); mTe = metrics(Pm(teM), Psim(teM), cap); mAl = metrics(Pm(alM), Psim(alM), cap);
    fprintf('Before calibration (all)  : RMSE = %.4f MW | R2 = %.4f\n', m0.RMSE, m0.R2);
    fprintf('After  calibration TRAIN  : RMSE = %.4f MW | nRMSE = %.2f %% | MBE = %.4f MW | R2 = %.4f\n', mTr.RMSE, mTr.nRMSE, mTr.MBE, mTr.R2);
    fprintf('After  calibration TEST   : RMSE = %.4f MW | nRMSE = %.2f %% | MAE = %.4f MW | MBE = %.4f MW | R2 = %.4f | CC = %.4f\n', mTe.RMSE, mTe.nRMSE, mTe.MAE, mTe.MBE, mTe.R2, mTe.CC);
    CalTbl = [CalTbl; table(string(sName), cap, kWp, Np, Np0, effStr, eta, gT, beta, calStatus, m0.RMSE, m0.R2, mTr.RMSE, mTr.R2, mTe.RMSE, mTe.nRMSE, mTe.MAE, mTe.MBE, mTe.R2, mTe.CC, mAl.R2, ...
        'VariableNames', {'Station','Capacity_MW','DC_kWp','ParallelStrings','ParallelStrings_Initial','EffStrings_Continuous','Eta_Fixed','ThermalFactor','IrradianceExponent','CalStatus','RMSE_Before','R2_Before','RMSE_Train','R2_Train','RMSE_Test','nRMSE_Test_pct','MAE_Test','MBE_Test','R2_Test','CC_Test','R2_All'})]; %#ok<AGROW>
    R5{i} = struct('Vpv', Vpv, 'Ipv', Ipv, 'Pdc', Pdc, 'Psim', Psim, 'Pmeas', Pm, 'Tcell', Tc, 'G', G, 'Ge', Ge, 'Ta', Ta, 'ws', ws, ...
                   'trM', trM, 'teM', teM, 'valid', valid, 'par', struct('eta', eta, 'gT', gT, 'beta', beta, 'cap', cap, 'Np', Np));
end
fprintf('\nCALIBRATION / VALIDATION SUMMARY - ALL STATIONS\n'); disp(CalTbl);
if height(CalTbl) > 0
    nFix = sum(CalTbl.CalStatus ~= "Fitted");
    fprintf('Stations fully fitted : %d of %d | stations with a fixed (physical) parameter : %d\n', height(CalTbl) - nFix, height(CalTbl), nFix);
    fprintf('gT range = [%.3f %.3f] | beta range = [%.3f %.3f] | mean test R2 = %.4f | test R2 range = [%.4f %.4f] | test nRMSE range = [%.2f %.2f] %%\n', min(CalTbl.ThermalFactor), max(CalTbl.ThermalFactor), min(CalTbl.IrradianceExponent), max(CalTbl.IrradianceExponent), mean(CalTbl.R2_Test, 'omitnan'), min(CalTbl.R2_Test), max(CalTbl.R2_Test), min(CalTbl.nRMSE_Test_pct), max(CalTbl.nRMSE_Test_pct));
    fprintf('NOTE: the PV model gives good-to-moderate, station-dependent agreement. Do not claim it accurately predicts all stations.\n');
    fprintf('NOTE: gT, beta and the effective array size are EFFECTIVE (lumped) parameters of a regional-data model, not measured module properties or installed capacity.\n');
end

%% PLOT 3 - PV POWER PREDICTION VS ACTUAL POWER
%% PLOT 4 - PV MODEL VALIDATION ERROR
d = displayStation;
if ~isempty(R5{d})
    Sx = Store(d); r = R5{d}; sName = Stations(d).Name; cap = r.par.cap;
    best = RepBest(d); if isnan(best), [~, ix] = max(r.G); best = Sx.dayId(ix); end
    m = Sx.dayId == best; if Sx.hasTime, tx = Sx.t(m); tLab5 = 'Time'; else, tx = find(m); tLab5 = 'Sample Index'; end
    te = r.teM;
    fig = newFig(); ax = axes(fig); hold(ax, 'on');
    plot(ax, tx, r.Pmeas(m), '-o', 'Color', [0.64 0.08 0.18], 'LineWidth', 2.5, 'MarkerSize', 5); plot(ax, tx, r.Psim(m), '-s', 'Color', [0.00 0.45 0.74], 'LineWidth', 2.5, 'MarkerSize', 5);
    xlabel(ax, tLab5); ylabel(ax, 'PV Power (MW)'); title(ax, sprintf('%s - PV Power Prediction vs Actual Power (Representative Day)', sName), 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
    lgd(ax, {'Actual', 'Predicted'}, 'best', fontName, fontSize, fontWt);
    savePlot(fig, outPlots, '03a_PV_Power_Prediction_vs_Actual_RepDay');
    fig = newFig(); ax = axes(fig); hold(ax, 'on');
    scatter(ax, r.Pmeas(te), r.Psim(te), 18, [0.00 0.45 0.74], 'filled', 'MarkerFaceAlpha', 0.35); plot(ax, [0 cap], [0 cap], 'k--', 'LineWidth', 2.5);
    mt = metrics(r.Pmeas(te), r.Psim(te), cap);
    xlabel(ax, 'Actual PV Power (MW)'); ylabel(ax, 'Predicted PV Power (MW)'); title(ax, sprintf('%s - Predicted vs Actual, Test Set (R^2 = %.3f, nRMSE = %.2f %%)', sName, mt.R2, mt.nRMSE), 'Interpreter', 'tex'); styleAx(ax, fontName, fontSize, fontWt);
    savePlot(fig, outPlots, '03b_PV_Power_Prediction_vs_Actual_Scatter_Test');

    fig = barSingle(CalTbl.nRMSE_Test_pct, cellstr(CalTbl.Station), 'PV Model Validation Error - nRMSE per Station (Test Set)', 'Station', 'nRMSE (% of capacity)', [], fontName, fontSize, fontWt, '%.2f');
    savePlot(fig, outPlots, '04a_PV_Model_Validation_Error_nRMSE');
    fig = newFig(); ax = axes(fig); eP = (r.Psim(te) - r.Pmeas(te)) / cap * 100; eP = eP(~isnan(eP));
    histogram(ax, eP, 60, 'FaceColor', [0.00 0.45 0.74], 'Normalization', 'probability');
    xlabel(ax, 'Prediction Error (% of capacity)'); ylabel(ax, 'Probability'); title(ax, sprintf('%s - PV Model Validation Error Distribution (Test Set)', sName), 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
    savePlot(fig, outPlots, '04b_PV_Model_Validation_Error_Distribution');
end

%% STEP 5.4 - CALIBRATED PV MODEL APPLIED TO SCENARIOS A-F
%   Controlled scenarios A-C are defined directly in the array plane (beta not applied); real profiles D-F use the station beta.
Sc5 = cell(1,6); Sc5Sum = table();
for n = 1:6
    if isempty(Sc{n}), continue; end
    if n <= 3, st = scenStations(1); else, st = scenStations(n-3); end
    if isempty(R5{st}), fprintf('Scenario %s skipped: no calibrated model for %s.\n', letters(n), Stations(st).Name); continue; end
    pr = R5{st}.par; s = Sc{n};
    G = s.Irradiance; Ta = s.Temperature; ws = s.WindSpeed; ws(isnan(ws)) = 1;
    if n <= 3, Ge = G; else, Ge = effG(G, pr.beta); end
    Tc = cellTemp(Ge, Ta, ws, PV, pr.gT);
    [Vpv, Ipv, Pdc] = simArray(Ge, Ta, ws, PV, pr.gT, pr.Np);
    Psim = min(pr.eta*Pdc, pr.cap);
    Sc5{n} = struct('Name', s.Name, 'Time_s', s.Time_s, 'tPlot', s.tPlot, 'tLabel', s.tLabel, 'G', G, 'Ta', Ta, 'Tcell', Tc, 'Vpv', Vpv, 'Ipv', Ipv, 'Pdc', Pdc, 'Psim', Psim, 'Ref', s.PVPower);
    mm = metrics(s.PVPower, Psim, pr.cap); isMeas = strcmp(s.Source, 'Measured');
    fprintf('%s | PV model of %s | Vpv max = %.1f V | Ppv max = %.3f MW\n', s.Name, Stations(st).Name, max(Vpv), max(Psim));
    Sc5Sum = [Sc5Sum; table(string(s.Name), string(Stations(st).Name), max(Vpv), mean(Vpv), max(Ipv), mean(Ipv), max(Psim), mean(Psim), max(Tc), ifelse(isMeas, mm.RMSE, NaN), ifelse(isMeas, mm.R2, NaN), ...
        'VariableNames', {'Scenario','ModelStation','MaxVpv_V','MeanVpv_V','MaxIpv_A','MeanIpv_A','MaxPsim_MW','MeanPsim_MW','MaxCellTemp_C','RMSE_vs_Measured_MW','R2_vs_Measured'})]; %#ok<AGROW>
end

%% PLOT 5 - MPP VOLTAGE VARIATION UNDER DYNAMIC CONDITIONS
for n = 1:6
    if isempty(Sc5{n}), continue; end
    fig = newFig(); ax = axes(fig);
    plot(ax, Sc5{n}.tPlot, Sc5{n}.Vpv, 'Color', scCol{n}, 'LineWidth', 2.5);
    xlabel(ax, Sc5{n}.tLabel); ylabel(ax, 'V_{mpp} (V)'); title(ax, sprintf('MPP Voltage Variation Under Dynamic Conditions - Scenario %s', letters(n)), 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
    savePlot(fig, outPlots, sprintf('05%s_MPP_Voltage_Variation_Scenario_%s', lower(letters(n)), letters(n)));
end

%% STEP 8 - MPP EXTRACTION (MPP REFERENCE DATASET)
fprintf('\n============================================================\nSTEP 8 - MPP EXTRACTION\n============================================================\n');
MPPRef = table(); MPPSum = table();
for i = 1:nSt
    if isempty(R5{i}), fprintf('%s skipped: no PV model result.\n', Stations(i).Name); continue; end
    r = R5{i}; Sx = Store(i); sName = Stations(i).Name; Np = r.par.Np;
    k = find(Sx.dayMask & r.Ge > 1);
    if isempty(k), continue; end
    nk = numel(k);
    [Pmod, Vmod, Imod] = mppSD(r.Ge(k), r.Tcell(k), PV);
    Vmpp = Vmod * PV.NserString; Impp = Imod * Np; Pmpp = Vmpp .* Impp;
    if Sx.hasTime, ts = Sx.t(k); ts = ts(:); else, ts = NaT(nk,1); end
    chk = max(abs(Pmpp/1e6 - r.Pdc(k)));
    fprintf('%s : %d daylight conditions | Pmpp max = %.3f MW | consistency with Step 5 (max abs diff) = %.2e MW\n', sName, nk, max(Pmpp)/1e6, chk);
    Ti = table(repmat("S" + string(i), nk, 1), repmat(string(sName), nk, 1), ts, k(:), r.G(k), r.Ta(k), r.Tcell(k), Vmpp, Impp, Pmpp, Pmpp/1e6, Vmod, Imod, Pmod, r.Pmeas(k), ...
        'VariableNames', {'Station','StationName','Timestamp','SampleIndex','Irradiance','Temperature','CellTemperature','Vmpp','Impp','Pmpp','Pmpp_MW','Vmpp_Module','Impp_Module','Pmpp_Module_W','MeasuredPower_MW'});
    MPPRef = [MPPRef; Ti]; %#ok<AGROW>
    MPPSum = [MPPSum; table("S" + string(i), string(sName), nk, mean(r.G(k)), max(r.G(k)), 100*mean(r.G(k) > 1000), mean(r.Ta(k)), mean(Vmpp), mean(Impp), mean(Pmpp)/1e6, max(Pmpp)/1e6, Np*PV.NserString*PV.Pmax/1e3, Np, chk, ...
        'VariableNames', {'Station','StationName','Samples','MeanIrradiance','MaxIrradiance','Pct_Samples_G_above_1000','MeanTemperature','MeanVmpp_V','MeanImpp_A','MeanPmpp_MW','MaxPmpp_MW','DC_kWp','ParallelStrings','ConsistencyDiff_MW'})]; %#ok<AGROW>
end
fprintf('\nMPP SUMMARY PER STATION\n'); disp(MPPSum);
fprintf('Total rows in MPP Reference Dataset : %d\n', height(MPPRef));
fprintf('NOTE: zero consistency difference proves internal consistency of the model, NOT experimental validation.\n');
MPPSumX = MPPSum;
if height(MPPRef) > 0
    MPPSumX = [MPPSum; table("ALL", "All stations", height(MPPRef), mean(MPPRef.Irradiance), max(MPPRef.Irradiance), 100*mean(MPPRef.Irradiance > 1000), mean(MPPRef.Temperature), mean(MPPRef.Vmpp), mean(MPPRef.Impp), mean(MPPRef.Pmpp_MW), max(MPPRef.Pmpp_MW), sum(MPPSum.DC_kWp), NaN, max(MPPSum.ConsistencyDiff_MW), 'VariableNames', MPPSum.Properties.VariableNames)];
    writetable(MPPRef, fullfile(outTab, 'MPP_Reference_Dataset.csv'));
    save(fullfile(outTab, 'MPP_Reference_Dataset.mat'), 'MPPRef', 'MPPSum', 'PV');
end

%% STEP 9 - GA-SMC MPPT DESIGN
%   e = (dI/dV + I/V)/(I/V)  (normalised INC error) | S = e + lambda*de/dt | D(k+1) = D(k) - K*Ts*sat(S/phi) | GA searches lambda, K, phi
%   J = mean over training profiles [w1*E_tracking + w2*T_convergence + w3*P_ripple] + 0.5*worst-profile tracking error + ripple penalty
fprintf('\n============================================================\nSTEP 9-18 - GA-SMC MPPT FRAMEWORK\n============================================================\n');
fprintf('SCALING CHAIN : China regional PV data (MW-scale stations) -> real irradiance / temperature / wind profiles\n');
fprintf('  -> calibrated single-diode PV model (effective, lumped parameters) -> representative environmental conditions\n');
fprintf('  -> SCALED reference array (%d x %d modules of 330 W = %.0f W STC) -> MPPT controller -> boost converter -> %d V battery bus -> UAV propulsion load.\n', mp.nSer, mp.nPar, mp.nSer*mp.nPar*330, mp.Vbat);
fprintf('  The MW-scale stations are NOT installed on the UAV: only their environmental profiles are reused.\n');
fprintf('Convergence definition used EVERYWHERE : first time P >= 95 %% of MPP and staying there for 20 ms.\n');
algCol = [0.00 0.45 0.74; 0.93 0.69 0.13; 0.64 0.08 0.18; 0.47 0.67 0.19; 0.49 0.18 0.56];
algName = {'P&O', 'INC', 'GA-SMC', 'GA-SMC + GMPP scan', 'SMC only'};
aS = mp.nSer; aP = mp.nPar;
gT9 = 1; if ~isempty(R5{displayStation}), gT9 = R5{displayStation}.par.gT; end
cfg0 = struct('dt', mp.dt, 'N', 0, 'L', mp.L, 'rL', mp.rL, 'Cin', mp.Cin, 'Dmin', mp.Dmin, 'Dmax', mp.Dmax, 'D0', mp.D0, 'Vbat', mp.Vbat, 'useBat', false, 'Pload', [], ...
              'TsSMC', mp.TsSMC, 'TsPO', mp.TsPO, 'stepD', mp.stepD, 'incEps', mp.incEps, 'scanT', 0.3, 'scanHold', 0.5, ...
              'Ncell', mp.Ncell, 'Cb', mp.Cb, 'Rb', mp.Rb, 'SOC0', mp.SOC0, 'socAccel', mp.socAccel);
segGtr = [1000; 400; 800; 200; 600; 1000]; TaSeg = [25; 35; 15; 30; 20; 40]; nSegTr = numel(segGtr);
Ntr = round(mp.gaTsim/mp.dt); ttr = (0:Ntr-1)' * mp.dt; sidx = min(floor(ttr/(mp.gaTsim/nSegTr)) + 1, nSegTr);
Gtr = segGtr(sidx); Tctr = cellTemp(Gtr, TaSeg(sidx), 1, PV, gT9);
Pref_tr = mppSD(Gtr, Tctr, PV) * aS * aP;
plTr = makeUniformPlant(Gtr, Tctr, PV, aS, aP); cfgTr = cfg0; cfgTr.N = Ntr; cfgTr.nSeg = nSegTr;
segG = [1000; 400; 800];
fprintf('\nBase training profile : %.2f s, G = %s W/m^2, reference MPP = %s W\n', mp.gaTsim, mat2str(segGtr'), mat2str(round(unique(Pref_tr, 'stable')')));

% ---- multi-condition GA training set ----
TP = {};
TP{end+1} = struct('pl', plTr, 'Pref', Pref_tr, 'N', Ntr, 'rippleBase', 1);
tBt = (0:0.1:16)'; GBt = interp1([0 4 6 10 12 16], [900 850 300 200 800 950], tBt, 'linear');
[pB_, PB_] = buildScenario((1:numel(tBt))', GBt, 25*ones(size(tBt)), ones(size(tBt)), mp.gaTsim, mp.dt, PV, gT9, aS, aP);
TP{end+1} = struct('pl', pB_, 'Pref', PB_, 'N', pB_.N, 'rippleBase', 1);
tCt = (0:0.1:30)'; TCt = interp1([0 10 20 30], [25 30 35 40], tCt, 'linear');
[pC_, PC_] = buildScenario((1:numel(tCt))', 1000*ones(size(tCt)), TCt, ones(size(tCt)), mp.gaTsim, mp.dt, PV, gT9, aS, aP);
TP{end+1} = struct('pl', pC_, 'Pref', PC_, 'N', pC_.N, 'rippleBase', 1);
for st = mp.trStations
    if isnan(RepBest(st)) || isempty(R5{st}), continue; end
    [~, Gd_, Tad_, wsd_] = repDayProfile(Store(st), RepBest(st));
    ii_ = find(Gd_ > daylightThreshold); if numel(ii_) < 10, continue; end
    ix_ = ii_(1):ii_(end);
    [pR_, PR_] = buildScenario((1:numel(ix_))', effG(Gd_(ix_), R5{st}.par.beta), Tad_(ix_), wsd_(ix_), mp.gaTsim, mp.dt, PV, R5{st}.par.gT, aS, aP);
    TP{end+1} = struct('pl', pR_, 'Pref', PR_, 'N', pR_.N, 'rippleBase', 1); %#ok<SAGROW>
end
for q = 1:numel(TP)                                   % ripple limit = lowest ripple of P&O, INC, default SMC on that profile
    cq = cfg0; cq.N = TP{q}.N; rb = [];
    for a = [1 2 5]
        o = runAlg(TP{q}.pl, a, cq, [], mp); mt = mpptMetrics(o.t, o.V .* o.I, TP{q}.Pref); rb(end+1) = mt.RipplePct; %#ok<AGROW>
    end
    TP{q}.rippleBase = max(min(rb), 1e-6);
end
fprintf('GA training set : %d profiles (training steps, cloud transition, temperature ramp, real days of stations %s)\n', numel(TP), mat2str(mp.trStations));
fitFun = @(x) gaFitness(x, TP, cfg0, mp);

rng(mp.gaSeed);
popN = mp.gaPop; nGen = mp.gaGen; pop = rand(popN, 3);
pop(1,:) = [mp.smcDefault(1)/mp.lamMax, ...
            (log10(mp.smcDefault(2)) - log10(mp.Krange(1)))/(log10(mp.Krange(2)) - log10(mp.Krange(1))), ...
            (log10(mp.smcDefault(3)) - log10(mp.phiRange(1)))/(log10(mp.phiRange(2)) - log10(mp.phiRange(1)))];   % seed = default SMC
pop(1,:) = min(max(pop(1,:), 0), 1);
Jp = zeros(popN, 1); for p = 1:popN, Jp(p) = fitFun(pop(p,:)); end
[Jp, oo] = sort(Jp); pop = pop(oo,:);
bestHist = zeros(nGen+1, 1); meanHist = zeros(nGen+1, 1); bestHist(1) = Jp(1); meanHist(1) = mean(min(Jp, 10));
fprintf('GA generation  0 / %d : best J = %.5f | mean J = %.5f\n', nGen, bestHist(1), meanHist(1));
for g = 1:nGen
    nEl = 2; newPop = pop(1:nEl, :); newJ = Jp(1:nEl);
    while size(newPop, 1) < popN
        ia = randi(popN, 3, 1); [~, ba] = min(Jp(ia)); pa = pop(ia(ba), :);
        ib = randi(popN, 3, 1); [~, bb] = min(Jp(ib)); pb = pop(ib(bb), :);
        c = pa + (rand(1,3)*1.6 - 0.3) .* (pb - pa);
        c = c + (rand(1,3) < 0.25) .* 0.1 .* randn(1,3);
        c = min(max(c, 0), 1);
        newPop(end+1, :) = c; newJ(end+1, 1) = fitFun(c); %#ok<AGROW>
    end
    [Jp, oo] = sort(newJ); pop = newPop(oo, :);
    bestHist(g+1) = Jp(1); meanHist(g+1) = mean(min(Jp, 10));
    fprintf('GA generation %2d / %d : best J = %.5f | mean J = %.5f\n', g, nGen, bestHist(g+1), meanHist(g+1));
end
smcOpt = decodeGene(pop(1,:), mp); [Jbest, compBest] = gaFitness(pop(1,:), TP, cfg0, mp);
fprintf('\nOPTIMISED SMC PARAMETERS : lambda = %.5f s | K = %.4f 1/s | boundary layer phi = %.5f\n', smcOpt(1), smcOpt(2), smcOpt(3));
fprintf('Fitness J = %.5f  (mean E_tracking = %.4f | mean T_convergence = %.4f | mean P_ripple = %.4f)\n', Jbest, compBest(1), compBest(2), compBest(3));
fprintf('Bounds check : lambda [0 %.3f] | K [%.1f %.1f] | phi [%.4f %.3f]\n', mp.lamMax, mp.Krange(1), mp.Krange(2), mp.phiRange(1), mp.phiRange(2));
fprintf('At-bound flags (gene near 0 or 1) [lambda K phi]: %s\n', mat2str(pop(1,:) < 0.01 | pop(1,:) > 0.99));
aList = [1 2 5 3]; TrTbl = table();
for a = aList
    o = runAlg(plTr, a, cfgTr, smcOpt, mp); mt = mpptMetrics(o.t, o.V .* o.I, Pref_tr);
    TrTbl = [TrTbl; table(string(algName{a}), mt.Eff, mt.TrkErrPct, mt.Conv_ms, mt.RipplePct, 'VariableNames', {'Method','Efficiency_pct','TrackingError_pct','Convergence_ms','Ripple_pct'})]; %#ok<AGROW>
end
fprintf('\nBASE TRAINING-PROFILE RESPONSE\n'); disp(TrTbl);
cSMC = TrTbl.Convergence_ms(TrTbl.Method == "SMC only"); cGA = TrTbl.Convergence_ms(TrTbl.Method == "GA-SMC");
fprintf('GA-SMC vs SMC-only convergence on the base training profile : %+.1f ms (differences below the 1 ms controller sampling time are not meaningful).\n', cGA - cSMC);

%% STEP 10 - DC-DC BOOST CONVERTER DESIGN
[~, Vm0, Im0] = mppSD(1000, 25, PV); Vin0 = Vm0*aS; Iin0 = Im0*aP; Dn = 1 - Vin0/mp.Vbat;
dIL = Vin0*Dn/(mp.L*mp.fsw); dVc = dIL/(8*mp.fsw*mp.Cin); LminDes = Vin0*Dn/(0.2*Iin0*mp.fsw);
fprintf('\nSTEP 10 - BOOST CONVERTER DESIGN\n');
fprintf('PV array %d x %d modules : Vmp = %.1f V | Imp = %.2f A | Pmp = %.0f W (STC)\n', aS, aP, Vin0, Iin0, Vin0*Iin0);
fprintf('Bus voltage = %.0f V | nominal duty = %.3f | L = %.2f mH | Cin = %.0f uF | rL = %.3f ohm | fsw = %.0f kHz\n', mp.Vbat, Dn, mp.L*1e3, mp.Cin*1e6, mp.rL, mp.fsw/1e3);
fprintf('Inductor ripple = %.3f A (%.1f %%) | Cin ripple = %.4f V | L for 20 %% ripple = %.2f mH\n', dIL, 100*dIL/Iin0, dVc, LminDes*1e3);

%% STEP 11 - DYNAMIC SCENARIOS A-F : P&O, INC, GA-SMC (and SMC-only for A-C), each compressed to mp.dayDur seconds
fprintf('\n============================================================\nSTEP 11 - DYNAMIC ENVIRONMENTAL SCENARIOS\n============================================================\n');
DynTbl = table(); DynEff = NaN(6,5); DynConv = NaN(6,5); DynRip = NaN(6,5); DynTrk = NaN(6,5); DynCE = NaN(6,5);
DynDuty = cell(6,5); DynT = cell(6,1);
for n = 1:6
    if isempty(Sc{n}), continue; end
    if n >= 4 && ~isempty(R5{scenStations(n-3)}), gTn = R5{scenStations(n-3)}.par.gT; betan = R5{scenStations(n-3)}.par.beta; else, gTn = gT9; betan = 1; end
    Gs = Sc{n}.Irradiance; idx = find(Gs > daylightThreshold); if numel(idx) < 10, fprintf('Scenario %s skipped (too few daylight samples).\n', letters(n)); continue; end
    ix = idx(1):idx(end);
    [plN, PrefN] = buildScenario((1:numel(ix))', effG(Gs(ix), betan), Sc{n}.Temperature(ix), Sc{n}.WindSpeed(ix), mp.dayDur, mp.dt, PV, gTn, aS, aP);
    cfgN = cfg0; cfgN.N = plN.N;
    if n <= 3, aL = [1 2 3 5]; else, aL = [1 2 3]; end
    for a = aL
        o = runAlg(plN, a, cfgN, smcOpt, mp); mt = mpptMetrics(o.t, o.V .* o.I, PrefN); ce = convEff(o);
        DynEff(n,a) = mt.Eff; DynConv(n,a) = mt.Conv_ms; DynRip(n,a) = mt.RipplePct; DynTrk(n,a) = mt.TrkErrPct; DynCE(n,a) = ce;
        DynDuty{n,a} = o.D(1:10:end); DynT{n} = o.t(1:10:end);
        DynTbl = [DynTbl; table(string(Sc{n}.Name), string(Sc{n}.Source), string(algName{a}), mt.Eff, mt.Conv_ms, mt.RipplePct, mt.TrkErrPct, ce, mean(PrefN), ...
            'VariableNames', {'Scenario','Source','Method','Efficiency_pct','Convergence_ms','Ripple_pct','TrackingError_pct','ConverterEfficiency_pct','MeanReferenceMPP_W'})]; %#ok<AGROW>
    end
    fprintf('Scenario %s | Efficiency: P&O = %.2f %% | INC = %.2f %% | GA-SMC = %.2f %%\n', letters(n), DynEff(n,1), DynEff(n,2), DynEff(n,3));
end

%% PLOT 8 - DUTY CYCLE VARIATION ANALYSIS
for n = 1:6
    if isempty(DynT{n}), continue; end
    fig = newFig(); ax = axes(fig); hold(ax, 'on');
    for a = 1:3, plot(ax, DynT{n}, DynDuty{n,a}, 'Color', algCol(a,:), 'LineWidth', 2); end
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Duty Cycle D'); title(ax, sprintf('Duty Cycle Variation Analysis - Scenario %s', letters(n)), 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
    lgd(ax, algName(1:3), 'best', fontName, fontSize - 4, fontWt);
    savePlot(fig, outPlots, sprintf('08%s_Duty_Cycle_Variation_Scenario_%s', lower(letters(n)), letters(n)));
end

%% PLOT 11 - DC-DC CONVERTER EFFICIENCY COMPARISON
vs = find(~isnan(DynCE(:,1)))';
catsS = arrayfun(@(q) ['Scenario ' letters(q)], vs, 'UniformOutput', false);
fig = barGrouped(DynCE(vs,1:3), catsS, algCol(1:3,:), 'DC-DC Converter Efficiency Comparison', 'Scenario', 'DC-DC Converter Efficiency (%)', algName(1:3), fontName, fontSize, fontWt, '%.2f', 0, true);
savePlot(fig, outPlots, '11_DCDC_Converter_Efficiency_Comparison');

%% STEP 12 - BATTERY CHARGING TEST (PV -> MPPT + DC/DC -> battery, standby load only)
socV = linspace(0, 1, 11)'; vocPack = mp.Ncell*(3.3 + 0.9*socV);
EnomWh = mp.Cb * mp.Ncell*(3.3 + 0.9*0.5);
fprintf('\nSTEP 12 - BATTERY : %d series cells | %.0f Ah | OCV %.1f ... %.1f V | pack resistance %.3f ohm | energy about %.2f kWh\n', mp.Ncell, mp.Cb, vocPack(1), vocPack(end), mp.Rb, EnomWh/1e3);
Nc = round(mp.chgDur/mp.dt); tcg = (0:Nc-1)' * mp.dt;
Gc = 850*ones(Nc,1); Tcc = cellTemp(Gc, 30, 1, PV, gT9); plC = makeUniformPlant(Gc, Tcc, PV, aS, aP);
cfgC = cfg0; cfgC.N = Nc; cfgC.useBat = true; cfgC.Pload = mp.pChg*ones(Nc,1); cfgC.SOC0 = mp.socChg0;
oC = cell(1,3); ChgTbl = table();
for a = 1:3
    oC{a} = runAlg(plC, a, cfgC, smcOpt, mp);
    Pch = -oC{a}.Vbat .* oC{a}.Ib; Ech = trapz(tcg, Pch)/3600 * mp.socAccel;
    ChgTbl = [ChgTbl; table(string(algName{a}), 100*oC{a}.SOC(1), 100*oC{a}.SOC(end), 100*(oC{a}.SOC(end) - oC{a}.SOC(1)), Ech, mean(Pch), ...
        'VariableNames', {'Method','SOC_Start_pct','SOC_End_pct','SOC_Gain_pct','ChargeEnergy_Wh','MeanChargingPower_W'})]; %#ok<AGROW>
end
fprintf('\nBATTERY CHARGING TEST (simulation-level result, time-accelerated SOC model x%d)\n', mp.socAccel); disp(ChgTbl);

%% PLOT 12 - BATTERY SOC VARIATION DURING CHARGING
fig = newFig(); ax = axes(fig); hold(ax, 'on'); for a = 1:3, plot(ax, tcg, oC{a}.SOC*100, 'Color', algCol(a,:), 'LineWidth', 2.5); end
xlabel(ax, 'Time (s)'); ylabel(ax, 'Battery SOC (%)'); title(ax, sprintf('Battery SOC Variation During Charging (SOC time scale x%d)', mp.socAccel), 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt); lgd(ax, algName(1:3), 'northwest', fontName, fontSize - 4, fontWt);
savePlot(fig, outPlots, '12a_Battery_SOC_Variation_During_Charging');
fig = newFig(); ax = axes(fig); hold(ax, 'on'); for a = 1:3, plot(ax, tcg, -oC{a}.Vbat .* oC{a}.Ib, 'Color', algCol(a,:), 'LineWidth', 2); end
xlabel(ax, 'Time (s)'); ylabel(ax, 'Battery Charging Power (W)'); title(ax, 'Battery Charging Power During Charging', 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt); lgd(ax, algName(1:3), 'best', fontName, fontSize - 4, fontWt);
savePlot(fig, outPlots, '12b_Battery_Charging_Power');

%% STEP 13 - UAV PROPULSION LOAD + INTEGRATED MISSION DEMO
Pto = uavElecPower(uav.tauTO, uav); Phv = uavElecPower(uav.tauHov, uav); Pcr = uavElecPower(uav.tauCruise, uav);
fprintf('\nSTEP 13 - UAV : mass %.1f kg | Takeoff = %.0f W | Hover = %.0f W | Cruise = %.0f W\n', uav.mass, Pto, Phv, Pcr);
Tm = 3.0; Nm = round(Tm/mp.dt); tm = (0:Nm-1)' * mp.dt;
Gm = 850*ones(Nm, 1); Tcm = cellTemp(Gm, 30, 1, PV, gT9); plM = makeUniformPlant(Gm, Tcm, PV, aS, aP);
tau = uav.tauCruise*ones(Nm, 1); tau(tm < Tm/6) = uav.tauTO; tau(tm >= Tm/6 & tm < Tm/2) = uav.tauHov;
Plm = uavElecPower(tau, uav);
cfgM = cfg0; cfgM.N = Nm; cfgM.useBat = true; cfgM.Pload = Plm; oM = simMPPT(plM, 3, smcOpt, cfgM);
PpvM = oM.V .* oM.I;
EpvWh = trapz(tm, PpvM)/3600 * mp.socAccel; ElWh = trapz(tm, Plm)/3600 * mp.socAccel;
fprintf('Mission demo (G = 850 W/m^2, SOC time scale x%d): PV energy = %.1f Wh | load energy = %.1f Wh | SOC %.4f -> %.4f\n', mp.socAccel, EpvWh, ElWh, oM.SOC(1), oM.SOC(end));
DemoTbl = table(["PV_Energy_Wh";"Load_Energy_Wh";"SOC_Start";"SOC_End"], [EpvWh; ElWh; oM.SOC(1); oM.SOC(end)], 'VariableNames', {'Parameter','Value'});

%% STEP 15 - 10-STATION MPPT COMPARISON (representative day compressed to mp.dayDur seconds)
EffM = NaN(nSt, 3); ConvM = NaN(nSt, 3); RipM = NaN(nSt, 3); TrkM = NaN(nSt, 3);
fprintf('\n============================================================\nSTEP 15 - 10-STATION MPPT COMPARISON\n============================================================\n');
for i = 1:nSt
    best = RepBest(i);
    if isnan(best) || isempty(R5{i}), fprintf('%s skipped (no representative day / PV model).\n', Stations(i).Name); continue; end
    [~, Gd, Tad, wsd] = repDayProfile(Store(i), best);
    idx = find(Gd > daylightThreshold); if numel(idx) < 10, fprintf('%s skipped (too few daylight samples).\n', Stations(i).Name); continue; end
    ix = idx(1):idx(end);
    [plS, PrefS] = buildScenario((1:numel(ix))', effG(Gd(ix), R5{i}.par.beta), Tad(ix), wsd(ix), mp.dayDur, mp.dt, PV, R5{i}.par.gT, aS, aP);
    cfgS = cfg0; cfgS.N = plS.N;
    for a = 1:3
        o = runAlg(plS, a, cfgS, smcOpt, mp); mt = mpptMetrics(o.t, o.V .* o.I, PrefS);
        EffM(i,a) = mt.Eff; ConvM(i,a) = mt.Conv_ms; RipM(i,a) = mt.RipplePct; TrkM(i,a) = mt.TrkErrPct;
    end
    fprintf('%-10s | Efficiency: P&O = %.2f %% | INC = %.2f %% | GA-SMC = %.2f %%\n', Stations(i).Name, EffM(i,1), EffM(i,2), EffM(i,3));
end
fprintf('Mean efficiency : P&O = %.2f %% | INC = %.2f %% | GA-SMC = %.2f %%\n', mean(EffM(:,1), 'omitnan'), mean(EffM(:,2), 'omitnan'), mean(EffM(:,3), 'omitnan'));
fprintf('Mean ripple     : P&O = %.3f %% | INC = %.3f %% | GA-SMC = %.3f %%\n', mean(RipM(:,1), 'omitnan'), mean(RipM(:,2), 'omitnan'), mean(RipM(:,3), 'omitnan'));
fprintf('Mean convergence: P&O = %.1f ms | INC = %.1f ms | GA-SMC = %.1f ms (same 20 ms-sustained 95 %% definition as everywhere else)\n', mean(ConvM(:,1), 'omitnan'), mean(ConvM(:,2), 'omitnan'), mean(ConvM(:,3), 'omitnan'));
gaBetter = EffM(:,3) > max(EffM(:,1:2), [], 2);
WinTbl = table(string(stNames(:)), EffM(:,3) - EffM(:,1), EffM(:,3) - EffM(:,2), gaBetter, 'VariableNames', {'Station','dEff_GASMC_minus_PO_pct','dEff_GASMC_minus_INC_pct','GASMC_Best_of_Three'});
fprintf('\nPER-STATION WIN/LOSS (efficiency)\n'); disp(WinTbl);
badSt = find(~gaBetter & ~any(isnan(EffM), 2))';
if isempty(badSt), fprintf('GA-SMC has the highest efficiency at every station.\n');
else, fprintf('GA-SMC is NOT the best at: %s -> report as station-dependent.\n', strjoin(stNames(badSt), ', ')); end
hold_out = setdiff(1:nSt, mp.trStations);
fprintf('Held-out stations (not used in GA training) : mean eff P&O = %.2f %% | INC = %.2f %% | GA-SMC = %.2f %%\n', mean(EffM(hold_out,1),'omitnan'), mean(EffM(hold_out,2),'omitnan'), mean(EffM(hold_out,3),'omitnan'));

%% PLOT 7 - MPPT CONVERGENCE TIME COMPARISON
fig = barGrouped(ConvM, stNames, algCol(1:3,:), 'MPPT Convergence Time Comparison', 'Station', 'Convergence Time to 95 % of MPP (ms)', algName(1:3), fontName, fontSize, fontWt, '%.1f', 0, false);
savePlot(fig, outPlots, '07_MPPT_Convergence_Time_Comparison');

%% STEP 14 - FLIGHT ENDURANCE
dtMin = 1; nEnd = mp.maxH*60; tHr = (0:nEnd-1)' * dtMin/60;
EnduTbl = NaN(nSt, 5); SolE = NaN(nSt, 1); caseName = {'No solar', 'P&O', 'INC', 'GA-SMC', 'Ideal MPPT'};
fprintf('\n============================================================\nSTEP 14 - FLIGHT ENDURANCE (max simulated %.0f h, simulation-based energy-management evaluation)\n============================================================\n', mp.maxH);
for i = 1:nSt
    best = RepBest(i);
    if isnan(best) || isempty(R5{i}) || any(isnan(EffM(i,:))), continue; end
    [hrs, Gd, Tad, wsd] = repDayProfile(Store(i), best);
    tod = mod(mp.startH + tHr, 24); [Gt, Tat, wst] = todProfile(hrs, Gd, Tad, wsd, tod);
    Gt = effG(Gt, R5{i}.par.beta);
    Tct = cellTemp(Gt, Tat, wst, PV, R5{i}.par.gT); PrefT = mppSD(Gt, Tct, PV) * aS * aP;
    scl = [0, EffM(i,1:3)/100, 1];
    for c = 1:5
        tE = enduranceSim(scl(c) * mp.etaConv * PrefT, dtMin, mp, uav);
        EnduTbl(i,c) = tE;
    end
    SolE(i) = sum(scl(4) * mp.etaConv * PrefT(1:min(1440, numel(PrefT)))) * dtMin/60;
    fprintf('%-10s | Endurance (h): no solar = %.2f | P&O = %.2f | INC = %.2f | GA-SMC = %.2f | ideal = %.2f | GA-SMC gain vs P&O = %+.1f min\n', Stations(i).Name, EnduTbl(i,:), 60*(EnduTbl(i,4) - EnduTbl(i,2)));
end

%% PLOT 13 - DRONE FLIGHT ENDURANCE COMPARISON
cEn = [0.5 0.5 0.5; algCol(1,:); algCol(2,:); algCol(3,:); 0 0 0];
fig = barGrouped(EnduTbl, stNames, cEn, 'Drone Flight Endurance Comparison', 'Station', 'Flight Endurance (h)', caseName, fontName, fontSize, fontWt, '%.2f', 0, false);
savePlot(fig, outPlots, '13_Drone_Flight_Endurance_Comparison');

%% STEP 16 - PARTIAL SHADING
Gsh = [1000 700 300]; Tcsh = cellTemp(Gsh, 25, 1, PV, gT9);
plSh = makeShadedPlant(Gsh, Tcsh, PV, 1);
Vl = (0:plSh.nLut-1)' * plSh.dV; Plut = Vl .* plSh.lut; [PGmpp, jG] = max(Plut); VGmpp = Vl(jG); pk = findPeaksSimple(Vl, Plut);
fprintf('\n============================================================\nSTEP 16 - PARTIAL SHADING (G = %s W/m^2)\n============================================================\n', mat2str(Gsh));
fprintf('Number of power peaks = %d | GMPP = %.1f W at %.1f V\n', numel(pk), PGmpp, VGmpp);
Tsh = 1.0; cfgSh = cfg0; cfgSh.N = round(Tsh/mp.dt); oSh = cell(1,4); shTbl = table(); PrefSh = PGmpp * ones(cfgSh.N, 1);
for a = 1:4
    oSh{a} = runAlg(plSh, a, cfgSh, smcOpt, mp); Psh = oSh{a}.V .* oSh{a}.I; nl = round(0.2/mp.dt); Pfin = mean(Psh(end-nl+1:end)); Vfin = mean(oSh{a}.V(end-nl+1:end));
    shTbl = [shTbl; table(string(algName{a}), Pfin, 100*Pfin/PGmpp, Vfin, 100*Pfin/PGmpp >= 98, oSh{a}.nScan, 'VariableNames', {'Method','FinalPower_W','PercentOfGMPP','FinalVoltage_V','ReachedGMPP','GlobalScans'})]; %#ok<AGROW>
end
fprintf('\nPARTIAL-SHADING RESULT (GMPP = %.1f W)\n', PGmpp); disp(shTbl);
fprintf('INTERPRETATION: GA-SMC alone (local tracker) reaches %.2f %% of GMPP; the architecture\n  "GA-optimised sliding-mode MPPT + global GMPP scan layer" reaches %.2f %% of GMPP.\n', shTbl.PercentOfGMPP(3), shTbl.PercentOfGMPP(4));

%% PLOT 6 - GLOBAL MAXIMUM POWER POINT TRACKING RESPONSE
fig = newFig(); ax = axes(fig); hold(ax, 'on'); plot(ax, oSh{1}.t, PrefSh, 'k--', 'LineWidth', 2.5);
for a = 1:4, plot(ax, oSh{a}.t, oSh{a}.V .* oSh{a}.I, 'Color', algCol(a,:), 'LineWidth', 2); end
xlabel(ax, 'Time (s)'); ylabel(ax, 'PV Power (W)'); title(ax, 'Global Maximum Power Point Tracking Response', 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
lgd(ax, [{'GMPP'}, algName(1:4)], 'southeast', fontName, fontSize - 4, fontWt);
savePlot(fig, outPlots, '06_Global_MPP_Tracking_Response');

%% PLOT 10 - MPPT PERFORMANCE UNDER PARTIAL SHADING
fig = newFig(); ax = axes(fig); hold(ax, 'on'); plot(ax, Vl, Plut, 'k', 'LineWidth', 2.5);
for a = 1:4, plot(ax, oSh{a}.V(1:20:end), oSh{a}.V(1:20:end) .* oSh{a}.I(1:20:end), '.', 'Color', algCol(a,:), 'MarkerSize', 6); end
plot(ax, Vl(pk), Plut(pk), 'kv', 'MarkerFaceColor', 'y', 'MarkerSize', 12); plot(ax, VGmpp, PGmpp, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 16);
xlabel(ax, 'Array Voltage (V)'); ylabel(ax, 'Array Power (W)'); title(ax, 'MPPT Performance Under Partial Shading - P-V Curve and Operating Trajectories', 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
lgd(ax, [{'P-V curve'}, algName(1:4), {'Local peaks', 'GMPP'}], 'northeast', fontName, fontSize - 6, fontWt);
savePlot(fig, outPlots, '10a_MPPT_Partial_Shading_PV_Trajectories');
fig = barSingle(shTbl.PercentOfGMPP, algName(1:4), 'MPPT Performance Under Partial Shading - Tracking Accuracy', 'MPPT Method', 'Final Power (% of GMPP)', algCol(1:4,:), fontName, fontSize, fontWt, '%.2f');
ax = fig.CurrentAxes; hold(ax, 'on'); yline(ax, 98, 'k--', 'LineWidth', 2);
savePlot(fig, outPlots, '10b_MPPT_Partial_Shading_Tracking_Accuracy');

%% STEP 17 - ABLATION STUDY (P&O | SMC only | GA-SMC) on scenarios A, B, C (from Step 11)
ablM = [1 5 3]; abl = table(); ablVal = NaN(3, 3, 4);
for n = 1:3
    for q = 1:3
        a = ablM(q); ablVal(n,q,:) = [DynConv(n,a) DynEff(n,a) DynRip(n,a) DynTrk(n,a)];
        abl = [abl; table(string(Sc{n}.Name), string(algName{a}), DynConv(n,a), DynEff(n,a), DynRip(n,a), DynTrk(n,a), 'VariableNames', {'Scenario','Method','Convergence_ms','Efficiency_pct','Ripple_pct','TrackingError_pct'})]; %#ok<AGROW>
    end
end
fprintf('\n============================================================\nSTEP 17 - ABLATION STUDY\n============================================================\n'); disp(abl);
mAbl = squeeze(mean(ablVal, 1, 'omitnan'));
AblMean = table(string(algName(ablM))', mAbl(:,1), mAbl(:,2), mAbl(:,3), mAbl(:,4), 'VariableNames', {'Method','Convergence_ms','Efficiency_pct','Ripple_pct','TrackingError_pct'});
fprintf('MEAN OVER SCENARIOS A-C\n'); disp(AblMean);
AblDelta = table(mAbl(3,2) - mAbl(2,2), mAbl(3,4) - mAbl(2,4), mAbl(3,3) - mAbl(2,3), mAbl(3,1) - mAbl(2,1), 'VariableNames', {'dEfficiency_pct','dTrackingError_pct','dRipple_pct','dConvergence_ms'});
fprintf('GA-SMC minus SMC-only : efficiency %+.2f %% | tracking error %+.2f %% | ripple %+.2f %% | convergence %+.1f ms\n', AblDelta.dEfficiency_pct, AblDelta.dTrackingError_pct, AblDelta.dRipple_pct, AblDelta.dConvergence_ms);
if AblDelta.dEfficiency_pct > 0, fprintf('GA optimisation raises mean efficiency by %.2f points over SMC-only on scenarios A-C.\n', AblDelta.dEfficiency_pct);
else, fprintf('CLAIM GUARD: GA optimisation does not raise mean efficiency over SMC-only on scenarios A-C (%+.2f points).\n', AblDelta.dEfficiency_pct); end
if AblDelta.dConvergence_ms > 1, fprintf('CLAIM GUARD: GA-SMC is %.1f ms SLOWER than SMC-only -> write "comparable convergence", not "faster".\n', AblDelta.dConvergence_ms);
elseif AblDelta.dConvergence_ms < -1, fprintf('GA-SMC convergence is %.1f ms faster than SMC-only on scenarios A-C.\n', -AblDelta.dConvergence_ms);
else, fprintf('GA-SMC and SMC-only convergence differ by less than the 1 ms sampling time (comparable).\n'); end
if AblDelta.dRipple_pct > 0, fprintf('CLAIM GUARD: GA-SMC ripple is higher than SMC-only -> report as a trade-off.\n');
else, fprintf('GA-SMC ripple is %.3f points lower than SMC-only on A-C.\n', -AblDelta.dRipple_pct); end

%% STEP 18 - MPPT EFFICIENCY UNDER DIFFERENT TEMPERATURE LEVELS (G profile 1000 -> 400 -> 800 W/m^2)
tmpLv = mp.tempLevels(:)'; nTL = numel(tmpLv); EffT = NaN(nTL, 3); PmeanT = NaN(nTL, 1); TempTbl = table();
NT = round(mp.tempDur/mp.dt); tT = (0:NT-1)' * mp.dt; segIdxT = min(floor(tT/(mp.tempDur/3)) + 1, 3); GT = segG(segIdxT);
ssMask = false(NT, 1); for sg = 1:3, ii = find(segIdxT == sg); ssMask(ii(round(0.5*numel(ii))+1:end)) = true; end
fprintf('\n============================================================\nSTEP 18 - MPPT TRACKING EFFICIENCY VS TEMPERATURE (steady state)\n============================================================\n');
for q = 1:nTL
    TcT = tmpLv(q)*ones(NT,1); PrefT = mppSD(GT, TcT, PV) * aS * aP; PmeanT(q) = mean(PrefT);
    plT = makeUniformPlant(GT, TcT, PV, aS, aP); cfgT = cfg0; cfgT.N = NT;
    for a = 1:3
        o = runAlg(plT, a, cfgT, smcOpt, mp); PT = o.V .* o.I; mt = mpptMetrics(o.t, PT, PrefT);
        EffT(q,a) = 100*sum(PT(ssMask))/sum(PrefT(ssMask));
        TempTbl = [TempTbl; table(tmpLv(q), string(algName{a}), EffT(q,a), mt.Eff, mt.Conv_ms, mt.RipplePct, mt.TrkErrPct, mean(PrefT), 'VariableNames', {'CellTemperature_C','Method','SteadyStateEfficiency_pct','OverallEfficiency_pct','Convergence_ms','Ripple_pct','TrackingError_pct','MeanReferenceMPP_W'})]; %#ok<AGROW>
    end
    fprintf('Tc = %2d C | Steady-state MPPT efficiency: P&O = %.2f %% | INC = %.2f %% | GA-SMC = %.2f %%\n', tmpLv(q), EffT(q,1), EffT(q,2), EffT(q,3));
end

%% PLOT 9 - MPPT EFFICIENCY UNDER DIFFERENT TEMPERATURE LEVELS
fig = newFig(); ax = axes(fig); hold(ax, 'on'); mk = {'-o', '-s', '-^'};
for a = 1:3, plot(ax, tmpLv, EffT(:,a), mk{a}, 'Color', algCol(a,:), 'LineWidth', 2.5, 'MarkerSize', 9, 'MarkerFaceColor', algCol(a,:)); end
xlabel(ax, 'Cell Temperature (°C)'); ylabel(ax, 'Steady-State MPPT Efficiency (%)'); title(ax, 'MPPT Tracking Efficiency Under Different Temperature Levels', 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt); lgd(ax, algName(1:3), 'best', fontName, fontSize - 4, fontWt);
savePlot(fig, outPlots, '09a_MPPT_Efficiency_Under_Different_Temperature_Levels');
fig = newFig(); ax = axes(fig); plot(ax, tmpLv, PmeanT, '-d', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.5, 'MarkerSize', 9, 'MarkerFaceColor', [0.85 0.33 0.10]);
xlabel(ax, 'Cell Temperature (°C)'); ylabel(ax, 'Mean Reference MPP Power (W)'); title(ax, 'Reference MPP Power vs Cell Temperature', 'Interpreter', 'none'); styleAx(ax, fontName, fontSize, fontWt);
savePlot(fig, outPlots, '09b_Reference_MPP_Power_vs_Temperature');

%% TABLES -> ONE EXCEL FILE
fprintf('\n============================================================\nWRITING EXCEL TABLES\n============================================================\n');
algs3 = string(algName(1:3));
vnB = {'Eff_PO_pct','Eff_INC_pct','Eff_GASMC_pct','Conv_PO_ms','Conv_INC_ms','Conv_GASMC_ms','Ripple_PO_pct','Ripple_INC_pct','Ripple_GASMC_pct','TrkErr_PO_pct','TrkErr_INC_pct','TrkErr_GASMC_pct'};
MB = [EffM, ConvM, RipM, TrkM];
BaseTbl = [table(string(stNames(:)), 'VariableNames', {'Station'}), array2table(MB, 'VariableNames', vnB)];
sumRows = {"Mean", mean(MB,1,'omitnan'); "Std", std(MB,0,1,'omitnan'); "Min (worst-case for efficiency)", min(MB,[],1,'omitnan'); "Max", max(MB,[],1,'omitnan')};
for q = 1:size(sumRows,1), BaseTbl = [BaseTbl; [table(sumRows{q,1}, 'VariableNames', {'Station'}), array2table(sumRows{q,2}, 'VariableNames', vnB)]]; end %#ok<AGROW>
MethTbl = table(algs3(:), mean(EffM,1,'omitnan')', mean(ConvM,1,'omitnan')', mean(RipM,1,'omitnan')', mean(TrkM,1,'omitnan')', 'VariableNames', {'Method','MeanEfficiency_pct','MeanConvergence_ms','MeanRipple_pct','MeanTrackingError_pct'});
sh = 'T1_Baseline_MPPT'; r0 = 1;
r0 = putBlock(xlFinal, sh, 'Baseline MPPT Comparison - 10 stations (representative day)', BaseTbl, r0);
r0 = putBlock(xlFinal, sh, 'Baseline MPPT Comparison - method averages', MethTbl, r0);
r0 = putBlock(xlFinal, sh, 'Per-station win/loss of GA-SMC (efficiency)', WinTbl, r0);
putBlock(xlFinal, sh, 'Base training profile response', TrTbl, r0);
sh = 'T2_Ablation_Study'; r0 = 1;
r0 = putBlock(xlFinal, sh, 'Ablation Study - per scenario (P&O | SMC only | GA-SMC)', abl, r0);
r0 = putBlock(xlFinal, sh, 'Ablation Study - mean over scenarios A-C', AblMean, r0);
putBlock(xlFinal, sh, 'GA-SMC minus SMC-only (effect of GA optimisation)', AblDelta, r0);
sh = 'T3_MPP_Extraction'; r0 = 1;
r0 = putBlock(xlFinal, sh, 'MPP Extraction Results - per station (array level, daylight samples)', MPPSumX, r0);
r0 = putBlock(xlFinal, sh, 'Single-diode model check at STC (module)', STCTbl, r0);
putBlock(xlFinal, sh, 'Single-diode model parameters', ParamTbl, r0);
sh = 'T4_Dynamic_Env_Performance'; r0 = 1;
r0 = putBlock(xlFinal, sh, 'Dynamic Environmental Condition Performance - scenarios A-F', DynTbl, r0);
r0 = putBlock(xlFinal, sh, 'Temperature-level performance (G: 1000 -> 400 -> 800 W/m^2)', TempTbl, r0);
putBlock(xlFinal, sh, 'Partial shading - GA-SMC alone vs GA-SMC + global scan layer', shTbl, r0);
DroneTbl = table(string(stNames(:)), EnduTbl(:,1), EnduTbl(:,2), EnduTbl(:,3), EnduTbl(:,4), EnduTbl(:,5), 60*(EnduTbl(:,4) - EnduTbl(:,2)), 60*(EnduTbl(:,4) - EnduTbl(:,3)), ...
    100*(EnduTbl(:,4) - EnduTbl(:,1))./EnduTbl(:,1), 100*EnduTbl(:,4)./EnduTbl(:,5), SolE, ...
    'VariableNames', {'Station','Endurance_NoSolar_h','Endurance_PO_h','Endurance_INC_h','Endurance_GASMC_h','Endurance_Ideal_h','Gain_GASMC_vs_PO_min','Gain_GASMC_vs_INC_min','Gain_GASMC_vs_NoSolar_pct','GASMC_pct_of_Ideal','SolarEnergy_24h_GASMC_Wh'});
Pavg = (Pto*mp.mission(1) + Pcr*mp.mission(2) + Phv*mp.mission(3)) / sum(mp.mission);
SysTbl = table(["UAV_mass_kg";"Rotors";"Rotor_diameter_m";"Takeoff_power_W";"Hover_power_W";"Cruise_power_W";"Mission_average_load_W";"Mission_cycle_takeoff_cruise_hover_min";"Battery_series_cells";"Battery_capacity_Ah";"Battery_energy_kWh";"Start_SOC";"Cutoff_SOC";"PV_array_modules_series";"PV_array_modules_parallel";"PV_array_Pmax_STC_W";"Converter_efficiency_endurance";"Max_simulated_endurance_h"], ...
    [uav.mass; uav.nRot; uav.dia; Pto; Phv; Pcr; Pavg; sum(mp.mission); mp.Ncell; mp.Cb; EnomWh/1e3; mp.SOC0; mp.SOCmin; aS; aP; Vin0*Iin0; mp.etaConv; mp.maxH], 'VariableNames', {'Parameter','Value'});
sh = 'T5_Drone_PowerSystem'; r0 = 1;
r0 = putBlock(xlFinal, sh, 'Simulation-based UAV energy-management evaluation - flight endurance per station', DroneTbl, r0);
r0 = putBlock(xlFinal, sh, 'UAV / battery / PV system parameters', SysTbl, r0);
r0 = putBlock(xlFinal, sh, 'Battery charging test (standby load, time-accelerated SOC model)', ChgTbl, r0);
putBlock(xlFinal, sh, 'Integrated mission demo (takeoff -> hover -> cruise)', DemoTbl, r0);
metList = {'Efficiency_pct', EffM, true; 'Convergence_ms', ConvM, false; 'Ripple_pct', RipM, false; 'TrackingError_pct', TrkM, false};
nCmp = 2*size(metList, 1); StatTbl = table();
for q = 1:size(metList, 1)
    M = metList{q,2};
    for a = 1:3
        s = statCompare(M(:,a), M(:,3), metList{q,3});
        if a == 3
            sg = "-"; ip = "Reference method"; nb = NaN; pB = NaN;
        else
            nb = s(12); pB = min(s(10)*nCmp, 1);
            sg = ifelse(~isnan(s(10)) && s(10) < 0.05, "Yes", "No");
            better = (metList{q,3} && s(9) > 0) || (~metList{q,3} && s(9) < 0);     % direction of the mean difference in favour of GA-SMC
            if isnan(s(10))
                ip = "Not testable";
            elseif s(10) >= 0.05
                if better, ip = "GA-SMC better mean, not significant"; else, ip = "GA-SMC worse mean, not significant"; end
            elseif better
                ip = "GA-SMC significantly better";
            else
                ip = "GA-SMC significantly worse";
            end
        end
        StatTbl = [StatTbl; table(string(metList{q,1}), algs3(a), s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), s(9), s(10), s(11), pB, nb, sg, string(ip), ...
            'VariableNames', {'Metric','Method','N','Mean','Std','Min','Median','Max','CI95_Low','CI95_High','MeanDiff_GASMC_minus_Method','p_ttest','p_signrank','p_ttest_Bonferroni','Stations_GASMC_Better','Significant_p05','Interpretation'})]; %#ok<AGROW>
    end
end
fprintf('\nSTATISTICAL PERFORMANCE COMPARISON\n'); disp(StatTbl);
sh = 'T6_Statistical_Comparison'; r0 = 1;
r0 = putBlock(xlFinal, sh, 'Statistical Performance Comparison - MPPT methods across stations (paired tests vs GA-SMC)', StatTbl, r0);
if height(CalTbl) > 0
    pvVars = {'RMSE_Test','nRMSE_Test_pct','MAE_Test','MBE_Test','R2_Test','CC_Test'}; PVStat = table();
    for q = 1:numel(pvVars)
        x = CalTbl.(pvVars{q});
        PVStat = [PVStat; table(string(pvVars{q}), numel(x), mean(x,'omitnan'), std(x,'omitnan'), min(x), median(x,'omitnan'), max(x), 'VariableNames', {'Metric','N','Mean','Std','Min','Median','Max'})]; %#ok<AGROW>
    end
    putBlock(xlFinal, sh, 'PV model validation statistics across stations (test set)', PVStat, r0);
    writetable(CalTbl, xlFinal, 'Sheet', 'Supp_PV_Validation');
end
% ---- Notes sheet : scaling chain and claim guidance (GA-SMC lines are generated from this run's numbers) ----
notes = {'SCALING CHAIN';
 'China regional PV data (MW-scale utility stations) -> real irradiance / temperature / wind profiles';
 '-> calibrated single-diode PV model (effective, lumped parameters gT, beta, array size) -> representative environmental conditions';
 sprintf('-> scaled reference PV array (%d x %d modules, %.0f W STC) -> MPPT controller -> boost converter -> %d V battery bus -> UAV propulsion load', aS, aP, aS*aP*330, mp.Vbat);
 'The MW-scale stations are NOT installed on the UAV; only their environmental profiles are reused.';
 '';
 'CLAIM GUIDANCE';
 'PV model: good-to-moderate, station-dependent agreement; gT, beta and effective array size are lumped parameters, not measured module properties or installed capacity.';
 sprintf('GA-SMC (this run): mean efficiency P&O %.2f %% | INC %.2f %% | GA-SMC %.2f %% over the 10 stations; GA-SMC best of three at %d of %d stations. Report exactly these numbers.', mean(EffM(:,1),'omitnan'), mean(EffM(:,2),'omitnan'), mean(EffM(:,3),'omitnan'), sum(gaBetter), nSt);
 sprintf('GA training stations: %s (held-out mean efficiency: P&O %.2f %% | INC %.2f %% | GA-SMC %.2f %%).', mat2str(mp.trStations), mean(EffM(hold_out,1),'omitnan'), mean(EffM(hold_out,2),'omitnan'), mean(EffM(hold_out,3),'omitnan'));
 'Ablation: report the GA-SMC minus SMC-only difference exactly as printed in T2 (efficiency, convergence, ripple).';
 'Ripple: report the paired-test result in T6_Statistical_Comparison as it is (higher, equal or lower than P&O/INC).';
 'Partial shading: GA-SMC alone is a local tracker; GA-SMC + global scan layer finds the GMPP.';
 'Battery / UAV: simulation-based energy-management evaluation with a time-accelerated SOC model; not experimental UAV flight validation.'};
writecell(notes, xlFinal, 'Sheet', 'Notes_Scaling_Claims', 'Range', 'A1');

save(fullfile(outTab, 'Final_Results.mat'), 'mp', 'uav', 'cc', 'smcOpt', 'bestHist', 'meanHist', 'CalTbl', 'MPPSumX', 'BaseTbl', 'WinTbl', 'DynTbl', 'TempTbl', 'abl', 'AblMean', 'DroneTbl', 'StatTbl', 'shTbl', 'ChgTbl');
fprintf('\nALL DONE.\nPlots  (PNG) : %s\nTables (Excel)  : %s\n', outPlots, xlFinal);

%% LOCAL FUNCTIONS
function savePlot(fig, folder, name)
    drawnow; exportgraphics(fig, fullfile(folder, [name '.png']), 'Resolution', 300);
end

function r = putBlock(f, sh, ttl, T, r)
    writecell({ttl}, f, 'Sheet', sh, 'Range', sprintf('A%d', r));
    writetable(T, f, 'Sheet', sh, 'Range', sprintf('A%d', r+1));
    r = r + height(T) + 4;
end

function e = convEff(o)
    Pin = sum(o.V .* o.I); Pout = sum(o.Vbat .* (1 - o.D) .* o.iL);
    e = min(100*Pout/max(Pin, eps), 100);
end

function s = statCompare(x, ref, hb)
    x = x(:); ref = ref(:); xv = x(~isnan(x)); n = numel(xv); s = NaN(1, 12);
    if n == 0, return; end
    s(1) = n; s(2) = mean(xv); s(3) = std(xv); s(4) = min(xv); s(5) = median(xv); s(6) = max(xv);
    if n > 1, h = tinv(0.975, n-1) * std(xv) / sqrt(n); s(7) = s(2) - h; s(8) = s(2) + h; else, s(7) = s(2); s(8) = s(2); end
    k = ~isnan(x) & ~isnan(ref);
    if any(k)
        s(9) = mean(ref(k) - x(k));
        if hb, s(12) = sum(ref(k) > x(k)); else, s(12) = sum(ref(k) < x(k)); end
    end
    if sum(k) > 2 && ~isequal(x(k), ref(k))
        try, [~, pp] = ttest(ref(k), x(k)); s(10) = pp; catch, end
        try, s(11) = signrank(ref(k), x(k)); catch, end
    end
end

function fig = barGrouped(M, cats, cols, ttl, xl, yl, lg, fn, fs, fw, fmt, rot, tight)
    % grouped bar chart: one colour per series, value printed on every bar; x tick labels and value labels use rotation 'rot' (0 = horizontal)
    fig = newFig(); ax = axes(fig); b = bar(ax, M, 'grouped');
    for q = 1:numel(b), b(q).FaceColor = cols(q,:); end
    set(ax, 'XTick', 1:numel(cats), 'XTickLabel', cats, 'XTickLabelRotation', 0, 'TickLabelInterpreter', 'none');
    xlabel(ax, xl); ylabel(ax, yl); title(ax, ttl, 'Interpreter', 'none'); styleAx(ax, fn, fs, fw);
    lo = min(M(:), [], 'omitnan'); hi = max(M(:), [], 'omitnan');
    if tight, lo = max(floor(lo) - 2, 0); else, lo = min(0, lo); end
    if ~isfinite(hi) || hi <= lo, hi = lo + 1; end
    ylim(ax, [lo, hi + 0.25*(hi - lo)]);
    addBarValues(ax, b, fmt, fs - 12, fw, rot);
    legend(ax, lg, 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontName', fn, 'FontSize', fs - 4, 'FontWeight', fw, 'Box', 'off', 'Interpreter', 'none');
end

function fig = barSingle(vals, cats, ttl, xl, yl, cols, fn, fs, fw, fmt)
    vals = vals(:)'; n = numel(vals); if isempty(cols), cols = jet(max(n, 2)); end
    fig = newFig(); ax = axes(fig); b = bar(ax, 1:n, vals, 'FaceColor', 'flat'); b.CData = cols(1:n, :);
    set(ax, 'XTick', 1:n, 'XTickLabel', cats, 'XTickLabelRotation', 0, 'TickLabelInterpreter', 'none');
    xlabel(ax, xl); ylabel(ax, yl); title(ax, ttl, 'Interpreter', 'none'); styleAx(ax, fn, fs, fw);
    hi = max(vals, [], 'omitnan'); lo = min(0, min(vals, [], 'omitnan')); if ~isfinite(hi) || hi <= lo, hi = lo + 1; end
    ylim(ax, [lo, hi + 0.15*(hi - lo)]);
    addBarValues(ax, b, fmt, fs - 6, fw, 0);
end

function addBarValues(ax, b, fmt, fs, fw, rot)
    for q = 1:numel(b)
        xe = b(q).XEndPoints; ye = b(q).YEndPoints; v = b(q).YData; k = ~isnan(v) & ~isnan(ye);
        if ~any(k), continue; end
        lab = arrayfun(@(z) sprintf(fmt, z), v(k), 'UniformOutput', false);
        if rot == 0, ha = 'center'; va = 'bottom'; else, ha = 'left'; va = 'middle'; end
        text(ax, xe(k), ye(k), lab, 'HorizontalAlignment', ha, 'VerticalAlignment', va, 'Rotation', rot, 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', fw, 'Interpreter', 'none');
    end
end

function T = readStation(f)
    T = readtable(f, 'VariableNamingRule', 'preserve');
    if ~hasCJK(T.Properties.VariableNames)
        for enc = ["GBK", "GB2312", "UTF-8"]
            try, T2 = readtable(f, 'VariableNamingRule', 'preserve', 'Encoding', char(enc)); catch, continue; end
            if hasCJK(T2.Properties.VariableNames), T = T2; return; end
        end
    end
end

function tf = hasCJK(names)
    tf = any(cellfun(@(s) any(double(char(s)) >= 19968 & double(char(s)) <= 40959), cellstr(string(names))));
end

function [t, ok] = parseTime(d)
    ok = true;
    if isdatetime(d), t = d; return; end
    if isnumeric(d), t = datetime(d, 'ConvertFrom', 'excel'); return; end
    s = strtrim(string(d)); t = NaT(size(s));
    fm = ["yyyy-MM-dd HH:mm:ss","yyyy-MM-dd HH:mm","yyyy/MM/dd HH:mm:ss","yyyy/MM/dd HH:mm","yyyy/M/d H:mm","yyyy/M/d H:mm:ss","yyyy-M-d H:mm","yyyy-M-d H:mm:ss","M/d/yyyy H:mm","d/M/yyyy H:mm","MM/dd/yyyy HH:mm","dd/MM/yyyy HH:mm","yyyyMMddHHmm","yyyyMMdd HHmmss"];
    for f = fm
        try, tt = datetime(s, 'InputFormat', char(f)); catch, continue; end
        if mean(isnat(tt)) < 0.01, t = tt; return; end
    end
    try, tt = datetime(s); if mean(isnat(tt)) < 0.01, t = tt; return; end, catch, end
    ok = false;
end

function [col, how] = detectCols(Num, vn, ln, pats)
    nV = numel(vn); col = zeros(1,5); how = strings(1,5);
    valid = cellfun(@(z) sum(~isnan(z)) > 10, Num);
    isT = contains(ln, ["date","时间","日期"]) | (contains(ln, "time") & ~contains(ln, "temp"));
    bad = ["direction","风向","angle","pressure","气压"];
    meas = ["lmd","实测","measured","actual","实际","real"]; nwp = ["nwp","forecast","predict","预报","预测"];
    for k = 1:5
        best = 0; bs = -inf;
        for q = 1:nV
            if ~valid(q) || isT(q) || any(col == q) || any(contains(ln(q), bad)), continue; end
            if ~any(contains(vn(q), pats{k}) | contains(ln(q), lower(pats{k}))), continue; end
            s = 2*any(contains(ln(q), meas)) - 2*any(contains(ln(q), nwp)) + (k == 1)*any(contains(ln(q), ["total","总"]));
            if s > bs, bs = s; best = q; end
        end
        if best > 0, col(k) = best; how(k) = "name"; end
    end
    for k = [1 3 5 2 4]
        if col(k) > 0, continue; end
        for q = 1:nV
            if ~valid(q) || isT(q) || any(col == q) || any(contains(ln(q), bad)), continue; end
            y = Num{q}(~isnan(Num{q})); mx = max(y); mn = min(y); mu = mean(y); sd = std(y); zf = mean(y <= 0.001);
            switch k
                case 1, ok = mx >= 300 && mx <= 2000 && zf > 0.2;
                case 3, ok = mn >= 0 && mx < 300 && zf > 0.2;
                case 5, ok = mn >= 0 && mx > 60 && mx <= 100.5 && zf < 0.05 && sd > 3;
                case 2, ok = mn >= -50 && mx <= 70 && zf < 0.1 && sd > 1;
                case 4, ok = mn >= 0 && mx <= 60 && mu < 15 && sd > 0.3;
            end
            if ok, col(k) = q; how(k) = "data"; break; end
        end
    end
end

function [best, F, ud] = pickRepDay(Y, dayId, thr)
    best = NaN; F = []; dayId = dayId(:); ud = 1:max(dayId);
    if isempty(Y{1}), return; end
    cnt = accumarray(dayId, 1)'; F = NaN(numel(ud), 7);
    for u = ud
        m = dayId == u; g = Y{1}(m); dg = g(g > thr);
        if numel(dg) < 0.15*max(cnt) || cnt(u) < 0.9*max(cnt), continue; end
        row = NaN(1,7); row(1) = mean(dg); row(2) = max(g); row(3) = std(dg);
        for k = 2:5, if ~isempty(Y{k}), row(k+2) = mean(Y{k}(m), 'omitnan'); end, end
        F(u,:) = row;
    end
    ok = ~isnan(F(:,1)); if ~any(ok), return; end
    mu = median(F(ok,:), 1, 'omitnan'); sg = std(F(ok,:), 0, 1, 'omitnan'); sg(sg == 0 | isnan(sg)) = 1;
    d = sqrt(sum(((F - mu) ./ sg).^2, 2, 'omitnan')); d(~ok) = Inf;
    [~, best] = min(d);
end

function S = makeSc(nm, src, t, G, Tc, P, W, H, tScale, tLab)
    S.Name = nm; S.Source = src; S.Time_s = t(:); S.Irradiance = G(:); S.Temperature = Tc(:); S.PVPower = P(:); S.WindSpeed = W(:); S.Humidity = H(:);
    S.tPlot = t(:) * tScale; S.tLabel = tLab;
    S.IrradianceTS = timeseries(G(:), t(:)); S.TemperatureTS = timeseries(Tc(:), t(:));
end

function x = toNum(x)
    if isdatetime(x) || isduration(x), x = NaN(numel(x), 1); return; end
    if ~isnumeric(x), x = str2double(string(x)); end
    x = double(x(:));
end

function fig = newFig()
    fig = figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'WindowState', 'maximized');
end

function styleAx(ax, fn, fs, fw)
    set(ax, 'FontName', fn, 'FontSize', fs, 'FontWeight', fw, 'LineWidth', 1.5, 'Box', 'on'); grid(ax, 'off');
    ax.XLabel.FontName = fn; ax.XLabel.FontSize = fs; ax.XLabel.FontWeight = fw;
    ax.YLabel.FontName = fn; ax.YLabel.FontSize = fs; ax.YLabel.FontWeight = fw;
    ax.Title.FontName = fn; ax.Title.FontSize = fs; ax.Title.FontWeight = fw;
end

function v = ifelse(c, a, b)
    if c, v = a; else, v = b; end
end

function S = fitSD(S)
    A = S.n * S.Ns * 8.617333e-5 * 298.15; bestF = inf; bestX = [log(0.3) log(300)];
    opt = optimset('TolX', 1e-10, 'TolFun', 1e-16, 'MaxIter', 4000, 'MaxFunEvals', 8000, 'Display', 'off');
    for r0 = [0.1 0.3 0.6]
        for h0 = [100 500 1500]
            [x, fv] = fminsearch(@(x) sdRes(x, S, A), [log(r0) log(h0)], opt);
            if fv < bestF, bestF = fv; bestX = x; end
        end
    end
    [~, S.Rs, S.RshRef, S.IphRef, S.I0Ref] = sdRes(bestX, S, A);
    if bestF > 1e-3, warning('Single-diode fit residual is large (%.2e). Check the module datasheet values.', bestF); end
end

function [f, Rs, Rsh, Iph, I0] = sdRes(x, S, A)
    Rs = exp(x(1)); Rsh = min(exp(x(2)), 5000);
    E1 = exp(S.Isc*Rs/A) - 1; E2 = exp(S.Voc/A) - 1;
    sol = [1 -E1; 1 -E2] \ [S.Isc*(1 + Rs/Rsh); S.Voc/Rsh]; Iph = sol(1); I0 = sol(2);
    Vd = S.Vmp + S.Imp*Rs; E = exp(Vd/A);
    r1 = Iph - S.Imp - I0*(E - 1) - Vd/Rsh;
    dIdV = -(I0/A*E + 1/Rsh) / (1 + I0*Rs/A*E + Rs/Rsh);
    r2 = S.Imp + S.Vmp*dIdV;
    f = (r1/S.Imp)^2 + (r2/S.Imp)^2;
    if I0 <= 0 || Iph <= 0, f = f + 1e3; end
end

function [Iph, I0, Rs, Rsh, A] = sdParams(G, Tc, S)
    kB = 8.617333e-5; Tk = Tc + 273.15; Tref = 298.15;
    A    = S.n * S.Ns * kB .* Tk;
    Iph  = (S.IphRef + S.Ki*(Tc - 25)) .* G / 1000;
    I0   = S.I0Ref .* (Tk/Tref).^3 .* exp(S.Eg/(S.n*kB) .* (1/Tref - 1./Tk));
    Rs   = S.Rs;
    Rsh  = S.RshRef * 1000 ./ max(G, 1);
end

function I = solveI(V, Iph, I0, Rs, Rsh, A)
    I = Iph + 0*V;
    for it = 1:35
        E  = exp(min((V + I*Rs)./A, 100));
        fI = Iph - I - I0.*(E - 1) - (V + I*Rs)./Rsh;
        df = -1 - I0.*Rs./A.*E - Rs./Rsh;
        I  = I - fI./df;
    end
end

function [Pm, Vm, Im] = mppSD(G, Tc, S)
    G = G(:); Tc = Tc(:); n = numel(G); Pm = zeros(n,1); Vm = zeros(n,1); Im = zeros(n,1);
    for a = 1:4000:n
        idx = a:min(a+3999, n); g = max(G(idx), 1); tc = Tc(idx);
        [Iph, I0, Rs, Rsh, A] = sdParams(g, tc, S);
        Voc = max(S.Voc + S.Kv*(tc - 25) + A.*log(g/1000), 1);
        V = Voc .* linspace(0, 1.05, 120);
        I = max(solveI(V, Iph, I0, Rs, Rsh, A), 0); P = V .* I;
        [pm, j] = max(P, [], 2); lin = sub2ind(size(P), (1:numel(idx))', j);
        vm = V(lin); im = I(lin);
        r = (1:numel(idx))'; jm = min(max(j, 2), size(P,2)-1);
        p1 = P(sub2ind(size(P), r, jm-1)); p2 = P(sub2ind(size(P), r, jm)); p3 = P(sub2ind(size(P), r, jm+1));
        off = 0.5*(p1 - p3) ./ (p1 - 2*p2 + p3); off(~isfinite(off) | abs(off) > 1) = 0;
        Vr = V(sub2ind(size(V), r, jm)) + off .* (V(sub2ind(size(V), r, jm+1)) - V(sub2ind(size(V), r, jm)));
        Ir = max(solveI(Vr, Iph, I0, Rs, Rsh, A), 0); Pr = Vr .* Ir; bt = Pr > pm;
        pm(bt) = Pr(bt); vm(bt) = Vr(bt); im(bt) = Ir(bt);
        Pm(idx) = pm; Vm(idx) = vm; Im(idx) = im;
    end
    low = G <= 1; Pm(low) = 0; Vm(low) = 0; Im(low) = 0;
end

function Tc = cellTemp(G, Ta, ws, S, gT)
    Tc = Ta + gT * G ./ (S.U0 + S.U1*ws);
end

function Ge = effG(G, beta)
    % effective (array-plane) irradiance: Ge = 1000*(G/1000)^beta  (beta = 1 -> no change; STC-anchored at 1000 W/m^2)
    Ge = 1000 * (max(G, 0)/1000).^beta;
end

function [Vpv, Ipv, Pdc] = simArray(G, Ta, ws, S, gT, Np)
    n = numel(G); Vpv = zeros(n,1); Ipv = zeros(n,1); Pdc = zeros(n,1);
    Tc = cellTemp(G, Ta, ws, S, gT); act = G > 1;
    if any(act)
        [~, vm, im] = mppSD(G(act), Tc(act), S);
        Vpv(act) = vm * S.NserString; Ipv(act) = im * Np; Pdc(act) = Vpv(act) .* Ipv(act) / 1e6;
    end
end

function Pu = unitPower(p, G, Ta, ws, S)
    % MPP power (MW) of ONE series string for parameters p = [gT beta]
    Ge = effG(G, p(2)); Tc = cellTemp(Ge, Ta, ws, S, p(1));
    Pu = mppSD(Ge, Tc, S) * S.NserString / 1e6;
end

function [c, k] = calCost(p, G, Ta, ws, Pm, S, cc)
    % normalised squared misfit + weak physical priors; k = eta*Np_effective solved in closed form (least squares)
    Pu = unitPower(p, G, Ta, ws, S);
    k = (Pu' * Pm) / max(Pu' * Pu, eps);
    nm = sum((Pm - k*Pu).^2) / max(sum(Pm.^2), eps);
    c = nm + cc.w * (((p(1) - cc.gT0)/cc.sgT)^2 + ((p(2) - cc.b0)/cc.sb)^2);
end

function [p, k, status] = calibrate(G, Ta, ws, Pm, S, cc)
    % 1) coarse grid  2) smooth-bounded fminsearch (sigmoid transform)  3) fallback: parameters that end at a bound are fixed to physical values
    lo = [cc.gT(1) cc.beta(1)]; hi = [cc.gT(2) cc.beta(2)]; tol = 0.02*(hi - lo);
    cost = @(q) calCost(q, G, Ta, ws, Pm, S, cc);
    best = inf; p0 = [cc.gT0 cc.b0];
    for a = linspace(lo(1), hi(1), 7)
        for b = linspace(lo(2), hi(2), 5)
            c = cost([a b]); if c < best, best = c; p0 = [a b]; end
        end
    end
    map = @(u) lo + (hi - lo) ./ (1 + exp(-u));
    q0 = min(max((p0 - lo)./(hi - lo), 0.02), 0.98); u0 = log(q0 ./ (1 - q0));
    u = fminsearch(@(u) cost(map(u)), u0, optimset('TolX', 1e-3, 'TolFun', 1e-7, 'MaxFunEvals', 80, 'MaxIter', 60, 'Display', 'off'));
    p = map(u); status = "Fitted";
    atB = (p < lo + tol) | (p > hi - tol);
    ob = optimset('TolX', 1e-3, 'Display', 'off');
    if atB(1)                                    % gT not identifiable -> fix it, refit beta
        p(1) = cc.gT0; status = "gT fixed";
        p(2) = fminbnd(@(b) cost([p(1) b]), lo(2), hi(2), ob);
        if (p(2) < lo(2) + tol(2)) || (p(2) > hi(2) - tol(2)), p(2) = 1; status = "gT+beta fixed"; end
    elseif atB(2)                                % beta not identifiable -> fix it to 1, refit gT
        p(2) = 1; status = "beta fixed";
        p(1) = fminbnd(@(g) cost([g 1]), lo(1), hi(1), ob);
        if (p(1) < lo(1) + tol(1)) || (p(1) > hi(1) - tol(1)), p(1) = cc.gT0; status = "gT+beta fixed"; end
    end
    [~, k] = cost(p);
end

function m = metrics(y, yh, cap)
    y = y(:); yh = yh(:); k = ~isnan(y) & ~isnan(yh); y = y(k); yh = yh(k);
    m = struct('N', numel(y), 'RMSE', NaN, 'MAE', NaN, 'MBE', NaN, 'nRMSE', NaN, 'R2', NaN, 'CC', NaN);
    if numel(y) < 3, return; end
    e = yh - y; m.RMSE = sqrt(mean(e.^2)); m.MAE = mean(abs(e)); m.MBE = mean(e); m.nRMSE = 100*m.RMSE/cap;
    m.R2 = 1 - sum(e.^2)/sum((y - mean(y)).^2); c = corrcoef(y, yh); m.CC = c(1,2);
end

function lgd(ax, lb, loc, fn, fs, fw)
    legend(ax, lb, 'Location', loc, 'FontName', fn, 'FontSize', fs, 'FontWeight', fw, 'Box', 'off', 'Interpreter', 'none');
end

function pl = makeUniformPlant(G, Tc, S, aS, aP)
    G = max(G(:), 1e-3); Tc = Tc(:) + 0*G;
    [Iph, I0, Rs, Rsh, A] = sdParams(G, Tc, S);
    pl = struct('mode', 1, 'nS', aS, 'nP', aP, 'N', numel(G), 'Iph', Iph(:), 'I0', I0(:), 'Rs', Rs, 'Rsh', Rsh(:), 'A', A(:), 'lut', [], 'dV', 1, 'nLut', 0);
end

function pl = makeShadedPlant(Gm, Tcm, S, aP)
    nS = numel(Gm); vb = 0.7; Isc = zeros(1, nS); Vc = cell(1, nS); Ic = cell(1, nS);
    for j = 1:nS
        [Iph, I0, Rs, Rsh, A] = sdParams(max(Gm(j), 1), Tcm(j), S);
        Voc = max(S.Voc + S.Kv*(Tcm(j) - 25) + A*log(max(Gm(j), 1)/1000), 1);
        V = linspace(0, 1.05*Voc, 600)'; I = max(solveI(V, Iph, I0, Rs, Rsh, A), 0);
        Vc{j} = V; Ic{j} = I; Isc(j) = I(1);
    end
    Ig = linspace(0, max(Isc)*1.001, 4000)'; Va = zeros(size(Ig));
    for j = 1:nS
        [Iu, iu] = unique(Ic{j}, 'last'); Vu = Vc{j}(iu);
        vj = interp1(Iu, Vu, Ig, 'linear', -vb);
        vj(Ig > Isc(j)) = -vb;
        Va = Va + max(vj, -vb);
    end
    [Vs, so] = sort(Va); Is = Ig(so); [Vs, iu2] = unique(Vs); Is = Is(iu2);
    nLut = 6001; Vl = linspace(0, Vs(end), nLut)'; Il = max(interp1(Vs, Is, Vl, 'linear', 0), 0) * aP;
    pl = struct('mode', 2, 'nS', nS, 'nP', aP, 'N', 0, 'Iph', [], 'I0', [], 'Rs', 0, 'Rsh', [], 'A', [], 'lut', Il, 'dV', Vl(2) - Vl(1), 'nLut', nLut);
end

function [I, Iw] = pvArr(V, k, pl, Iw)
    if pl.mode == 1
        Vm = V / pl.nS; Iph = pl.Iph(k); I0 = pl.I0(k); Rs = pl.Rs; Rsh = pl.Rsh(k); A = pl.A(k); I = Iw;
        for it = 1:25
            Ex = exp(min((Vm + I*Rs)/A, 100));
            f  = Iph - I - I0*(Ex - 1) - (Vm + I*Rs)/Rsh;
            df = -1 - I0*Rs/A*Ex - Rs/Rsh;
            In = I - f/df;
            if abs(In - I) < 1e-9, I = In; break; end
            I = In;
        end
        Iw = I; I = max(I, 0) * pl.nP;
    else
        x = V / pl.dV + 1;
        if V <= 0, I = pl.lut(1); elseif x >= pl.nLut, I = 0; else, i0 = floor(x); fr = x - i0; I = pl.lut(i0)*(1 - fr) + pl.lut(i0+1)*fr; end
    end
end

function [pl, Pref, G, Tc] = buildScenario(tsrc, Gsrc, Tasrc, wssrc, T, dt, S, gT, aS, aP)
    N = round(T/dt); ts = (0:N-1)' * dt;
    tn = (tsrc(:) - tsrc(1)) / max(tsrc(end) - tsrc(1), eps) * T; [tn, iu] = unique(tn);
    G = max(resampleV(tn, Gsrc, iu, ts, 0), 0); Ta = resampleV(tn, Tasrc, iu, ts, 25); ws = max(resampleV(tn, wssrc, iu, ts, 1), 0);
    Tc = cellTemp(G, Ta, ws, S, gT); Pref = mppSD(G, Tc, S) * aS * aP; pl = makeUniformPlant(G, Tc, S, aS, aP);
end

function y = resampleV(tn, v, iu, ts, dflt)
    if isempty(v), y = dflt*ones(numel(ts), 1); return; end
    v = v(:); if isscalar(v), y = v*ones(numel(ts), 1); return; end
    v = fillmissing(v, 'linear', 'EndValues', 'nearest');
    if all(isnan(v)), y = dflt*ones(numel(ts), 1); return; end
    y = interp1(tn, v(iu), ts, 'linear', 'extrap');
end

function [hrs, G, Ta, ws] = repDayProfile(Sx, best)
    m = Sx.dayId == best;
    G = max(fillmissing(Sx.Y{1}(m), 'linear', 'EndValues', 'nearest'), 0);
    Ta = 25*ones(size(G)); if ~isempty(Sx.Y{2}), t_ = fillmissing(Sx.Y{2}(m), 'linear', 'EndValues', 'nearest'); if ~all(isnan(t_)), Ta = t_; end, end
    ws = ones(size(G));    if ~isempty(Sx.Y{4}), w_ = fillmissing(Sx.Y{4}(m), 'linear', 'EndValues', 'nearest'); if ~all(isnan(w_)), ws = max(w_, 0); end, end
    if Sx.hasTime, tt = Sx.t(m); hrs = hour(tt) + minute(tt)/60 + second(tt)/3600; else, hrs = (0:numel(G)-1)' * 0.25; end
    hrs = hrs(:);
end

function [Gt, Tat, wst] = todProfile(hrs, G, Ta, ws, tod)
    [hh, iu] = unique([0; hrs(:); 24]);
    g = [0; G(:); 0]; ta = [Ta(1); Ta(:); Ta(end)]; w = [ws(1); ws(:); ws(end)];
    Gt = interp1(hh, g(iu), tod, 'linear', 0); Tat = interp1(hh, ta(iu), tod, 'linear', 'extrap'); wst = interp1(hh, w(iu), tod, 'linear', 'extrap');
end

function v = batVoc(SOC, cfg)
    v = cfg.Ncell * (3.3 + 0.9*SOC);
end

function P = uavElecPower(tau, u)
    A = u.nRot * pi * (u.dia/2)^2; Th = u.mass * 9.81;
    Pid = Th^1.5 / sqrt(2*u.rho*A);
    P = Pid / (u.FoM*u.etaMotor*u.etaESC) * tau.^1.5;
end

function prm = decodeGene(x, mp)
    prm = [mp.lamMax*x(1), ...
           10^(log10(mp.Krange(1))   + x(2)*(log10(mp.Krange(2))   - log10(mp.Krange(1)))), ...
           10^(log10(mp.phiRange(1)) + x(3)*(log10(mp.phiRange(2)) - log10(mp.phiRange(1))))];
end

function [J, comp] = gaFitness(x, TP, cfg0, mp)
    % multi-profile fitness: mean [tracking error, normalised convergence, ripple] + worst-profile tracking error + ripple-limit penalty
    prm = decodeGene(x, mp); nP = numel(TP); C = zeros(nP, 3); pen = 0; comp = [NaN NaN NaN];
    for q = 1:nP
        cfg = cfg0; cfg.N = TP{q}.N;
        o = simMPPT(TP{q}.pl, 3, prm, cfg); P = o.V .* o.I;
        if any(~isfinite(P)), J = 1e3; return; end
        mt = mpptMetrics(o.t, P, TP{q}.Pref);
        C(q,:) = [mt.TrkErrPct/100, mt.Conv_ms/(1e3*o.t(end)), mt.RipplePct/100];
        pen = pen + mp.ripPen * max(0, mt.RipplePct/TP{q}.rippleBase - mp.rippleCap);
    end
    comp = mean(C, 1);
    J = mp.gaW * comp' + 0.5*max(C(:,1)) + pen/nP;
end

function o = runAlg(pl, a, cfg, smcOpt, mp)
    switch a
        case {1, 2}, o = simMPPT(pl, a, [], cfg);
        case 3,      o = simMPPT(pl, 3, smcOpt, cfg);
        case 4,      o = simMPPT(pl, 4, smcOpt, cfg);
        case 5,      o = simMPPT(pl, 3, mp.smcDefault, cfg);
    end
end

function out = simMPPT(pl, alg, prm, cfg)
    dt = cfg.dt; N = cfg.N;
    if alg >= 3, Ts = cfg.TsSMC; else, Ts = cfg.TsPO; end
    nc = max(round(Ts/dt), 1); Tsc = nc*dt; scanN = max(round(cfg.scanT/dt), 1); holdN = round(cfg.scanHold/dt);
    Vr = zeros(N,1); Ir = zeros(N,1); Dr = zeros(N,1); Lr = zeros(N,1); Vbr = zeros(N,1); Ibr = zeros(N,1); Sr = zeros(N,1);
    D = cfg.D0; SOC = cfg.SOC0; Ib = 0; Vb = cfg.Vbat; if cfg.useBat, Vb = batVoc(SOC, cfg); end
    Vp = Vb*(1 - D);
    if pl.mode == 1, Iw = pl.Iph(1); else, Iw = 0; end
    [Ip, Iw] = pvArr(Vp, 1, pl, Iw); iL = Ip;
    Pprev = NaN; Vprev = NaN; Iprev = NaN; dirD = -1; eprev = 0;
    scanning = (alg == 4); scanK0 = 1; Pbest = 0; Dbest = D; Pscan = 0; lastScan = -inf; nScan = 0;
    for k = 1:N
        [Ip, Iw] = pvArr(Vp, k, pl, Iw);
        Ibst = (1 - D)*iL;
        if cfg.useBat
            Pl = cfg.Pload(k); a = batVoc(SOC, cfg) + cfg.Rb*Ibst; dsc = a*a - 4*cfg.Rb*Pl; if dsc < 0, dsc = 0; end
            Vb = 0.5*(a + sqrt(dsc)); Ib = Pl/Vb - Ibst;
            SOC = min(max(SOC - cfg.socAccel*Ib*dt/(3600*cfg.Cb), 0), 1);
        end
        Vr(k) = Vp; Ir(k) = Ip; Dr(k) = D; Lr(k) = iL; Vbr(k) = Vb; Ibr(k) = Ib; Sr(k) = SOC;
        iLn = iL + dt*(Vp - cfg.rL*iL - (1 - D)*Vb)/cfg.L; if iLn < 0, iLn = 0; end
        Vp = Vp + dt*(Ip - iLn)/cfg.Cin; if Vp < 0, Vp = 0; end
        iL = iLn;
        P = Vr(k)*Ir(k);
        if scanning
            kk = k - scanK0;
            if P > Pbest, Pbest = P; Dbest = Dr(k); end
            if kk >= scanN
                scanning = false; D = Dbest; Pscan = Pbest; lastScan = k; nScan = nScan + 1; Vprev = NaN; eprev = 0;
            else
                D = cfg.Dmax - (cfg.Dmax - cfg.Dmin)*(kk + 1)/scanN;
            end
        elseif mod(k, nc) == 0
            V = Vr(k); I = Ir(k);
            if alg == 4 && Pscan > 0 && (k - lastScan) > holdN && (P < 0.7*Pscan || P > 1.3*Pscan)
                scanning = true; scanK0 = k; Pbest = 0; Dbest = D;
            else
                switch alg
                    case 1
                        if ~isnan(Pprev) && P < Pprev, dirD = -dirD; end
                        Pprev = P; D = D + dirD*cfg.stepD;
                    case 2
                        if isnan(Vprev)
                            D = D - cfg.stepD;
                        else
                            dV = V - Vprev; dI = I - Iprev;
                            if abs(dV) < 1e-3
                                if dI > 1e-3, D = D - cfg.stepD; elseif dI < -1e-3, D = D + cfg.stepD; end
                            else
                                g = dI/dV + I/max(V, 1e-3);
                                if g > cfg.incEps, D = D - cfg.stepD; elseif g < -cfg.incEps, D = D + cfg.stepD; end
                            end
                        end
                        Vprev = V; Iprev = I;
                    otherwise
                        if isnan(Vprev)
                            D = D - cfg.stepD; e = 0;
                        else
                            dV = V - Vprev; dI = I - Iprev;
                            if abs(dV) > 1e-3
                                gI = I/max(V, 1e-3);
                                e = (dI/dV + gI) / max(gI, 1e-3);      % normalised incremental-conductance error (independent of irradiance level)
                                e = max(min(e, 3), -3);                 % clamp against measurement spikes
                            elseif abs(dI) > 1e-3
                                e = 0.1*sign(dI);
                            else
                                e = eprev;
                            end
                            Ssw = e + prm(1)*(e - eprev)/Tsc;
                            D = D - prm(2)*Tsc*max(min(Ssw/prm(3), 1), -1);
                        end
                        eprev = e; Vprev = V; Iprev = I;
                end
            end
            D = min(max(D, cfg.Dmin), cfg.Dmax);
        end
    end
    out = struct('t', (0:N-1)'*dt, 'V', Vr, 'I', Ir, 'D', Dr, 'iL', Lr, 'Vbat', Vbr, 'Ib', Ibr, 'SOC', Sr, 'nScan', nScan);
end

function m = mpptMetrics(t, P, Pref)
    % single, consistent metric definitions used in EVERY experiment
    %   Eff       : energy ratio sum(P)/sum(Pref)
    %   Conv_ms   : first time P >= 95 % of Pref (Pref > 1 W) that is SUSTAINED for 20 ms
    %   RipplePct : RMS of P minus its 50 ms moving mean, relative to mean(Pref)
    t = t(:); P = P(:); Pref = Pref(:); dt = t(2) - t(1);
    m.Eff = 100*sum(P)/max(sum(Pref), eps);
    m.TrkErrPct = 100*mean(abs(Pref - P))/max(mean(Pref), eps);
    m.RMSE_W = sqrt(mean((Pref - P).^2));
    okc = double(P >= 0.95*Pref & Pref > 1); h = max(round(0.02/dt), 1);
    if h > 1, sm = movsum(okc, [0 h-1]); i95 = find(sm >= h, 1, 'first'); else, i95 = find(okc > 0, 1, 'first'); end
    if isempty(i95), m.Conv_ms = 1e3*t(end); else, m.Conv_ms = 1e3*t(i95); end
    w = max(round(0.05/dt), 3); hp = P - movmean(P, w);
    m.RipplePct = 100*sqrt(mean(hp.^2))/max(mean(Pref), eps);
end

function [tEndH, SOCh] = enduranceSim(Pbus, dtMin, mp, uav)
    n = numel(Pbus); dtH = dtMin/60; SOC = mp.SOC0; SOCh = NaN(n, 1); tEndH = (n-1)*dtH;
    Pto = uavElecPower(uav.tauTO, uav); Pcr = uavElecPower(uav.tauCruise, uav); Phv = uavElecPower(uav.tauHov, uav);
    for k = 1:n
        mn = mod((k-1)*dtMin, sum(mp.mission));
        if mn < mp.mission(1), Pl = Pto; elseif mn < mp.mission(1) + mp.mission(2), Pl = Pcr; else, Pl = Phv; end
        Voc = mp.Ncell*(3.3 + 0.9*SOC); Pnet = Pl - Pbus(k); dsc = Voc^2 - 4*mp.Rb*Pnet;
        if dsc < 0, SOCh(k) = SOC; tEndH = (k-1)*dtH; return; end
        Ib = (Voc - sqrt(dsc))/(2*mp.Rb);
        SOC = min(SOC - Ib*dtH/mp.Cb, 1); SOCh(k) = SOC;
        if SOC <= mp.SOCmin, tEndH = k*dtH; return; end
    end
end

function pk = findPeaksSimple(V, P)
    Ps = movmean(P, 25); idx = find(Ps(2:end-1) > Ps(1:end-2) & Ps(2:end-1) >= Ps(3:end)) + 1;
    idx = idx(Ps(idx) > 0.05*max(Ps)); [~, o] = sort(Ps(idx), 'descend'); idx = idx(o); pk = [];
    for q = 1:numel(idx)
        if isempty(pk) || all(abs(V(idx(q)) - V(pk)) > 0.08*max(V)), pk(end+1) = idx(q); end %#ok<AGROW>
    end
    pk = sort(pk);
end