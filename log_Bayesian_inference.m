D=importdata('WBC_data.xlsx');
n = size(D.textdata, 1);
patientID = D.textdata(2:n,1);
time = D.data(:,3);
wbc = D.data(:,4);
n = n - 1;

T = importdata('Dx_to_Tx_data.xlsx');

nTx = size(T.textdata, 1);
Tx_patientID = T.textdata(2:nTx,1);
Tx = T.data;
nTx = nTx - 1;

npatients = 1;
for i = 2:n
    if ~strcmp(patientID(i), patientID(i-1))
        npatients = npatients + 1;
    end
end

S_patientID = Tx_patientID;

start_ind = zeros(npatients, 1);
end_ind = zeros(npatients, 1);
start_ind(1) = 1;
end_ind(1) = 1;
ipatient = 1;
if (~strcmp(patientID(1), Tx_patientID(1))) | (~strcmp(patientID(1), S_patientID(1)))
    error(['IDs do not match at ' Tx_patientID(1)]);
end
for i = 2:n
    if ~strcmp(patientID(i), patientID(i-1))
        ipatient = ipatient + 1;
        if (~strcmp(patientID(1), Tx_patientID(1))) | (~strcmp(patientID(1), S_patientID(1)))
            error(['IDs do not match at ' Tx_patientID(ipatient)]);
        end
        start_ind(ipatient) = i;
        end_ind(ipatient) = i;    
    else
        if isnan(Tx(ipatient)) || time(i) < Tx(ipatient)
            end_ind(ipatient) = i;
        end
    end
end

Results = cell(npatients, 10);

n_samples = 10000;
n_skip = 2000;
n_delta = 10;
k_bar = log(10^6);
X0_low = 1;
X0_bar = 300;
r_low = 0;
r_bar = 5;
kappa = 72;
Psi = kappa * 0.1266^2;

ig_alpha = kappa / 2;
ig_beta = Psi / 2;
inv_gamma_prior = @(x) ig_beta^ig_alpha / gamma(ig_alpha) * x.^(-ig_alpha-1) .* exp(-ig_beta./x);

n_hist_bins = 50;

nrows = 4;
ncols = 4;
npage = nrows;

nfigs = ceil(npatients/npage);
param_labels = ["log_{10}(K)", "X_0", "r", "\sigma^2"];

all_samples = cell(1, npatients);

for fig_ind = 1:nfigs
    f = figure;

    if fig_ind < nfigs
        nsubs = npage;
    else
        nsubs = npatients - npage*(nfigs-1);
    end
    
    for isub = 1:nsubs
        i = npage*(fig_ind-1)+isub;
        pid = patientID{start_ind(i)};
        
        Results{i, 1} = pid;
        Results{i, 2} = "excl";
        for ir = 3:10
            Results{i, ir} = "";
        end
        
        if end_ind(i) - start_ind(i) < 3
            continue
        end

        times = time(start_ind(i):end_ind(i));
        values = wbc(start_ind(i):end_ind(i));

        i, length(times)

        X0_init = values(1);
        k_init = 1 * values(length(values));
        r_init = 0.5;
        sig2_init = 0.1^2;

        samples = logistic_model_gibbs_sampler(n_samples, n_skip, n_delta, times, values, k_bar, X0_low, X0_bar, r_low, r_bar, kappa, Psi, k_init, X0_init, r_init, sig2_init);
        
        all_samples{i} = samples;
        
        for k_sub = 1:4
            k = k_sub;
            if k == 2 || k == 3
                k = 5 - k;
            end
            h = subplot(nrows, ncols, (isub-1)*nsubs+k_sub);
            histogram(samples(:,k), n_hist_bins, 'Normalization', 'pdf');
            if k == 1
                h = histogram(log10(exp(1))*samples(:,k), n_hist_bins, 'Normalization', 'pdf');
                ylabel("P" + num2str(i));
                hold on;
                line([0 6], [1./6 1./6], 'Color', 'g', 'LineWidth', 1);
                xlim([0 6]);
                probk = sum(samples(:,k) < log(10^3)) / n_samples;
                text(0.2, h.Parent.YLim(1) + 0.8 * (h.Parent.YLim(2) - h.Parent.YLim(1)), num2str(probk), 'Color', 'r', 'FontSize', 8);
                Results{i, 2} = num2str(probk);
            elseif k == 2
                hold on;
                line([0 300], [1./300 1./300], 'Color', 'g', 'LineWidth', 1);
                xlim([0 150]);
            elseif k == 3
                hold on;
                line([0 5], [1./5 1./5], 'Color', 'g', 'LineWidth', 1);
                xlim([0 2]);
            else
                hold on;
                fplot(inv_gamma_prior, [0 0.12], 'g');
                xlim([0 0.12]);
            end
            xlabel(param_labels(k));
        end

        medlowhigh = zeros(1, 12);
        for k = 1:4
            if k == 1
                [med, low, high] = median_with_errors(exp(samples(:, k)), 0.15, 0.85);
            else
                [med, low, high] = median_with_errors(samples(:, k), 0.15, 0.85);
            end
            Results{i, 2*k+1} = num2str(med);
            Results{i, 2*k+2} = "(" + num2str(low) + ", " + num2str(high) + ")";
        end

        i
    end
    
    saveas(f, ['posterior_histograms_' int2str(fig_ind) '.pdf'], 'pdf')
end


formatSpec = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
fileID = fopen('results.txt','w');
for i = 1:npatients
    fprintf(fileID,formatSpec,Results{i,:});
end
fclose(fileID);

save('all_samples.mat', 'all_samples');
    

