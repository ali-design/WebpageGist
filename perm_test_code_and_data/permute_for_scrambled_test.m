% This script loads the data from the web page gist paper, and runs
% permutation tests on the data from all experiments. Set "nrep" to however
% many simulations desired. This code takes a while to run, decrease "nrep"
% to make it go faster (but with lower statistical power).

% Search the document for "IN THE PAPER" to find out where the reported
% values are finally computed.

% Shaiyan Keshvari, 8/23/2017, shaiyan@mit.edu

% load gist data and scrambled data
load results_normalText_334subjects
responseStruct_gist = responseStruct;

load results_scrambled_334subjects
responseStruct_scrambled = responseStruct;

% intialize
gist_data = [];
scrambled_data = [];
gist_mat = [];
scrambled_mat = [];

% for each subject in each experiment
for jj = 1:10
    
    % get the subject's data and save in gist_data
    responseStruct_curr = responseStruct_gist(((jj-1)*334 + 1):(jj*334),:);
    gist_data{jj} = responseStruct_curr;
    
    % for each trial of the current subject, save it in gist_data again
    for kk = 1:size(responseStruct_curr,1)
        gist_data{jj}{kk,1} = str2num(responseStruct_curr{kk,1}(4:end));
    end
    
    % sort the gist data
    gist_data{jj} = sortrows(gist_data{jj},1);
    
    % add the data into a matrix. Each row is a subject's responses
    gist_mat = [gist_mat cell2mat(gist_data{jj}(:,4))];
    
    % get scrambled letter data
    responseStruct_curr = responseStruct_scrambled(((jj-1)*334 + 1):(jj*334),:);
    scrambled_data{jj} = responseStruct_curr;
    for kk = 1:size(responseStruct_curr,1)
        scrambled_data{jj}{kk,1} = str2num(responseStruct_curr{kk,1}(4:end));
    end
    scrambled_data{jj} = sortrows(scrambled_data{jj},1);
    scrambled_mat = [scrambled_mat cell2mat(scrambled_data{jj}(:,4))];
        
end

% more intialization
gist_unique = zeros(size(scrambled_mat,1),1);
scrambled_unique = gist_unique;
rand_unique = gist_unique;
rand_mat = zeros(size(scrambled_mat));

% generate a column-wise random permutation of gist_mat
for ii = 1:size(rand_mat,2)
    rand_mat(:,ii) = gist_mat(randperm(size(gist_mat,1)),ii);
end

% find number of unique values in each row of gist_mat
for kk = 1:length(gist_unique)
    gist_unique(kk)=length(unique(gist_mat(kk,:)));
    scrambled_unique(kk)=length(unique(scrambled_mat(kk,:)));
    rand_unique(kk)=length(unique(rand_mat(kk,:)));
end

% get histogram of unique values
gist_hist = hist(gist_unique,1:10);
scrambled_hist = hist(scrambled_unique,1:10);
rand_hist = hist(rand_unique,1:10);

% show the histogram
% bar(1:10,[gist_hist ;scrambled_hist ;rand_hist]')
% legend('gist','scrambled','gist rand permuted')

%% Let's get some p-values

% holds pc (percent correct) per category
pc_per_category_gist = zeros(1,10);
pc_per_category_scrambled = zeros(1,10);

% holds the correct answer for each subject for each trial
correct_resp_gist = cell2mat(gist_data{1}(:,3));
correct_resp_scrambled = cell2mat(scrambled_data{1}(:,3));
correct_resp_gist_repmat = repmat(correct_resp_gist,1,10);
correct_resp_scrambled_repmat = repmat(correct_resp_scrambled,1,10);
pc_per_category_gist_stdev = zeros(10,10);
pc_per_category_scrambled_stdev = zeros(10,10);

% for each category
for ii = 1:10

    % get correct responses per category
    pc_per_category_gist(ii) = sum((gist_mat(:)==correct_resp_gist_repmat(:)) & (correct_resp_gist_repmat(:)==ii))/sum((correct_resp_gist_repmat(:)==ii));
    pc_per_category_scrambled(ii) = sum((scrambled_mat(:)==correct_resp_scrambled_repmat(:)) & (correct_resp_scrambled_repmat(:)==ii))/sum((correct_resp_scrambled_repmat(:)==ii));
    
    % get pc for each subject and category
    % for each subject
    % COLUMNS ARE SUBJECTS
    for jj = 1:10
        curr_subj = zeros(size(gist_mat));
        curr_subj(:,jj) = 1;
        curr_subj = curr_subj == 1;
        pc_per_category_gist_stdev(ii,jj) = sum((gist_mat(curr_subj)==correct_resp_gist_repmat(curr_subj)) & (correct_resp_gist_repmat(curr_subj)==ii))/sum((correct_resp_gist_repmat(curr_subj)==ii));
        pc_per_category_scrambled_stdev(ii,jj) = sum((scrambled_mat(curr_subj)==correct_resp_scrambled_repmat(curr_subj)) & (correct_resp_scrambled_repmat(curr_subj)==ii))/sum((correct_resp_scrambled_repmat(curr_subj)==ii));
    end
end

% get correct responses per subject
pc_per_subj_gist = sum(gist_mat==correct_resp_gist_repmat,1)./size(gist_mat,1);
pc_per_subj_scrambled = sum(scrambled_mat==correct_resp_scrambled_repmat,1)./size(scrambled_mat,1);

%% RUN THE SIMULATION
% number of permutations
num_reps = 5000;

% holds each simulated experiment
gist_mat_randperm = zeros(size(gist_mat));
scrambled_mat_randperm = gist_mat_randperm;

% hold percent correct per simulation
randperm_pc_per_category_gist = zeros(10,1);
randperm_pc_per_category_scrambled = zeros(10,1);
randperm_pc_per_subj_gist = zeros(10,1);
randperm_pc_per_subj_scrambled = zeros(10,1);

% hold percent correct over all simulations
all_gist_rand_pc_per_category = zeros(num_reps,10);
all_scrambled_rand_pc_per_category = zeros(num_reps,10);
all_gist_rand_pc_per_subj = zeros(num_reps,10);
all_scrambled_rand_pc_per_subj= zeros(num_reps,10);

% for each rep
for kk = 1:num_reps
    
    % randomly permute data
    % for each subject (columns) 
    for nn = 1:size(gist_mat,2)
        % permute the gist data and the scrambled data
        gist_mat_randperm(:,nn) = gist_mat(randperm(size(gist_mat,1)),nn);
        scrambled_mat_randperm(:,nn) = scrambled_mat(randperm(size(gist_mat,1)),nn);
    end
    
    % get correct responses. rows are trials, columns are subjects
    gist_correct_tmp = correct_resp_gist_repmat == gist_mat_randperm;
    scrambled_correct_tmp = correct_resp_scrambled_repmat == scrambled_mat_randperm;
    
    % for each category
    for ii = 1:10
        
        % pc for each category, averaged over subjects
        randperm_pc_per_category_gist(ii) = sum(gist_correct_tmp(:) & (correct_resp_gist_repmat(:)==ii))/sum(correct_resp_gist_repmat(:)==ii);
        randperm_pc_per_category_scrambled(ii) = sum(scrambled_correct_tmp(:) & (correct_resp_scrambled_repmat(:)==ii))/sum(correct_resp_scrambled_repmat(:)==ii);
        
    end
    
    % pc for each subject, averaged over categories
    randperm_pc_per_subj_gist = (sum(gist_correct_tmp,1)./size(gist_correct_tmp,1))';
    randperm_pc_per_subj_scrambled = (sum(scrambled_correct_tmp,1)./size(scrambled_correct_tmp,1))';

    
    % save the results
    all_gist_rand_pc_per_category(kk,:) = randperm_pc_per_category_gist;
    all_scrambled_rand_pc_per_category(kk,:) = randperm_pc_per_category_scrambled;
    all_gist_rand_pc_per_subj(kk,:) = randperm_pc_per_subj_gist;
    all_scrambled_rand_pc_per_subj(kk,:) = randperm_pc_per_subj_scrambled;
    
end

% get p-values per category
gist_p_value_category = squeeze(mean(bsxfun(@ge,all_gist_rand_pc_per_category,pc_per_category_gist),1));
scrambled_p_value_category = squeeze(mean(bsxfun(@ge,all_scrambled_rand_pc_per_category,pc_per_category_scrambled),1));

% get p-values per subject
gist_p_value_subj = squeeze(mean(bsxfun(@ge,all_gist_rand_pc_per_subj,pc_per_subj_gist),1));
scrambled_p_value_subj = squeeze(mean(bsxfun(@ge,all_scrambled_rand_pc_per_subj,pc_per_subj_scrambled),1));

% Bonferroni correction - IN THE PAPER
gist_p_value_category_BF_corrected(gist_p_value_category==0) = 2*10*1/num_reps;
scrambled_p_value_category_BF_corrected(scrambled_p_value_category==0) = 2*10*1/num_reps;

%% now do stuff for the ads. Columns are 1-category, 2-response, 3-correct
load('all_response_exp_adNoAd-Menu')
resp_mat = zeros(200,10);
correct_cat_mat = zeros(200,10);  
is_correct_mat = zeros(200,10);  

% initialize data (it was out of order in the original data, for the ads)
for jj = 1:10
    
    tmp =  responseStruct(((jj-1)*200 + 1):(jj*200),:);
    tmp = sortrows(tmp,1);
    resp_mat(:,jj) = cell2mat(tmp(:,4)); % subject resp
    correct_cat_mat(:,jj) = cell2mat(tmp(:,3)); % correct category
    is_correct_mat(:,jj) = cell2mat(tmp(:,5))==1; % is correct resp
    
end
is_correct_mat_fixed=is_correct_mat;
is_correct_mat_fixed(1:50,:) = 1- is_correct_mat(1:50,:);
is_correct_mat_fixed(151:200,:) = 1- is_correct_mat(151:200,:);
ad_correct = [is_correct_mat_fixed(1:50,:); is_correct_mat_fixed(151:200,:)];
menu_correct = is_correct_mat_fixed(51:150,:);
ad_ground_truth = [ones(50,10); 2*ones(50,10)];
menu_ground_truth = ad_ground_truth;


% 1 if responded present, 2 if responded absent
ad_responses = ad_correct;
tmp = ad_correct(1:51,:); % ad present trials
tmp(ad_correct(1:51,:)==0)=2; % mark miss trials 
ad_responses(1:51,:) = tmp;
tmp = ad_correct(51:100,:); % ad absent trials
tmp(ad_correct(51:100,:)==1)=2; % mark correct reject
tmp(ad_correct(51:100,:)==0)=1; % mark false alarm
ad_responses(51:100,:) = tmp;

% menu
menu_responses = menu_correct;
tmp = menu_correct(1:51,:); % menu side
tmp(menu_correct(1:51,:)==0)=2; % mark top responses
menu_responses(1:51,:) = tmp;
tmp = menu_correct(51:100,:); % menu top trials
tmp(menu_correct(51:100,:)==1)=2; % mark top responses
tmp(menu_correct(51:100,:)==0)=1; % mark left responses
menu_responses(51:100,:) = tmp;

% run simulation
nrep = 500000;
[~,ad_perms] = sort(rand([size(ad_responses) nrep]));
[~,menu_perms] = sort(rand([size(menu_responses) nrep]));

ad_randomized_tmp = zeros(size(ad_perms));
menu_randomized_tmp = zeros(size(menu_perms));
ad_randomized = zeros([size(ad_perms,1) size(ad_perms,2) 10*nrep]);
menu_randomized = ad_randomized;

for kk = 1:10
    for ii = 1:size(menu_responses,1)
        for jj = 1:size(menu_responses,2)
            ad_randomized_tmp(ii,jj,:) = ad_responses(ad_perms(ii,jj,:),jj);
            menu_randomized_tmp(ii,jj,:) = menu_responses(menu_perms(ii,jj,:),jj);
        end
    end
    ad_randomized(:,:,(nrep*(kk-1)+1):(nrep*kk)) = ad_randomized_tmp;
    menu_randomized(:,:,(nrep*(kk-1)+1):(nrep*kk)) = menu_randomized_tmp;
end

% compute performance, real and simulated
perms_correct_ad = bsxfun(@eq,ad_ground_truth,ad_randomized);
overall_correct_ad = squeeze(sum(sum(perms_correct_ad,2),1)/numel(ad_ground_truth));

perms_correct_menu = bsxfun(@eq,menu_ground_truth,menu_randomized);
overall_correct_menu = squeeze(sum(sum(perms_correct_menu,2),1)/numel(menu_ground_truth));

real_perf_ad = sum(sum(ad_ground_truth==ad_responses,1),2)/numel(ad_ground_truth);
real_perf_menu = sum(sum(menu_ground_truth==menu_responses,1),2)/numel(menu_ground_truth);
real_diff = real_perf_ad - real_perf_menu;

% get p-values for above chance performance - IN THE PAPER
p_ad_v_chance =  2*mean(overall_correct_ad > real_perf_ad);
p_ad_v_chance(p_ad_v_chance==0) = 2/nrep;
p_menu_v_chance = 2*mean(overall_correct_menu > real_perf_menu);
p_menu_v_chance(p_menu_v_chance==0) = 2/nrep;


% get p-values for whether ad is better than menu - IN THE PAPER
p_ad_greater = 2*mean((overall_correct_ad-overall_correct_menu)<=real_diff);
p_ad_greater(p_ad_greater==0) = 2/nrep;

% do it per subject
overall_correct_ad_per_subj = squeeze(sum(perms_correct_ad,1)/size(ad_ground_truth,1));
overall_correct_menu_per_subj = squeeze(sum(perms_correct_menu,1)/size(menu_ground_truth,1));
real_diff_per_subj = sum(ad_ground_truth==ad_responses,1)/size(ad_ground_truth,1)...
    -sum(menu_ground_truth==menu_responses,1)/size(menu_ground_truth,1);

% Bonferroni correction
p_ad_greater_per_subj = 2*10*mean(bsxfun(@le,overall_correct_ad_per_subj-overall_correct_menu_per_subj,real_diff_per_subj'),2);

correct_mat= bsxfun(@eq,B,A);
C = mean(correct_mat,1);
% reshape responses 
num_reps = 5000;
sim_diff_per_cat = zeros(num_reps,10);

% get real differences
for ii = 1:10
    real_diff_per_cat(ii) =  mean(gist_mat(correct_resp_gist_repmat==ii)==ii)-mean(scrambled_mat(correct_resp_scrambled_repmat==ii)==ii);
end

%% Get p-value for intact vs scrambled over categories ands subjects

% for each rep, do the scrambled vs intact for each category
for kk = 1:num_reps
    
    % randomly permute data
    % for each subject (columns)
    
    for ii = 1:10 %over categories
        tmp_gist = gist_mat(correct_resp_gist_repmat==ii);
        tmp_scramble = scrambled_mat(correct_resp_gist_repmat==ii);
        
        tmp_all = [tmp_gist;tmp_scramble];
        tmp_all(:) = tmp_all(randperm(length(tmp_all)));
        
        gist_correct = mean(tmp_all(1:(length(tmp_all)/2))==ii);
        scrambled_correct = mean(tmp_all((length(tmp_all)/2+1):end)==ii);
        
        sim_diff_per_cat(kk,ii) = gist_correct - scrambled_correct;
        
    end    

end

% compute how many were greater and bonferroni correct
p_value_diff_correct = 2*10*mean(bsxfun(@ge,sim_diff_per_cat,real_diff_per_cat),1);

% do it over categories
% real difference between experiments
num_reps = 5000;
real_diff_overall = mean(gist_mat(:)==correct_resp_gist_repmat(:)) - mean(scrambled_mat(:)==correct_resp_scrambled_repmat(:));

% simulated
sim_diff_overall = zeros(num_reps,1);
tmp_all = [gist_mat(:);scrambled_mat(:)] == [correct_resp_gist_repmat(:);correct_resp_scrambled_repmat(:)];

for kk = 1:num_reps
    
        tmp_all(:) = tmp_all(randperm(length(tmp_all)));
        
        gist_correct = mean(tmp_all(1:(length(tmp_all)/2)));
        scrambled_correct = mean(tmp_all((length(tmp_all)/2+1):end));
        
        sim_diff_overall(kk) = gist_correct - scrambled_correct;
end

% Bonferonni-corrected value - IN THE PAPER
p_value_diff_overall = 2*mean(sim_diff_overall>real_diff_overall);
p_value_diff_overall(p_value_diff_overall==0) = 2/nrep;