% prototype attention-based anomaly discovery

clear all; close all;

%%

train_file = '/Users/DoerLBH/Dropbox/git/AttentiveCancerExplorer/data/selected_loc_train.txt';
test_file = '/Users/DoerLBH/Dropbox/git/AttentiveCancerExplorer/data/selected_loc_test.txt';

train_locs = importdata(train_file);
test_locs = importdata(test_file);
read_length = 100;
BRAF_START = 140719327;
BRAF_END = 140924764;

%%

train_cover = zeros(length(train_locs),read_length);

parfor t = 1:length(train_locs)
    train_cover(t,:) = train_locs(t):train_locs(t)+read_length-1;
end

test_cover = zeros(length(test_locs),read_length);

parfor t = 1:length(test_locs)
    test_cover(t,:) = test_locs(t):test_locs(t)+read_length-1;
end

%%  

train_flatten = reshape(train_cover,[1,length(train_locs)*read_length]);
test_flatten = reshape(test_cover,[1,length(test_locs)*read_length]);

%% 

figure(1)
h1 = histogram(train_flatten - BRAF_START); hold on
h2 = histogram(test_flatten - BRAF_START);
xlim([0 BRAF_END-BRAF_START+1])
legend({'training','testing'});
title('read coverage in BRAF range')
