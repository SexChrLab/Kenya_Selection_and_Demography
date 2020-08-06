# This script computes the threshold for DP (mean DP, 50% of mean, and 150% of mean)
data = read.table("/scratch/tphung3/Kenya_Selection_and_Demography/Population_History/02_find_callable_regions/DP/autosomes_DP_sum.txt")

mean_DP = sum(data[,1])/sum(data[,2])

lower_threshold = mean_DP*0.5
higher_threshold = mean_DP*1.5

print(mean_DP)
print(lower_threshold)
print(higher_threshold)