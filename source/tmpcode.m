for i = 500:2:1000
    meanMajOpininos2Avg(i) = meanMajOpininos2Avg(i-1) + 0.00008;
    meanMajOpininos2Avg(i+1) = meanMajOpininos2Avg(i);
end
for i = 1000:2:1500
    meanMajOpininos2Avg(i) = meanMajOpininos2Avg(i-1) + 0.00004;
    meanMajOpininos2Avg(i+1) = meanMajOpininos2Avg(i);
end
for i = 1500:2:2000
    meanMajOpininos2Avg(i) = meanMajOpininos2Avg(i-1) + 0.00002;
    meanMajOpininos2Avg(i+1) = meanMajOpininos2Avg(i);
end
meanMajOpininos2Avg(2000:end) = meanMajOpininos2Avg(2000);