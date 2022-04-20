% From the previous analysis, A good plot was found over the log plot
% of:

N = [1 10 100 1000 10000 100000 1000000]
a = [1120 1048 991 809 643 560 425]

a1 = -.0712
b1 = 3.0935

%x = 10^b*(10.^(a1.*log(N)))
(10^b)*(10.^(a1*log(N)))
