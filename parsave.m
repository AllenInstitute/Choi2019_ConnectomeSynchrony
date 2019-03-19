function parsave(fname, phi,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   Save phi and t during the parallel integration of the differential eqs.

%   Code used for simulations in Choi & Mihalas (2019)
%   Written by Hannah Choi, 2019 (hannahch@uw.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(fname,'phi','t')
end