function paperfig()
%PAPERFIG draws the figures of the paper [1].
% WARNING : due to extremly fine discretization in some figures, the
% execution of this function is VERY SLOW... You had better to execute only
% a part of this. You can also change parameters to get similar figures
% with lower resolution in a more reasonable time.
%
% The functions called are all located in the folder 'figfun'


% Fig. 2 : Evaluation of the procedure which determine the differentiation order.
fig_diff();

% Fig 3 : the test-signal and its component
[s s1 s2 s3] = gen_tests(1);

% Fig 4 : Sensitivity to the spline order
sensitiv_spo();

% Fig 5 : Sensitivity to the initial guess
sensitiv_init();

% Fig 6 : Orthogonality index
sensitiv_io();

% Fig 7 : separation plots in the (a,f)-plane
separTest();  %%% WARNING NOT WORKING!!

% Fig 8 : 



end

