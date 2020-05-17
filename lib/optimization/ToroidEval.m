function [] = ToroidEval(x,D,fn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to evaluate the toroidal core inductor metrics   
%                                                                         %
% Written by Vinicius Nascimento                                          %
% viniciuscn1@gmail.com
%                                                                         %
% Current revision date: 08-25-2019                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% if fn not parsed, plotting is disabled 
if(nargin<3)
    fn = 0;
end

% call fitness function in evaluation mode
f = ToroidFit(x,D,fn);

% display constraints that were not satisfied
if any(f.c<1)                   % if any constraint was not satisfied
    fprintf('\nList of constraints not satisfied: --------------------\n');
    for k = 1:length(f.c)       % loop through each constraint
        if f.c(k)<1             % if k-th constraint was not satisfied
            PrintConstraint(k); % print non satisfied constraint
        end
    end
end
end

function [] = PrintConstraint(c)
    switch c
        case 1
            fprintf('1. Wire length constraint was not satisfied.\n');
        case 2
            fprintf('2. Maximum packing factor constraint was not satisfied.\n');
        case 3
            fprintf('3. Maximum aspect ratio allowed was exceeded.\n');
        case 4
            fprintf('4. Maximum mass allowed was exceeded.\n');
        case 5
            fprintf('5. Minimum incremental inductance was not met.\n');
        case 6
            fprintf('6. Maximum allowed flux density was exceeded.\n');
        case 7
            fprintf('7. TEC did not converge.\n');
        case 8
            fprintf('8. Maximum allowed losses was exceeded.\n');
        case 9
            fprintf('9. Maximum allowed winding temperature was exceeded.\n');
        case 10
            fprintf('10. Maximum allowed core temperature was exceeded.\n');
        otherwise
            fprintf('A non-listed constraint was not satisfied.\n');
    end
end