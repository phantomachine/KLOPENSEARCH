function [choice,exitflag] = choicekilljoy(question_str,title_str)

%% CHOICEKILLJOY.M
%% 
%% ========================================================================
%% PURPOSE: To produce GUI with Kill or Continue buttons
%%
%% INPUT: question_str : 'Question String'
%%        title_str    : 'Title of Dialog Box String'
%% 
%% OUTPUT: choice      :  string values: 'Kill' or 'Skip' or 'Continue'
%%         exitflag    : 1 or 2 or 3
%%
%% Usage Example:
%% --------------       
%%    question_str = 'How now brown cow?';
%%    title_str = 'Title of dialog box';
%%    [choice,killjoyflag] = choicekilljoy(question_str,title_str);
%%
%%    if strcmp(choice,'Kill') % killjoyflag == 2
%%        return
%%    elseif strcmp(choice,'Continue') % killjoyflag == 1
%%        [output_do] = dosomething(input_do);
%%    end
%%
%% See also QUESTDLG (c) 2009 T.Kam Email: mortheus@gmail.com
%% ========================================================================

choice = questdlg(question_str, ...
                    title_str, ...
                    'Kill','Skip','Continue','Continue'); 
                                       % default button highlight: Continue

% Handle response
switch choice
    case 'Continue'
        disp([choice, ' ', title_str])
        exitflag = 1;
    case 'Kill'
        disp([choice, ' ', title_str])
        exitflag = 2;
    case 'Skip'
        disp([choice, ' ', title_str])
        exitflag = 3;
end
