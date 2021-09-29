function y=iif(condition,a,b)
    % iif: one-line IF for MATLAB. 
    % In case `condition` is true, returns `a`.
    % In case `condition` is false, returns `b`.
    % Note that if `a` or `b` is an expression, both of them
    % are evaluated before running into `iif`, not only the one returned
    % based on the condition. This can have unexpected side effects and/or 
    % performance effects, as it differs from most implementations in other 
    % programming languages. As a result, this function is a substitute that
    % you can use to make some code lines shorter, but take into
    % consideration the performance effects and use `if` instead wherever 
    % possible.
    if condition, y=a; else, y=b; end
