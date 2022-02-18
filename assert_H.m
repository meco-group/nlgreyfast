function assert_H(x,rows)
%assert that *x* is a horizontally long and vertically short matrix, and the vertical size equals *rows*.
    assert(size(x,2)==rows)
end
