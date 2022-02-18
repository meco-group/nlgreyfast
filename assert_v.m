function assert_v(x)
%assert that *x* is a vertical vector (its horizontal size <= 1) 
    assert(size(x,2)<=1)
end
