function assert_h(x)
%assert that *x* is a horizontal vector (its vertical size <= 1) 
    assert(size(x,1)<=1)
end
