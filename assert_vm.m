function assert_vm(x,cols)  
%assert that *x* is a vertically tall and horizontally skinny matrix, and the horizontal size equals *cols*.
    assert(size(x,2)==cols)
end
