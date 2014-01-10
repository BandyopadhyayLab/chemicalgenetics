function newarray = extract384from96(vect)


    a=reshape([1:384],16,24);
    newmap = extract96from384(a);


    newarray = ones(16,24)-1;
    for i=[1:384]
        newarray(newmap(i)) = vect(i);
    end
end
