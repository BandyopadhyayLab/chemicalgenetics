% Function to take a 384 well plate and extract a matrix of 96 rows and 4
% columns correspondign to 4x replicates.
%

function [rowvector]=extract96from384(sdata)
    % every other column
    a = sdata(:,[1:2:24]);

    q1 = a([1:2:16],:);
    q3 =  a([2:2:16],:);


    a =sdata(:,[2:2:24]);
    q2 = a([1:2:16],:);
    q4 =  a([2:2:16],:);

    a=reshape(q1',1,96)';
    b=reshape(q2',1,96)';
    c=reshape(q3',1,96)';
    d=reshape(q4',1,96)';


    rowvector = [a b c d];

   
end




