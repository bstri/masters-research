function[HollowSym] = hollowSymmetricCreator(n)
    A = randi(10, n);
    L = tril(A);
    HollowSym = L + L' - 2*(A.*eye(n));
end
% After running this function, a good way to copy this
% matrix over to a text file is to run the command
% num2str(HollowSym,'%3d')
% and then copy the result from the command window
% OR
% you can right click on the HollowSym variable in the Workspace,
% and select Open Selection
% From there you can ctrl+a and copy