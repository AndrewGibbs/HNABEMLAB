function meshEl = wholeMeshSide(L)
%creates single mesh element on an entire side
    meshEl.distL = 0;
    meshEl.distR = 0;
    meshEl.distSideL = 0;
    meshEl.distSideR = 0;
    meshEl.width = L;
end

