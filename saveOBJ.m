function [] = saveOBJ(frameIdx)

global Frame labels V0i;

ColorList = generateColorList();
for i=1:size(V0i{frameIdx},1)
    Color(i,:) = ColorList( labels(i), : );
end
succ2 = write_to_obj(V0i{frameIdx}, Color, 'refinecluster.obj');


end