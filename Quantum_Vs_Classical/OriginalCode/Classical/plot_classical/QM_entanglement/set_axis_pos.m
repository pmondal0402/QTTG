function[] = set_axis_pos(ax, left, top, right, bot)
pos = get( ax, 'Position' )

% Move vertically
if left > 0
  pos(1) = left ; 
end

if top > 0
  pos(2) = top ; 
end

if right > 0
  pos(3) = right ; 
end

if bot > 0
  pos(4) = bot ; 
end

if top > 0 || left > 0 || right > 0 || bot > 0
  set(ax, 'Position', pos)
end

