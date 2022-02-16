NAME Problem2
ROWS
  N obj
  L c1
  L c2
  L c3
COLUMNS
  x1 obj 1 c1 -1
  x1 c2 1
  x2 obj 1 c1 1
  x2 c3 1

RHS
  rhs c1 -1 c2 0.3 c3 -0.5
BOUNDS
  FR bind x1
  FR bind x2
ENDATA