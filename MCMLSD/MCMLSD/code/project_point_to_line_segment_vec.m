function [q, d] = project_point_to_line_segment_vec(A,B,p)
  % returns q the closest point to p on the line segment from A to B 
%   q = zeros(size(p));
  % vector from A to B
  % Vectorize input
  [m, n] = size(p);
  A_v = repmat(A,m,1);
  B_v = repmat(B, m,1);
  AB = (B_v-A_v);
  % squared distance from A to B
  AB_squared = dot(AB(1,:),AB(1,:));
  if(AB_squared == 0)
    % A and B are the same point
    q = A_v;
  else
    % vector from A to p
    Ap = (p-A_v);
    t = dot(Ap,AB, 2)/AB_squared;
    q = A_v + repmat(t,1,2).* AB;
    q(t < 0.0,:) = A_v(t < 0.0,:);
    q(t > 1.0,:) = B_v(t > 1.0,:);
  end
  diff = q-p;
  d = sqrt(diff(:,1).^2 + diff(:,2).^2);
end