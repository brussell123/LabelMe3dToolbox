function mesh = getSceneMesh(annotation)
% Computes triangular mesh given 3D scene.
%
% Inputs:
% annotation
%
% Outputs:
% mesh.vertices - 3xN 3D vertices (single)
% mesh.faces - 3xM triangle faces (int32)
% mesh.objndx - 1xM object indices for each face (int32)

% Parameters:
% TO DO: this should eventually be used to bound 3D scene:
vBB = [-10^5 10^5 -10^5 10^5 -100 10^8]; % volume bounding box (remove objects outside this volume)

mesh.vertices = [];
mesh.faces = [];
mesh.objndx = [];
for i = 1:length(annotation.object)
  if isfield(annotation.object(i).world3d,'polygon3d')
    switch annotation.object(i).world3d.type
     case {'standingplanes','groundplane'}
      % Get set of planes and 3D points:
      PI = getPlanes(annotation.object(i).world3d.plane);
      X = getLMpolygon3D(annotation.object(i).world3d.polygon3d);
      faces = [];
      for j = 1:size(PI,2)
        % Get 3D polygon points that live on plane:
        if strcmp(annotation.object(i).world3d.type,'standingplanes')
          n = zeros(1,length(annotation.object(i).world3d.polygon3d.pt));
          for k = 1:length(annotation.object(i).world3d.polygon3d.pt)
            for ii = 1:1:length(annotation.object(i).world3d.polygon3d.pt(k).planeindex)
              if str2num(annotation.object(i).world3d.polygon3d.pt(k).planeindex(ii).index)==(j-1)
                n(k) = 1;
              end
            end
          end
          n = find(n==1);
        else
          n = [1:length(annotation.object(i).world3d.polygon3d.pt)];
        end

        pp = PI(:,j);
        R = eye(4);
        if abs(pp(1)) < 1e-4
          % Rotate axes to avoid numerical instabilities:
          R = [0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
          pp = R*pp;
          if abs(pp(1)) < 1e-4
            pp = R*pp;
            R = R*R;
          end
        end

        % Get 2D to 3D plane transformation:
        T = [-pp(2)/pp(1) -pp(3)/pp(1) -pp(4)/pp(1); eye(3)];
        
        % Frontal-rectify plane points:
        Xj = X(:,n);
        xx = (T'*T)\T'*R*[Xj; ones(1,size(Xj,2))];
        yy = xx(2,:)./xx(3,:);
        xx = xx(1,:)./xx(3,:);

% $$$         keyboard;
        
        % Get Delaunay triangulation of plane polygon points:
        tri = delaunay(xx,yy,{'Qbb','Qc','QJ','Pp'})';

        % Determine which triangles are inside polygon:
        cx = mean(reshape(xx(tri),size(tri)),1);
        cy = mean(reshape(yy(tri),size(tri)),1);
        nInside = inpolygon(cx,cy,xx,yy);
        
% $$$         if size(PI,2)>1
% $$$           keyboard;
% $$$         end

        doDisplay = 0;
        if doDisplay%i==6%doDisplay
          figure;
          triplot(tri',xx,yy);
          hold on;
          plot([xx xx(1)],[yy yy(1)],'r');
          plot(xx,yy,'r+');
          plot(cx,cy,'g.');
        end
        
        % Remove faces that are outside of polygon:
        faces = [faces reshape(int32(n(tri(:,nInside))),3,sum(nInside))];
      end
      
      % Subdivide mesh:
      % pass single(X), faces
      vertices = single(X);
      
      mesh.faces = [mesh.faces faces+size(mesh.vertices,2)];
      mesh.vertices = [mesh.vertices vertices];
      mesh.objndx = [mesh.objndx int32(i*ones(1,size(faces,2)))];
    end
  end
end

return;

% $$$ % Parameters:
% $$$ sampDist = 50; % Sample distances (in centimeters)
% $$$ 
% $$$ mesh.vertices = [];
% $$$ mesh.faces = [];
% $$$ mesh.objndx = [];
% $$$ for i = 1:length(annotation.object)
% $$$   if isfield(annotation.object(i).world3d,'polygon3d')
% $$$     switch annotation.object(i).world3d.type
% $$$      case {'standingplanes','groundplane'}
% $$$       % Get set of planes and 3D points:
% $$$       PI = getPlanes(annotation.object(i).world3d.plane);
% $$$       X = getLMpolygon3D(annotation.object(i).world3d.polygon3d);
% $$$ 
% $$$       % Sample additional points along polygon boundary:
% $$$       X = [X X(:,1)];
% $$$       Xs = [];
% $$$       for j = 1:size(X,2)-1
% $$$         dd = sqrt(sum((X(:,j)-X(:,j+1)).^2))+eps;
% $$$         alpha = linspace(0,1,ceil(dd/sampDist)+1);
% $$$         alpha = alpha(2:end-1);
% $$$         Xs = [Xs [(1-alpha)*X(1,j)+alpha*X(1,j+1); (1-alpha)*X(2,j)+alpha*X(2,j+1); (1-alpha)*X(3,j)+alpha*X(3,j+1)]];
% $$$       end
% $$$       X = X(:,1:end-1);
% $$$ 
% $$$       keyboard;
% $$$       
% $$$       figure;
% $$$       plot3(X(1,:),X(2,:),X(3,:),'r.');
% $$$       hold on;
% $$$       plot3(Xs(1,:),Xs(2,:),Xs(3,:),'b.');
% $$$       
% $$$       % Sample along standing plane boundaries:
% $$$       
% $$$       
% $$$       % Iterate over planes:
% $$$       for j = 1:size(PI,2)
% $$$         % Get 3D polygon points that live on plane:
% $$$         if strcmp(annotation.object(i).world3d,'standingplanes')
% $$$           n = zeros(1,length(annotation.object(i).world3d.polygon3d.pt));
% $$$           for k = 1:length(annotation.object(i).world3d.polygon3d.pt)
% $$$             for ii = 1:1:length(annotation.object(i).world3d.polygon3d.pt(k).planeindex)
% $$$               if str2num(annotation.object(i).world3d.polygon3d.pt(k).planeindex(ii).index)==(j-1)
% $$$                 n(k) = 1;
% $$$               end
% $$$             end
% $$$           end
% $$$           n = find(n==1);
% $$$         else
% $$$           n = [1:length(annotation.object(i).world3d.polygon3d.pt)];
% $$$         end
% $$$         
% $$$         keyboard;
% $$$       
% $$$         figure;
% $$$         plot3(X(1,n),X(2,n),X(3,n));
% $$$       
% $$$ 
% $$$         pp = PI(:,j);
% $$$         R = eye(4);
% $$$         if abs(pp(1)) < 1e-4
% $$$           % Rotate axes to avoid numerical instabilities:
% $$$           R = [0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
% $$$           pp = R*pp;
% $$$           if abs(pp(1)) < 1e-4
% $$$             pp = R*pp;
% $$$             R = R*R;
% $$$           end
% $$$         end
% $$$ 
% $$$         % Get 2D to 3D plane transformation:
% $$$         T = [-pp(2)/pp(1) -pp(3)/pp(1) -pp(4)/pp(1); eye(3)];
% $$$         
% $$$         % Frontal-rectify plane points:
% $$$         Xj = X(:,n);
% $$$         xx = (T'*T)\T'*R*[Xj; ones(1,size(Xj,2))];
% $$$         xx = [xx(1,:)./xx(3,:); xx(2,:)./xx(3,:)];
% $$$       
% $$$         figure;
% $$$         plot(xx(1,:),xx(2,:));
% $$$ 
% $$$         xs = (T'*T)\T'*R*[Xs; ones(1,size(Xs,2))];
% $$$         xs = [xs(1,:)./xs(3,:); xs(2,:)./xs(3,:)];
% $$$         figure;
% $$$         plot(xx(1,:),xx(2,:),'r.');
% $$$         hold on;
% $$$         plot(xs(1,:),xs(2,:),'b.');
% $$$         
% $$$         % Convert back to 3D:
% $$$         XX = R'*T*[xx; ones(1,size(xx,2))];
% $$$         XX = [XX(1,:)./XX(4,:); XX(2,:)./XX(4,:); XX(3,:)./XX(4,:)];
% $$$       end
% $$$     end
% $$$   end
% $$$ end
% $$$ 
