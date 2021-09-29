classdef idnlgreyfast < matlab.mixin.CustomDisplay
% class for nlgreyfast identification results
   properties
      param_est_v
      x0_est_v
      diag_data
   end
   methods
      function r = testFunction(obj)
         disp(obj.param_est_v)
      end
   end
end