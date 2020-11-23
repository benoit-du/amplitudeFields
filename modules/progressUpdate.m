function progressUpdate(idx,N)
%%% 22-11-20        first revision
%%% Benoit Duchet, University of Oxford

% progress bar inspired by http://undocumentedmatlab.com/blog/command-window-text-manipulation/

   msg = sprintf('Percent done: %3.1f', 100*idx/N);
   if idx == 1
       reverseStr = '';
   else
       reverseStr = repmat(sprintf('\b'), 1, length(msg));
   end
   fprintf([reverseStr, msg]);
end
