classdef StrSymbolic < handle
  %% STRSYMBOLIC This class implements symbolic operations but using strings
  % instead of using the symbolic toolbox. The functionallity developed in this
  % class is faster than the symbolic toolbox, but it is at the expense of using 
  % strings and not all the functionallity of the symbolic toolbox will be
  % implemented here. Use this class only for optimization purposes.
  %

  methods (Static)
    function [out] = symvar(in)
      %% SYMVAR Find symbolic variables in string input.
      %
      % param: in [char] String input.
      %
      % return: out {[char]} String variables.

      out = {};

      aux = in;

      ind = false(size(aux));
      % Search for everything but words.
      % Words cannot start with a number but they can contain numbers.
      ind(regexp(aux,'\W')) = true;
      aux(ind)= ' ';
      % Remove numbers that follow an empty space.
      while 1
        ind_old = ind;
        ind(regexp(aux,'(?<=\s)\d')) = true;
        aux(ind)= ' ';
        offset = 0;
        if ~sum(ind ~= ind_old)
          break;
        end
      end

      for i = regexp(aux,'\w*')
        words = (split(aux(i:end)));
        word = words{1};

        % Remove builtin variables.
        if (exist(word,'builtin')==5)
          continue;
        end
        
        % Remove names of files with .m extension.
        if exist(word) == 2
            continue;
        end

        % Remove the time variable.
        if (strcmp(word,'t'))
          continue;
        end
        
        % Remove not valid var names.
        if ~isvarname(word)
            continue;
        end

        % Remove the exponencial.
        if word(1) == 'e' 
          if length(word) == 1
            continue;
          end
          if ~isnan(str2double(word(2:end)))
            continue;
          end
        end

        % Check if var was already find.
        findName = false;
        for j = 1:length(out)
          if strcmp(out{j},word)
            findName = true;
            break;
          end
        end

        if ~findName
          out{end+1} = word;
        end
      end

    end % symvar

    function [out] = subs(s,old,new)
      %% SUBS Symbolic substitution.
      %
      % param: s {[char]} String with the expression.
      %      : old {[char]} Variable to be substituted in s.
      %      : new {[char]} New value for old.
      %
      % return: out {[char]} String with the substitution done.

      aux = s;

      ind = false(size(aux));
      % Search for everything but words.
      % Words cannot start with a number but they can contain numbers.
      ind(regexp(aux,'\W')) = true;
      aux(ind)= ' ';
      % Remove numbers that follow an empty space.
      while 1
        ind_old = ind;
        ind(regexp(aux,'(?<=\s)\d')) = true;
        aux(ind)= ' ';
        offset = 0;
        if ~sum(ind ~= ind_old)
          break;
        end
      end
      
      [newStr, match] = split(aux);
      
      j = length(s);
      
      for i = length(newStr):-1:1
          jNext = j - length(newStr{i}) +1;
          
          subsVar = s(jNext:j);
          
          ind = strcmp(subsVar,old);
          
          if any(ind)
              s = [s(1:jNext-1) new{ind} s(j+1:end)];           
          end
          
          if i > 1
            j = j - length(newStr{i}) - length(match{i-1});
          else
            j = j - length(newStr{i});  
          end
          
      end
      
      out = s;

    end % subs

  end % methods

end % classdef
