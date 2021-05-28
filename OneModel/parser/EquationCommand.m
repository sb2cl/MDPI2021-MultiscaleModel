classdef EquationCommand < LineCommand

  properties 
    % [char] Name used for the command. Name is auto-included to keywords.
    name = 'equation';
    % {[char]} struct with the list of keywords that must be reserved for this command.
    keywords = {};
    % [char] End sequence of the command.
    endWith = ';';
    
  end % properties 

  methods

    function [] = execute(obj, raw)
      %% EXECUTE Execute the command.
      %
      % param: raw  Raw text from the OneModel file.
      %
      % return: true if the argument is complete.

      % Remove intros.
      raw = raw(raw~=newline);
      arg = obj.getArgument(raw);
      [name,options] = obj.getOptions(arg);

      if isempty(options{1})
        options{1} = name;
        name = '';
      end

      fprintf(obj.mcp.fout,'\t\t\te = EquationClass(obj,''%s'');\n',name);      

      try
        fprintf(obj.mcp.fout,'\t\t\te.eqn = ''%s'';\n',options{1});
      catch
        error('eqn is not defined in the options.');
      end

      for i=2:length(options)
        % Skip empty options.
        if isempty(options{i})
          continue
        end

        fprintf(obj.mcp.fout,'\t\t\te.%s;\n',options{i});

      end

      fprintf(obj.mcp.fout,'\t\t\tobj.addEquation(e);\n');

      end % execute

  end % methods

end % classdef


