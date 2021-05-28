classdef ClassCommand < LineCommand

  properties 
    % [char] Name used for the command. Name is auto-included to keywords.
    name = 'class';
    % {[char]} struct with the list of keywords that must be reserved for this command.
    keywords = {};
    % [char] End sequence of the command. It is defined dynamically.
    endWith; 

  end % properties 

  methods

    function [obj] = ClassCommand(mcp)
      %% Constructor of ClassCommand.
      %
      % param: mcp  OneModelParser object.

      obj = obj@LineCommand(mcp);

      % We want to execute this command with "use" command.
      obj.execUse = true;

      obj.introEnd = false;

    end % LineCommand


    function out = isComplete(obj, raw)
      %% ISCOMPLETE Does the command have everything it needs to run?
      %
      % param: raw Raw text from the OneModel file.
      %
      % return: true if the argument is complete.

      [tokens] = regexp(raw,'\s*class\s\s*(\w*)\s*','tokens');

      if isempty(tokens)
        out = false;
        return;
      end

      obj.endWith = ['\s*end\s\s*' tokens{1}{1} '\s*;'];

      out = isComplete@LineCommand(obj, raw);

    end


    function [] = execute(obj, raw)
      %% EXECUTE Execute the command.
      %
      % param: raw  Raw text from the OneModel file.
      %
      % return: true if the argument is complete.

      [tokens] = regexp(raw,'\s*class\s\s*(\w*)\s*','tokens');

      obj.mcp.className{end+1} = tokens{1}{1};

      expr = ['\s*class\s\s*' tokens{1}{1} '\s*([\s\S]+)end\s\s*' tokens{1}{1} ';'];

      [tokens] = regexp(raw,expr,'tokens');

      obj.mcp.classCode{end+1} = tokens{1}{1};

    end % execute

  end % methods

end % classdef


