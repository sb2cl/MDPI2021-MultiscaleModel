classdef OneModelParser < handle
  %% MODELCLASSPARSER This class parses files .mc into OneModel models.
  %

  properties
    % [char] Name of the model to parse.
    filename
    % [char] Base name of the model without extension.
    basename
    % [char] Name for the OneModel model.
    nameM
    % File descriptor of the ouput model file.
    fout
    % {Command} List of command defined for the parser.
    commands = {}
    % {filename} Filename of extended files. This is to avoid infinite recursion
    % of files
    filenameExtended = {}
    % {char} Cell with all the path to the files that make up the model.
    dependenciesPath = {}
    % {[char]} List of name of defined classes.
    className = {}
    % {[char]} The OneModel code of the defined classes.
    classCode = {}

  end % properties

  methods 

    function [obj] = OneModelParser(filename)
      %% Constructor of OneModelParser.
      %
      % param: filename File to parse into a OneModel.

      obj.filename = filename;

      [folder, baseFileNameNoExt, extension] = fileparts(obj.filename);

      obj.basename = baseFileNameNoExt;

      obj.nameM = [baseFileNameNoExt '.m'];

      % Add here all the commands of the parser.
      obj.commands = {
        TestCommand(obj),
        VariableCommand(obj),
        ParameterCommand(obj),
        EquationCommand(obj),
        ImportCommand(obj),
        MatlabCodeCommand(obj),
        SimOptionsCommand(obj),
        ClassCommand(obj),
        NamespaceCommand(obj),
        ObjectCommand(obj),
        UseCommand(obj),
        ExtendsCommand(obj),
        EndIfCommand(obj),
        IfCommand(obj),
        ConnectCommand(obj),
        ChangeAtributeCommand(obj)
      };

    end % OneModelParser

    function [] = parse(obj)
      %% PARSE Parse the file into the a OneModel.
      %
      % return: void

      % Open the model.
      obj.avoidRecursion(obj.filename);
      fid = fopen(obj.filename);

      % Check if build directory exists.
      if ~exist('build', 'dir')
       mkdir('build')
      end

      addpath('build');

      % Open the generated model.

      obj.fout = fopen(['build/' obj.nameM],'w');

      obj.addHeader(obj.fout);
      obj.executeFileLines(fid,obj.fout);
      obj.addFooter(obj.fout);

      fclose(fid);
      fclose(obj.fout);

    end % parse

    function [] = executeFileLines(obj,fid,fout,opt)
      %% EXECUTEFILELINES 
      %
      % param: fid File descriptor to read.
      %      : fout File descriptor to write.
      %      : opt [char] options for the changing the behavior.
      %
      % return: void

      if nargin < 4
        opt = '';
      end

      % Read a line of the model.
      tline = fgets(fid);

      % Aux buffer for preparing the lines to be executed.
      aux = '';

      % Variable to store the current command to execute.
      cmd = [];

      while ischar(tline)
        % Remove commented text in lines.
        tline = obj.removeComments(tline);

        % Remove empty lines.
        if isempty(tline)
          tline = fgets(fid);
          continue;
        end

        % Process char by char the line readed.
        for i = 1:length(tline)
          aux(end+1) = tline(i);

          % If we did not find a command.
          if isempty(cmd)
            % Look if any of the commands in the list is found.
            for i = 1:length(obj.commands)
              % If any of the commands if found.
              if obj.commands{i}.findCommand(aux)
                % Save it and stop looking for more commands.
                cmd = obj.commands{i};
                break;
              end
            end
          end

          % If we have found a command
          if ~isempty(cmd)
            % Collect all the argument data needed for the command in aux.

            % Does the command all it needs to be executed? 
            if cmd.isComplete(aux)
              if ~strcmp(opt,'execUse') || (strcmp(opt,'execUse') && cmd.execUse)
                % Execute the comand.
                cmd.execute(aux);

                if cmd.introEnd
                  fprintf(fout,'\n');
                end

              end

                % Reset the aux for new lines.
                aux = '';
                % Reset the cmd for new commands.
                cmd = [];
              
            end
          end
        end

        % Read the next line in the document.
        tline = fgets(fid);
      end

      % Check if something is remaining in aux after finishing the model.
      if ~all(isspace(aux)) && ~isempty(aux)
        try
          feval(cmd,obj,aux,-1);
        catch
          error('The following code did not correspond to any OneModel command: \n--- Start of code error ---\n%s\n--- End of code error ---',aux);
        end
      end
    end % executeFileLines

    function [out] = removeComments(obj,tline)
      %% REMOVECOMMENTS Remove commented text of the model text.
      %
      % param: tline Line where to find and remove comments.
      %
      % return: out Line without comments.

      ind = strfind(tline,'%');

      if ~isempty(ind)
        tline = tline(1:ind(1)-1);
      end

      out = tline;

    end % removeComments

    function [] = print(obj,str)
      %% PRINT Prints 'str' into the output file.
      %
      % param: str String to be printed.
      
      fprintf(obj.fout,['\t\t\t' str]);
      
    end % print

    function [] = addHeader(obj,fout)
      %% ADDHEADER Add the header to the Matlab Class.
      %
      % param: fout File output.
      %
      % return: void

      fprintf(fout,'classdef %s < OneModel\n',obj.basename);
      fprintf(fout,'\t%% This code was generated by OneModel %s\n', OneModel.version);
      fprintf(fout,'\tmethods\n');
      fprintf(fout,'\t\tfunction [obj] = %s(opts)\n',obj.basename);

    end % addHeader

    function [] = addFooter(obj,fout)
      %% ADDFOOTER Add the footer to the Matlab Class.
      %
      % param: fout File output.
      %
      % return: void

      fprintf(fout,'\t\t\tobj.checkValidModel();\n');
      fprintf(obj.fout,'\t\tend\n\n');
      fprintf(fout,'\tend\n');


      fprintf(fout,'\tmethods(Static)\n');

      fprintf(fout,'\t\tfunction [out] = isUpToDate()\n');
      fprintf(fout,'\t\t\tdependenciesPath = {...\n');
      for i = 1:length(obj.dependenciesPath)
      fprintf(fout,'\t\t\t\t''%s''...\n',obj.dependenciesPath{i});
      end
      fprintf(fout,'\t\t\t};\n');
      fprintf(fout,'\t\t\tout = %s.checkUpToDate(dependenciesPath);\n',obj.basename);
      fprintf(fout,'\t\tend\n');
      fprintf(fout,'\tend\n');

      fprintf(fout,'end\n');

    end % addFooter

    function [] = avoidRecursion(obj,filename)
      %% AVOIDRECURSION Check if the filename was already used to avoid
      % recursion.
      %
      % param: filename Filename to check.
      %
      % return: void

      [pathstr,name,ext] = fileparts(filename);

      if any(strcmp(obj.filenameExtended,name),'all')
        % error('The file "%s" was already included in the model. The parsing of the model was stoped to avoid infinite recursion.',name);
        disp(sprintf('The file "%s" was already included in the model.',name));
        return;
      end

      obj.filenameExtended{end+1} = name;
      obj.dependenciesPath{end+1} = filename;

    end % avoidRecursion

  end %methods

end % classdef
