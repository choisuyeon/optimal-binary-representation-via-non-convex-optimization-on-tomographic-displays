function [currentStringCommand] = print2cmd(currentStringCommand, message)
    message = strrep(message, '%', '%%');
    message = strrep(message, '\', '\\');
    fprintf([num2str(currentStringCommand), message]);    
    currentStringCommand = repmat(sprintf('\b'), 1, length(message));    
end

