function nextStartRow = writeExcelTable(T, filePath, sheetName, startRow, heading)
    % Print heading
    if ~isempty(heading)
        headingRange = sprintf('A%d', startRow);
        writecell({heading}, filePath, 'Sheet', sheetName, 'Range', headingRange);
        startRow = startRow+1;
    end
    % Specify the range for writing the table
    range = sprintf('A%d', startRow);
    writetable(T, filePath, 'Sheet', sheetName, 'Range', range);

    % Calculate the next starting row for the next table (current table size + 2 for spacing)
    nextStartRow = startRow + size(T, 1) + 2;
end