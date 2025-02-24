classdef Excel_utils
	% Summary of this class goes here
	%   Useful utilities for interacting with Excel via ActiveX.  Requires Windows operating system.
	% Sample class structure (if you want to add new properties or methods).
	% 	properties
	% 		raw;
	% 	end
	
	% 	methods(Static)
	% 		function obj = LoadFromExcel(excelFilepath)
	% 			[typ, desc, fmt] = xlsfinfo(excelFilepath);
	% 			[num, txt, raw]= xlsread(excelFilepath, desc{1});
	% 			obj = raw; % Return to caller
	% 			obj.raw = raw; % Can also, optionally, make it a property of this class.
	% 		end % function LoadFromExcel
	% 	end
	
	% This is a static method.  This means that you DO NOT have to do something like this:
	% excelUtilityClass = Excel_utils();
	% DON'T do that.  You simply call the class by itself, like the sample call below,
	% rather than doing something like excelUtilityClass.someMethod(arguments....)

	%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
	% SAMPLE USAGE           SAMPLE USAGE           SAMPLE USAGE           SAMPLE USAGE           SAMPLE USAGE           SAMPLE USAGE           SAMPLE USAGE
	% So, let's see how to actually use these methods in your MATLAB code.
	% Most of these static methods require an Excel object that you have to INSTANTIATE in advance before you call them.
	% Then you can CALL the function, and then you can SAVE your workbook, and finally DESTROY the Excel object.
		% Sample call
		% 		Excel = actxserver('Excel.Application');
		% 		excelWorkbook = Excel.Workbooks.Open(excelFullFileName);
		% 		Excel_utils.DeleteEmptyExcelSheets(Excel); % Delete any sheets (like the default "Sheet1" that don't have any cells filled on them.
		%		Excel_utils.AutoSizeAllSheets(Excel);
		% 		Excel.ActiveWorkbook.Save;
		% 		delete(Excel);
		% 		clear('Excel')
	% If you plan on calling multiple functions, don't instantiate and destroy for every function that you call, though you can - it won't hurt anything, it will just take longer to execute.
	% Only instantiate before you call the first method, then call as many methods as you want, then destroy only after you are done calling all methods.
	%--------------------------------------------------------------------------------------------------------------------------------------------------------------------

	
% 	Methods for class Excel_utils:
% 	Static methods:
% 
% 	ActivateSheet                          DeleteEmptyExcelSheets                 FormatDecimalPlaces                    LeftAlignSheet                         
% 	AlignCells                             DeleteExcelSheets                      FormatLeftBorder                       WrapText                               
% 	AutoSizeAllSheets                      DuplicateExcelSheet                    FormatRightBorder                      
% 	CenterCellsAndAutoSizeColumns          FormatBottomBorder                     GetNumberOfExcelSheets                 
% 	CenterCellsAndAutoSizeSpecificColumns  FormatCellColor                        GoToNextRowInColumn                    
% 	ClearCells                             FormatCellFont                         InsertComments      	
% 	AutoSizeColumns
	methods(Static)
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% GetNumberOfExcelSheets: returns the number of worksheets in the active workbook as a scalar, and their names in a cell array.
		% Sample call
		% 	Excel = actxserver('Excel.Application');
		%   Excel_utils.DuplicateExcelSheet(Excel, 'Results1', 'Results2');
		%   numberOfSheets = Excel_utils.GetNumberOfExcelSheets(Excel);
		function [numberOfSheets, sheetNames] = GetNumberOfExcelSheets(Excel)
			try
				worksheets = Excel.Sheets;
				numberOfSheets = worksheets.Count;
				sheetNames = cell(numberOfSheets, 1);
				for k = 1 : numberOfSheets
					sheetNames{k} = worksheets.Item(k).Name;
				end
			catch ME
				errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
					ME.stack(1).name, ME.stack(1).line, ME.message);
				fprintf(1, '%s\n', errorMessage);
				WarnUser(errorMessage);
			end
			return; % from GetNumberOfExcelSheets()
		end % of the GetNumberOfExcelSheets() method.

		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% DeleteEmptyExcelSheets: deletes all empty sheets in the active workbook.
		% This function loops through all sheets and deletes those sheets that are empty.
		% Can be used to clean a newly created xls-file after all results have been saved in it.
		% Sample call
		% 		Excel = actxserver('Excel.Application');
		% 		excelWorkbook = Excel.Workbooks.Open(excelFullFileName);
		% 		Excel_utils.DeleteEmptyExcelSheets(Excel);
		% 		Excel.ActiveWorkbook.Save;
		% 		Excel.Quit;
		% 		delete(Excel);
		% 		clear('Excel')
		function DeleteEmptyExcelSheets(Excel)
			try
				% 	Excel = actxserver('Excel.Application');
				% 	excelWorkbook = Excel.Workbooks.Open(fileName);
				
				% Run Yair's program http://www.mathworks.com/matlabcentral/fileexchange/17935-uiinspect-display-methods-properties-callbacks-of-an-object
				% to see what methods and properties the Excel object has.
				% 	uiinspect(Excel);
				
				worksheets = Excel.Sheets;
				sheetIndex = 1;
				sheetIndex2 = 1;
				initialNumberOfSheets = worksheets.Count;
				% Prevent beeps from sounding if we try to delete a non-empty worksheet.
				originalSoundSetting = Excel.EnableSound;
				Excel.EnableSound = false;
				% Tell it to not ask you for confirmation to delete the sheet(s).
				Excel.DisplayAlerts = false;
				
				% Loop over all sheets
				while sheetIndex2 <= initialNumberOfSheets
					% Saves the current number of sheets in the workbook.
					preDeleteSheetCount = worksheets.count;
					% Check whether the current worksheet is the last one. As there always
					% need to be at least one worksheet in an xls-file the last sheet must
					% not be deleted.
					if or(sheetIndex>1,initialNumberOfSheets-sheetIndex2>0)
						% worksheets.Item(sheetIndex).UsedRange.Count is the number of used cells.
						% This will be 1 for an empty sheet.  It may also be one for certain other
						% cases but in those cases, it will beep and not actually delete the sheet.
						if worksheets.Item(sheetIndex).UsedRange.Count == 1
							worksheets.Item(sheetIndex).Delete;
						end
					end
					% Check whether the number of sheets has changed. If this is not the
					% case the counter "sheetIndex" is increased by one.
					postDeleteSheetCount = worksheets.count;
					if preDeleteSheetCount == postDeleteSheetCount
						% If this sheet was not empty, and wasn't deleted, move on to the next sheet.
						sheetIndex = sheetIndex + 1;
					else
						% sheetIndex stays the same.  It's not incremented because the current sheet got deleted,
						% and all the other sheets shift down in their sheet number.  So now sheetIndex will
						% point to the same number which is the next sheet in line for checking.
					end
					sheetIndex2 = sheetIndex2 + 1; % prevent endless loop...
				end
				Excel.EnableSound = originalSoundSetting;
			catch ME
				errorMessage = sprintf('Error in function DeleteEmptyExcelSheets.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
				WarnUser(errorMessage);
			end
			return; % from DeleteEmptyExcelSheets
		end % of the DeleteEmptyExcelSheets() method.
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% DeleteExcelSheets: deletes sheets in the active workbook that have their name specified.
		% This function loops through all sheets and deletes those sheets whose name is in the caSheetNames cell array list.
		% Excel_utils.DeleteExcelSheets(Excel, {'Results1', 'Results2'});
		function DeleteExcelSheets(Excel, caSheetNames)
			try
				% 	Excel = actxserver('Excel.Application');
				% 	excelWorkbook = Excel.workbooks.Open(fileName);
				
				% Run Yair's program http://www.mathworks.com/matlabcentral/fileexchange/17935-uiinspect-display-methods-properties-callbacks-of-an-object
				% to see what methods and properties the Excel object has.
				% 	uiinspect(Excel);
				
				worksheets = Excel.sheets;
				initialNumberOfSheets = worksheets.Count;
				% Prevent beeps from sounding if we try to delete a non-empty worksheet.
				originalSoundSetting = Excel.EnableSound;
				Excel.EnableSound = false;
				% Tell it to not ask you for confirmation to delete the sheet(s).
				Excel.DisplayAlerts = false;
				
				% Loop over all the names
				for k = 1 : length(caSheetNames)
					% Get the current number of sheets in the workbook.
					preDeleteSheetCount = worksheets.count;
					% There must always be at least one worksheet in an xls-file, so the last sheet must not be deleted.
					if preDeleteSheetCount <= 1
						break;
					end
					% Loop over all the currently existing sheets, looking for this name.
					for sheetIndex = 1 : preDeleteSheetCount
						% Activate the worksheet.  (Perhaps unnecessary.)
% 						worksheets.Item(sheetIndex).Activate;
						% Get the name of the worksheet with this sheet index.
						thisName = worksheets.Item(sheetIndex).Name;
						% See if this name is in the caSheetNames list.
						itsInTheList = ismember(thisName, caSheetNames);

						% If it's in the list, delete it.
						if itsInTheList
							worksheets.Item(sheetIndex).Delete;
% 							postDeleteSheetCount = worksheets.count;
% 							fprintf('%d sheets left in workbook\n', postDeleteSheetCount);
							break;
						end
					end
				end
				% End up with the first sheet activated.
				worksheets.Item(1).Activate;
				Excel.EnableSound = originalSoundSetting;
			catch ME
				errorMessage = sprintf('Error in function DeleteExcelSheets.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
% 				WarnUser(errorMessage);
			end
			return; % from DeleteExcelSheets
		end % of the DeleteExcelSheets() method.
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% DuplicateExcelSheet: Duplicates the specified sheets in the active workbook and gives it the specified new name.
		% Sample call:
		% 	Excel = actxserver('Excel.Application');
		% 	excelWorkbook = Excel.workbooks.Open(fileName);
		%   Excel_utils.DuplicateExcelSheet(Excel, 'Results1', 'Results2');
		%	Duplicate the 'Results' workbook.
		%	Excel_utils.DuplicateExcelSheet(Excel, 'Results', sheetName);
		%	%numberOfSheets = Excel_utils.GetNumberOfExcelSheets(Excel);
		%	Save the workbook with the new sheet in it.
		%	Excel.ActiveWorkbook.Save;
		%	% Shut down Excel.
		%	Excel.Quit;
		%	delete(Excel);
		%	clear('Excel');
% 		Sheets("Results 1").Select
% 		Sheets("Results 1").Copy After:=Sheets(2)
% 		Sheets("Results 1 (2)").Select
% 		Sheets("Results 1 (2)").Name = "Results 2"
		function DuplicateExcelSheet(Excel, sourceSheetName, newSheetName)
			try
				Sheets = Excel.sheets;
% 				Excel.Visible = true;
				for sheetIndex = 1 : Sheets.count
					% Activate the worksheet.  (Perhaps unnecessary.)
					% 						worksheets.Item(sheetIndex).Activate;
					% Get the name of the worksheet with this sheet index.
					thisName = Sheets.Item(sheetIndex).Name;
					if strcmpi(thisName, sourceSheetName)
						% We found the sheet to copy.
						Sheets.Item(sheetIndex).Activate;
						% Run code from Mathworks technical support, on 11/9/2018, to duplicate a sheet.
                        MathWorks = get(Sheets, 'Item', sheetIndex);
                        MathWorks.Copy([], MathWorks);
                        Sheets.Item(sheetIndex+1).Name = newSheetName;
% 						copiedSheetName = sprintf('%s (2)', sourceSheetName);	% For example "Results 1 (2)"
% 						Sheets(copiedSheetName).Select
% 						Sheets(copiedSheetName).Name = newSheetName;
					end
				end
			catch ME
				errorMessage = sprintf('Error in function DuplicateExcelSheet.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
% 				User(errorMessage);
			end
			return; % from DuplicateExcelSheet
		end % of the DuplicateExcelSheet() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Add comments to cells on sheet.
		% Sometimes this throws exception #0x800A03EC on the second and subsequent images.  It looks like this:
		% "Error: Object returned error code: 0x800A03EC"
		% It is because of trying to insert a comment for a worksheet cell when a comment already exists for that worksheet cell.
		% So in that case, rather than deleting the comment and then inserting it, I'll just let it throw the exception
		% but I won't pop up any warning message for the user.
		function InsertComments(Excel, caComments, sheetNumber, startingRow, startingColumn)
			try
				worksheets = Excel.sheets;
				% 		thisSheet = get(worksheets, 'Item', sheetNumber);
				thisSheet = Excel.ActiveSheet;
				thisSheetsName = Excel.ActiveSheet.Name;  % For info only.
				numberOfComments = size(caComments, 1);  % # rows
				for columnNumber = 1 : numberOfComments
					columnLetterCode = cell2mat(ExcelCol(startingColumn + columnNumber - 1));
					% Get the comment for this row.
					myComment = sprintf('%s', caComments{columnNumber});
					% Get a reference to the cell at this row in column A.
					cellReference = sprintf('%s%d', columnLetterCode, startingRow);
					theCell = thisSheet.Range(cellReference);
					% You need to clear any existing comment or else the AddComment method will throw an exception.
					theCell.ClearComments();
					% Add the comment to the cell.
					theCell.AddComment(myComment);
				end
				
			catch ME
				errorMessage = sprintf('Error in function InsertComments.\n\nError Message:\n%s', ME.message);
				fprintf(errorMessage);
% 				WarnUser(errorMessage);
			end
			return; % from InsertComments
		end % of the InsertComments() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% places is the number of decimal places to the right of the decimal point, like 0.000 for 3.
		% excelRange is the range over which you want to apply the formatting, like 'A1..C5'
		% Example call:
		% Excel_utils.FormatDecimalPlaces(Excel, 3, 'A1..C5');
		function FormatDecimalPlaces(Excel, places, excelRange)
			try
				if places == 0
					formatString = '0';
				else
					formatString = '0.';
					for p = 1 : places
						% Append additional zeros.
						formatString = sprintf('%s0', formatString);
					end
				end
				
				% Select the range
				Excel.Range(excelRange).Select;
				
				% Format cells to the specified number of decimal places.
				Excel.Selection.NumberFormat = formatString;
			catch ME
				errorMessage = sprintf('Error in function FormatDecimalPlaces.\nThe Error Message:\n%s', ME.message);
				fprintf(errorMessage);
% 				WarnUser(errorMessage);
			end
			return; % from FormatDecimalPlaces
		end % of the FormatDecimalPlaces() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Turns on text wrap for cells in excelRange.
		% Example call:  Excel_utils.WrapText(Excel, 'A1..A12', true);
		function WrapText(Excel, excelRange, trueOrFalse)
			try
				% Select the range
				Excel.Range(excelRange).Select;
				
				% Turn wrapping on or off
				Excel.Selection.WrapText = trueOrFalse;
			catch ME
				errorMessage = sprintf('Error in function WrapText.\nThe Error Message:\n%s', ME.message);
				fprintf(errorMessage);
% 				WarnUser(errorMessage);
			end
			return; % from WrapText
		end % of the WrapText() method.
		

		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% borders is a collections of all. if you want, you can set one
		% particular border as,
		%
		% my_border = get(borders, 'Item', <item>);
		% set(my_border, 'ColorIndex', 3);
		% set(my_border, 'LineStyle', 9);
		%
		% where, <item> can be,
		% 1 - all vertical but not rightmost
		% 2 - all vertical but not leftmost
		% 3 - all horizontal but not bottommost
		% 4 - all horizontal but not topmost
		% 5 - all diagonal down
		% 6 - all diagonal up
		% 7 - leftmost only
		% 8 - topmost only
		% 9 - bottommost only
		% 10 - rightmost only
		% 11 - all inner vertical
		% 12 - all inner horizontal
		%
		% so, you can choose your own side.
		function FormatLeftBorder(sheetReference, columnNumbers, startingRow, endingRow)
			try
				numberOfColumns = length(columnNumbers);
				for col = 1 : numberOfColumns
					% Put a thick black line along the left edge of column columnNumber
					columnLetterCode = cell2mat(ExcelCol(columnNumbers(col)));
					cellReference = sprintf('%s%d:%s%d', columnLetterCode, startingRow, columnLetterCode, endingRow);
					theCell = sheetReference.Range(cellReference);
					borders = get(theCell, 'Borders');
					% Get just the left most border.
					leftBorder = get(borders, 'Item', 7);
					% Set it's style.
					set(leftBorder, 'LineStyle', 1);
					% Set it's weight.
					set(leftBorder, 'Weight', 4);
				end
				
			catch ME
				errorMessage = sprintf('Error in function FormatLeftBorder.\n\nError Message:\n%s', ME.message);
% 				WarnUser(errorMessage);
			end
			return; % from FormatLeftBorder
		end % of the FormatLeftBorder() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Format the right border with a thick line.
		% Sample call:
		%         Put a right border on the first column, the third column, and the last column (ONLY) going down all "numberOfRows" rows.
		%         Excel_utils.FormatRightBorder(Excel.ActiveSheet, [1, 3, numberOfColumns], 1, numberOfRows);
		function FormatRightBorder(sheetReference, columnNumbers, startingRow, endingRow)
			try
				numberOfColumns = length(columnNumbers);
				for col = 1 : numberOfColumns
					% Put a thick black line along the left edge of column columnNumber
					columnLetterCode = cell2mat(ExcelCol(columnNumbers(col)));
					cellReference = sprintf('%s%d:%s%d', columnLetterCode, startingRow, columnLetterCode, endingRow);
					theCell = sheetReference.Range(cellReference);
					borders = get(theCell, 'Borders');
					% Get just the left most border.
					leftBorder = get(borders, 'Item', 10);
					% Set it's style.
					set(leftBorder, 'LineStyle', 1);
					% Set it's weight.
					set(leftBorder, 'Weight', 4);
				end
				
			catch ME
				errorMessage = sprintf('Error in function FormatRightBorder.\n\nError Message:\n%s', ME.message);
% 				WarnUser(errorMessage);
			end
			return; % from FormatLeftBorder
		end % of the FormatRightBorder() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Format the bottom border with a line of specified thickness.
		% Sample call:
		%         Put a thick bottom border on the first row, the third row, and the last row (ONLY) going down all "numberOfRows" rows.
		%         Excel_utils.FormatBottomBorder(Excel.ActiveSheet, [1, 3, numberOfRows], 1, numberOfColumns, 4);
		% Weights can be from 1 to 4. (I think)
		function FormatBottomBorder(sheetReference, rowNumbers, startingCol, endingCol, weight)
			try
				numberOfRows = length(rowNumbers);
				for row = 1 : numberOfRows
					% Put a thick black line along the bottom edge of row rowNumbers(row)
					column1Letter = cell2mat(ExcelCol(startingCol));
					column2Letter = cell2mat(ExcelCol(endingCol));
					cellReference = sprintf('%s%d:%s%d', column1Letter, rowNumbers(row), column2Letter, rowNumbers(row));
					theCell = sheetReference.Range(cellReference);
					borders = get(theCell, 'Borders');
					% Get just the bottom most border.
					theBorder = get(borders, 'Item', 9);
					% Set it's style.
					set(theBorder, 'LineStyle', 1);
					% Set it's weight.
					set(theBorder, 'Weight', weight);
				end
				
			catch ME
				errorMessage = sprintf('Error in function FormatBottomBorder.\n\nError Message:\n%s', ME.message);
% 				WarnUser(errorMessage);
			end
			return; % from FormatBottomBorder
		end % of the FormatBottomBorder() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Selects all cells in the current worksheet and auto-sizes all the columns
		% and vertically and horizontally aligns all the cell contents.
		% Leaves with cell A1 selected.
		% Example call:  Excel_utils.CenterCellsAndAutoSizeColumns(Excel);
		function CenterCellsAndAutoSizeColumns(Excel)
			try
				% Select the entire spreadsheet.
				Excel.Cells.Select;
				% Auto fit all the columns.
				Excel.Cells.EntireColumn.AutoFit;
				% Center align the cell contents.
				Excel.Selection.HorizontalAlignment = 3;
				Excel.Selection.VerticalAlignment = 2;
				% Put "cursor" or active cell at A1, the upper left cell.
				Excel.Range('A1').Select;
			catch ME
				errorMessage = sprintf('Error in function CenterCellsAndAutoSizeColumns.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
% 				WarnUser(errorMessage);
			end
			return; % from CenterCellsAndAutoSizeColumns
		end % of the CenterCellsAndAutoSizeColumns() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Selects all cells in the current worksheet and auto-sizes all the columns.
		% Leaves with cell A1 selected.
		% Example call:  Excel_utils.AutoSizeColumns(Excel);
		function AutoSizeColumns(Excel)
			try
				% Select the entire spreadsheet.
				Excel.Cells.Select;
				% Auto fit all the columns.
				Excel.Cells.EntireColumn.AutoFit;
				% Put "cursor" or active cell at A1, the upper left cell.
				Excel.Range('A1').Select;
			catch ME
				errorMessage = sprintf('Error in function AutoSizeColumns.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
% 				WarnUser(errorMessage);
			end
			return; % from AutoSizeColumns
		end % of the AutoSizeColumns() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Selects specified columns in the current worksheet and auto-sizes just those columns
		% and vertically and horizontally aligns all the cell contents.
		% Leaves with cell A1 selected.
		% Example call:  Excel_utils.CenterCellsAndAutoSizeColumns(Excel, [1, 3]); % Autosize columns 1 and 3 only
		function CenterCellsAndAutoSizeSpecificColumns(Excel, columnNumbers)
			try
				numberOfColumns = length(columnNumbers);
				for col = 1 : numberOfColumns
					% Turn the column number into a letter code, for example 1=>A, and 27=>AB.
					columnLetterCode = cell2mat(ExcelCol(columnNumbers(col)));
					% Get the cell reference of this column.
					% The cell reference is just the column letter followed by a colon followed by the column letter again, like 'A:A'.
					cellReference = sprintf('%s:%s', columnLetterCode, columnLetterCode);
					% Auto fit only this particular column.
					Excel.Range(cellReference).Columns.AutoFit
					% Center align the cell contents.
					Excel.Selection.HorizontalAlignment = 3;
					Excel.Selection.VerticalAlignment = 2;
				end
				% Put "cursor" or active cell at A1, the upper left cell.
				Excel.Range('A1').Select;
			catch ME
				errorMessage = sprintf('Error in function CenterCellsAndAutoSizeColumns.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
% 				WarnUser(errorMessage);
			end
			return; % from CenterCellsAndAutoSizeSpecificColumns
		end % of the CenterCellsAndAutoSizeSpecificColumns() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Loops over all sheets in a workbook, auto-sizing columns and center-aligning all cells.
		function AutoSizeAllSheets(Excel)
			try
				% 	Excel = actxserver('Excel.Application');
				% 	excelWorkbook = Excel.workbooks.Open(fileName);
				%	Excel_utils.AutoSizeAllSheets(Excel);
				worksheets = Excel.sheets;
				numSheets = worksheets.Count;
				
				% Loop over all sheets
				for currentSheet = 1 : numSheets
					thisSheet = get(worksheets, 'Item', currentSheet);
					invoke(thisSheet, 'Activate');
					% Center data in cells, and auto-size all columns.
					Excel_utils.CenterCellsAndAutoSizeColumns(Excel)
				end
			catch ME
				errorMessage = sprintf('Error in function AutoSizeAllSheets.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
% 				WarnUser(errorMessage);
			end
			return; % from AutoSizeAllSheets
		end % of the AutoSizeAllSheets() method.
		

		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Tries to Activate the sheet.  Can pass in the number of the sheet (1,2,3,etc.) or name ('Results').
		function ActivateSheet(Excel, sheetNameOrNumber)
			try
				worksheets = Excel.sheets;
				numSheets = worksheets.Count;
				
				if isnumeric(sheetNameOrNumber)
					thisSheet = get(worksheets, 'Item', sheetNameOrNumber);
					thisSheet.Activate;
				else
					% Loop over all sheets looking for sheetname.
					for currentSheet = 1 : numSheets
						thisSheet = get(worksheets, 'Item', currentSheet);
						thisSheetName = strtrim(thisSheet.Name);
						if strcmpi(thisSheetName, sheetNameOrNumber)
							% Found the sheet we were looking for.  Activate it.
							thisSheet.Activate;
							break; % No need to keep looking.
						end
					end
				end
			catch ME
				errorMessage = sprintf('Error in function ActivateSheet.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
				WarnUser(errorMessage);
			end
			return; % from ActivateSheet
		end % of the ActivateSheet() method.
		

		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Left-align the specified sheet only.
		function LeftAlignSheet(Excel, sheetNumber)
			try
				sheetNumber3 = get(Excel.sheets, 'Item', sheetNumber);
				sheetNumber3.Activate;
				% Select the entire spreadsheet.
				Excel.Cells.Select;
				% Auto fit all the columns.
				% 		Excel.Cells.EntireColumn.AutoFit;
				% Left align the cell contents.
				Excel.Selection.HorizontalAlignment = 1;
				Excel.Selection.VerticalAlignment = 2;
				% Put "cursor" or active cell at A1, the upper left cell.
				Excel.Range('A1').Select;
			catch ME
				errorMessage = sprintf('Error in function LeftAlignSheet.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
				WarnUser(errorMessage);
			end
			return; % from LeftAlignSheet
		end % of the LeftAlignSheet() method.
		

		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Selects all cells in the current worksheet in the specified range.
		% Horizontally aligns all the specified cell range.
		% horizAlign = 1 for general (left alignment text but right align numbers).
		% horizAlign = 2 for left alignment.
		% horizAlign = 3 for center alignment.
		% horizAlign = 4 for right alignment.
		% Leaves with cell A1 selected.
		% Sample call Excel_utils.AlignCells(Excel, cellReference, 2, autoFit);
		function AlignCells(Excel, cellReference, horizAlign, autoFit)
			try
				% Select the range
				Excel.Range(cellReference).Select;
				% Align the cell contents.
				Excel.Selection.HorizontalAlignment = horizAlign;
				Excel.Selection.VerticalAlignment = 2;
				if autoFit
					% Auto fit all the columns.
					Excel.Cells.EntireColumn.AutoFit;
				end
				% Put "cursor" or active cell at A1, the upper left cell.
				Excel.Range('A1').Select;
			catch ME
				errorMessage = sprintf('Error in function AlignCells.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
				WarnUser(errorMessage);
			end % from AlignCells
			return;
		end % of the AlignCells() method.

		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Clears/erases cells in the current worksheet in the specified range.
		% Example call:
		% Excel_utils.ClearCells(Excel, 'A1..C5');
		function ClearCells(Excel, cellReference)
			try
				% Select the range
				Excel.Range(cellReference).Select;
				% Clear the cell contents.
				Excel.Selection.Clear;
				% Put "cursor" or active cell at A1, the upper left cell.
				Excel.Range('A1').Select;
			catch ME
				errorMessage = sprintf('Error in function ClearCells.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
				WarnUser(errorMessage);
			end % from ClearCells
			return;
		end % of the ClearCells() method.
		

		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Fills the background color of the Excel cell with the specified color.
		% Sample call:
		% Set the color of cells "A1..C1" of the active sheet to Yellow
		% cellReference = 'A1..C1';
		% Excel_utils.FormatCellColor(Excel, cellReference, 6);
		function FormatCellColor(Excel, cellReference, cellFillColorIndex)
			try
				Excel.ActiveSheet.Range(cellReference).Interior.ColorIndex = cellFillColorIndex;
			catch ME
				errorMessage = sprintf('Error in function FormatCellColor.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
				WarnUser(errorMessage);
			end % from FormatCellColor
			return;
		end % of the FormatCellColor() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Sample call to make font Calibri, size 12, blue, and bold.
		% Excel_utils.FormatCellFont(Excel, cellReference, 'Calibri', 12, true, 'b')
		function FormatCellFont(Excel, cellReference, fontName, fontSize, boldFace, fontColor)
			try
% 				worksheets = Excel.sheets;
				% 		thisSheet = get(worksheets, 'Item', sheetNumber);
% 				thisSheet = Excel.ActiveSheet;
% 				thisSheetsName = Excel.ActiveSheet.Name;  % For info only.
				% Select the range
				Excel.Range(cellReference).Select;
	% 			theCell = thisSheet.Range(cellReference);
				% Set the horizontal alignment to left justified.
% 				Excel.Selection.HorizontalAlignment = 2;
				% Set the font.
				Excel.Selection.Font.Name = fontName;
				% Set the font to bold.
				Excel.Selection.Font.Bold = boldFace;
				% Set the font size to 12 points.
				Excel.Selection.Font.Size = fontSize;
				if fontColor == 'b'
					% Set the font color to blue.
					Excel.Selection.Font.Color = -65536;
				else
					% Set the font color to black.
					Excel.Selection.Font.Color = 0;
				end

			catch ME
				errorMessage = sprintf('Error in function FormatCellFont.\n\nError Message:\n%s', ME.message);
				fprintf(errorMessage);
				WarnUser(errorMessage);
			end
			return; % from FormatCellFont
		end % of the FormatCellFont() method.
		
		
		%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
		% Returns the next empty cell in column after row 1.  Basically it puts the active cell in row 1
		% and types control-(down arrow) to put you in the last row.  Then it adds 1 to get to the next available row.
		% Sample call:
		% nextRow = Excel_utils.GoToNextRowInColumn(Excel, 'A')
		function nextRow = GoToNextRowInColumn(Excel, column)
			try
				nextRow = -1;
				if isnumeric(column)
					% If they passed in a number, convert it to a column letter.
					column = char(ExcelCol(column));
				end
				% Make a reference to the very last cell in this column.
				cellReference = sprintf('%s1048576', column);
				Excel.Range(cellReference).Select;
				currentCell = Excel.Selection;
				bottomCell = currentCell.End(3); % Control-up arrow.  We should be in row 1 now.
				% Well we're kind of in that row but not really until we select it.
				bottomRow = bottomCell.Row;
				cellReference = sprintf('%s%d', column, bottomRow);
				Excel.Range(cellReference).Select;
				bottomCell = Excel.Selection;
				bottomRow = bottomCell.Row;  % This should be the last row
				% If this cell is empty, then it's the next row.
				% If this cell has something in it, then the next row is one row below it.
				cellContents = Excel.ActiveCell.Value;  % Get cell contents - the value (number of string that's in it).
				% If the cell is empty, cellContents will be a NaN.
				if isnan(cellContents)
					% Row 1 is empty.  Next row should be 1.
					nextRow = bottomRow;  % Don't add 1 since it was empty (the top row already).
				else
					% Row 1 is not empty.  Next row should be row 2.
					nextRow = bottomRow + 1;  % Will add 1 to get row 1 as the next row.
				end
			catch ME
				errorMessage = sprintf('Error in function GoToNextRowInColumn.\n\nError Message:\n%s', ME.message);
				fprintf('%s\n', errorMessage);
				WarnUser(errorMessage);
			end
			return; % from GoToNextRowInColumn() method.			
		end % of the GoToNextRowInColumn() method.
		
		
	end % End of the methods definitions.
	
end % of the Excel_utils class.
