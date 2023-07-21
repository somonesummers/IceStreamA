clear
cd data/
files = dir;
directoryNames = {files.name};
directoryNames = directoryNames(~ismember(directoryNames,{'._.DS_Store','.','..','.DS_Store','PreApril2023'}));
currDir = pwd;
cd ..
for i = 1:numel(directoryNames)
    filename = append(currDir,'/',directoryNames{i});
%     disp(filename)
   try
        data = load(filename);
        fix = 0;
        try
            fields = {'bed','bed_interp','surface','surf_interp','thick','thick_interp'};
            data = rmfield(data,fields);
            fix = fix + 1;
        catch
            % no action if Goll fields already lacking
        end
        try
            fields = {'dhdt_neg','dhdt_pos','thick_neg','thick_pos','thick_interp'};
            data = rmfield(data,fields);
            fix = fix + 1;
        catch
            % no action if DhDt fields already lacking
        end

        if(fix == 1)
            save(filename,'-struct','data');
            disp("Succesfully reduced " + filename + " fix = " + fix)
        else
            disp("No changes made for " + filename + " fix = " + fix)
        end
   catch
        disp("Could not load item " + filename + ". Skipping")
   end
end
