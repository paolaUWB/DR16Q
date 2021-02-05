# Normalized-Spectra

Welcome to Github!

Steps to contribute to this project:

1. Create a Github account. Once the account is created, the owner of the repository will need to give you access.

2. Once you have access to the repository, Press on the "Pull requests" on the top left tabs of the repository, and press the green button "New pull request on the top right, this allows any changes that you make to be on a separate branch that will not affect the master file until it has been approved. 

3. To make changes, you can choose to use other programs like Git Desktop for an easy interface or use Github itself by clicking on the file you want to make changes to in your branch, click the pencil at the top right, and paste the changed contents into the textbox.
  - Another way of doing this is pressing the "Push" button and it will allow you to paste the updated files into the branch.


## General Notes

-We run this code under Spyder IDE, because Spyder is a powerful scientific environment written in Python. It's designed by and for scientists, engineers, and data analysts. https://www.spyder-ide.org/

-Running on Visual Studio Code: Download Anaconda Navigator and launch the VS Code via Anaconda. This will you give a conda base interpreter.

-Any major change should be added to the README file.


### Normalization File

-Follow the GitHub steps above and open the file "normalization.py" using Spyder IDE or Visual Studio Code.
    
-There are 6760 spectra (It can be more later on this project). If you don't generate all 6760 files, find variables "STARTS_FROM" and "ENDS_AT" and define a range.

-The "data_types" and "utility_functions" files must be in your directory for import

-The file currently named "sorted_norm.csv" contains all the necessary information to read.

- At the top of the file, the constant variables are defined by Astrophysicist. Any change of these constants SHOULD be discussed with the client.

-File currently has DR9 "Data Release 9" extension. In the future, it will be changed to the recent Data Release from the SDSS Database.

-After you run this code, all your graphs will be added to a pdf file in your directory. Based on your range, you will see graphs in pdf files.

-If you clone it correctly, then just RUN the file; then check the pdf files you just created by running the program.

-There is a test file in the directory. If you change something wrong, the test file will catch the different results. If you change constant variables, this will cause different output. If changes are correct then update the test file with new results for future tests.


### Absorption File

-Follow the GitHub steps above and open the file "absorption.py" using Spyder IDE or Visual Studio Code.

-File currently has DR9 "Data Release 9" extension. In the future, it will be changed to the recent Data Release from SDSS Database
 
-The Absorption file reads the Normalization files. After you generate the normalization file for each spectra, the Absorption file will read that and generate necessary graphs and information. Currently the format of the processed Normalization file is : "spectra name + norm + .dr9"

-Based on "BALNICITY_INDEX_LIMIT", the file will import different results. You can see the results in the same directory under the "Absortion_cleaning" document. Also the same "BALNICITY_INDEX_LIMIT" constant will generate different graph output because it is defining the narrowness of the output.

- At the top of the file, the constant variables are defined by Astrophysicist. Any change of these constants SHOULD discuss with the client.