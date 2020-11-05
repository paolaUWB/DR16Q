# Normalized-Spectre

Welcome to Github!

Steps to contribute to this project:

1. Create a Github account. Once the account is created, the owner of the repositiory will need to give you access.

2. Once you have access to the repository, Press on the "Pull requests" on the top left tabs of the repository, and press the green button "New pull request on the top right , this allows any changes that you make to be on a seperate branch that will not affect the master file until it has been approved. 

3. To make an changes, you can choose to use other programs like Git Desktop for an easy interface or use Github itself by clicking on the file you want to make changes to in your branch, click the pencil at the top right, and paste the changed contents into the textbox.
  - Another way of doing this is pressing the "Push" button and it will allow you to paste the updated files into the branch.


## Code Notes

### Normalization File

-Please check comments very carefully.

-Follow the github steps above and open the file "norm_Sean_v20181029.py" using Spyder IDE / Visual Studio Code. That is the main code for normalization.

-We run this code under Spyder IDE, beacuse Spyder is a powerful scientific environment written in Python. It's designed by and for scientists, engineers and data analysts. It offers a unique combination of the advanced editing, analysis, debugging, and profiling functionality of a comprehensive development tool with the data exploration, interactive execution, deep inspection, and beautiful visualization capabilities of a scientific package. https://www.spyder-ide.org/

-Running on Visual Studio Code: Download Anaconda Navigator and launch the VS Code via Anaconda. This will you give an conda base interpreter. This is the easiest way to run Normalization file.
    
-There are 6760 spectrum (It can be more later on this project). If you don't generate all 6760 files, find variables "starts_from" and "ends_at" and define a range.

-After you run this code, all your graphs will be added to pdf file in your directory. Based on your range, you will see graphs in pdf files.

-You can see the final initinal parameters in additional file "Final_Initial_Parameters.txt" under same directory.

-If the program catches any error with the powerlaw, it will save it to the file "powerlaw_did_not_work.txt". Check that file under same directory.

-If you clone it correctly, then just RUN the file; then check the pdf files you just created with running the program.

-There is a test file in the directory. If you change something wrong, the test file will catch the different results. If you change constant variables, this will cause different output. If changes are correct, than update the test file with new results for future tests.

-Any major change should be added to README file.
 
