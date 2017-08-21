#!/usr/bin/env python3

def yesAnswer(question):  #asks the question passed in and returns True if the answer is yes, False if the answer is no, and keeps the user in a loop until one of those is given.  Also useful for walking students through basic logical python functions
    answer = False  #initializes the answer variable to false.  Not absolutely necessary, since it should be undefined at this point and test to false, but explicit is always better than implicit
    while not answer:  #enters the loop and stays in it until answer is equal to True
        print (question + ' (Y/N)')  #Asks the question contained in the argument passed into this subroutine
        answer = input('>>') #sets answer equal to some value input by the user
        if str(answer) == 'y' or str(answer) == 'Y' or str(answer).upper() == "YES":  #checks if the answer is a valid yes answer
            return True  #sends back a value of True because of the yes answer
        elif str(answer) == 'n' or str(answer) == 'N' or str(answer).upper() == "NO": #checks to see if the answer is a valid form of no
            return False  #sends back a value of False because it was not a yes answer
        else: #if the answer is not a value indicating a yes or no
            print ('Invalid response.')
            answer = False #set ansewr to false so the loop will continue until a satisfactory answer is given
            
class CheckArgs(object):
    
    def __init__(self):
        import os
        defaultInstallPath = os.environ["SCRATCH"] + os.sep + "veraT"
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument ("-d", "--installDirectory", help = "Print out the user manual and sample command lines.", default = defaultInstallPath)
        rawArgs = parser.parse_args()
        installDirectory = rawArgs.installDirectory
        if not os.path.isdir(installDirectory):
            os.makedirs(installDirectory)
        else:
            if not yesAnswer("Install directory %s already exists. Continue and overwrite any conflicts?"):
                raise FileExistsError("Install directory already exists and user does not wish to overwrite.")
        self.installDirectory = installDirectory

def main(args = False):
    import os
    import shutil
    if not args:
        args = CheckArgs()
    veraTRoot = "/u/nobackup/bound.T.hunter/mweinstein/"
    smallFolders = ["runners", "pipelines"]
    largeFolders = ["bin", "references"]
    for folder in smallFolders:
        source = veraTRoot + folder
        destination = args.installDirectory
        shutil.copytree(source, destination)
    for folder in largeFolders:
        source = veraTRoot + folder
        destination = args.installDirectory
        jobID = "cp" + folder
        command = 'echo "cp -r %s %s" | /u/systems/UGE8.0.1vm/bin/lx-amd64/qsub -cwd -V -N %s -l h_data=4G,time=24:00:00 -m bea' %(source, destination, jobID)
        print("COMMAND: %s" %command)
        os.system(command)
        
if __name__ == '__main__':
    main()
        