import subprocess

bashCommand = "cd /H/FU/FLWO/"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
bashCommand = "htr HTR111-006"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()