import subprocess

o = subprocess.Popen(["vrp_release_and_delays.exe", "problemFile.txt"])
o.wait()