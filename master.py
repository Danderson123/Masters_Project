import subprocess

with open('D39V.fasta','r') as f:
    f = f.read()

# Set up the echo command and direct the output to a pipe
subprocess.run(['prokka', 'D39V.fasta'], capture_output = True)

