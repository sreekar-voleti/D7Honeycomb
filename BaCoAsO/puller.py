import sys

ty = input ("File (a) or Folder (b)?")

if ty == "a":

    file_path = input("Path to file: ")
    destination = input("Where to store file: ")

    f = open("puller.sh" , "w")

    f.write("scp sreekarv@niagara.computecanada.ca:/scratch/a/aparamek/sreekarv/" + file_path + " " + destination)

    f.close()

elif ty == "b":

    folder_name = input("Name of folder to transfer: ")
    destination = input("Where to store folder: ")

    f = open("puller.sh" , "w")

    f.write("scp -r sreekarv@niagara.computecanada.ca:/scratch/a/aparamek/sreekarv/" + folder_name + " " + destination)

    f.close()

else:
    print("You need in input 'a' or 'b', man")
