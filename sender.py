import sys

ty = input ("File (a) or Folder (b)?")

if ty == "a":

    file_name = input("Name of file to transfer: ")
    folder_name = input("Name of folder to transfer to: ")

    f = open("sender.sh" , "w+")

    f.write("scp "+ file_name + " sreekarv@niagara.computecanada.ca:/scratch/a/aparamek/sreekarv/" + folder_name)

    f.close()

elif ty == "b":

    file_name = input("Name of folder to transfer: ")
    folder_name = input("Name of folder to transfer to: ")

    f = open("sender.sh" , "w+")

    f.write("scp -r "+ file_name + " sreekarv@niagara.computecanada.ca:/scratch/a/aparamek/sreekarv/" + folder_name)

    f.close()

else:
    print("You need in input 'a' or 'b', man")



# ++++++
