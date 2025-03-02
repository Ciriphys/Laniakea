import platform
import shutil
import os

# TODO : Change Project Name

def XCodeSetPCHLocation():
    filedata = ""
    with open("Build/Laniakea.xcodeproj/project.pbxproj", "r") as file:
        filedata = file.read()
    
    filedata = filedata.replace("GCC_PREFIX_HEADER = lnkpch.h", "GCC_PREFIX_HEADER = Include/lnkpch.h")

    with open("Build/Laniakea.xcodeproj/project.pbxproj", "w+") as file:
        file.write(filedata)

def GenerateDarwinDocs():
    os.system("./Vendor/doxygen/macOS/doxygen Doxyfile")
    os.chdir("Docs/PDF/")
    os.system("make pdf")
    os.system("mv refman.pdf ../Manual.pdf")
    os.chdir("../../")
    os.system("rm Docs/PDF/*")
    os.system("mv Docs/Manual.pdf Docs/PDF/")

def GenerateWinDocs():
    os.system("\"Vendor\\doxygen\\Windows\\doxygen.exe\" Doxyfile")
    os.chdir("Docs\\PDF")
    os.system("make.bat")
    os.system("move refman.pdf ..\\Manual.pdf")
    os.chdir("..\\..")
    
    try:
        shutil.rmtree("Docs\\PDF")
    except PermissionError:
        print("Docs\\PDF could not be deleted because it is used by a process. \nTerminate it and relaunch the script.")
    finally: 
        os.system("mkdir Docs\\PDF")
        os.system("move Docs\\Manual.pdf Docs\\PDF\\")

if __name__ == "__main__":
    system = platform.system()
    command = ""
    selection = ""

    if system == "Darwin":
        selection = input("Enter the action (gmake2 ~ Makefile | xcode4 ~ XCode): ")
        command = f"./Vendor/premake5/macOS/premake5 {selection}"
    elif system == "Windows":
        command = str("\"Vendor\\premake5\\Windows\\premake5.exe\" vs2022")
    elif system == "Linux":
        command = "./Vendor/premake5/Linux/premake5 gmake2"
    else:
        print("Unidentified system. Halting.")
        exit(1)

    print(command)
    os.system(command)
    
    if selection == "xcode4":
        XCodeSetPCHLocation()

    choice = input("Generate documentation? [y/N] ")

    if choice == "y" or choice == "Y":
        if system == "Darwin":
                GenerateDarwinDocs()
        elif system == "Windows":
                GenerateWinDocs()
        elif system == "Linux":
                print("Linux is not yet supported! Halting.")
        else:
            print("Unidentified system. Halting.")
            exit(1)
    else:
        exit(0)