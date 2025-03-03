import platform
import shutil
import os

def RemoveFile(path):
    if os.path.isfile(path):
        os.remove(path)
        print(f"Cleaned {path}.")
    else:
        print(f"No file named {path} was found.")
    return

def WrapperRemoveTree(path):
    try:
        shutil.rmtree(path)
        print(f"Cleaned {path}.")
    except FileNotFoundError:
        print(f"No directory named {path} was found!")
    except PermissionError:
        print(f"{path} could not be deleted because it is used by a process. \nTerminate it and relaunch the script to remove {path}.")
    return

def DarwinCleanMakefile():
    RemoveFile("Makefile")
    RemoveFile("Build/Core/Makefile")
    RemoveFile("Build/Sandbox/Makefile")
    WrapperRemoveTree(".vscode/")
    return

def DarwinCleanXCode():
    workspace = f"{projectName1}.xcworkspace/"
    project = f"Build/Core/{projectName1}.xcodeproj/"
    WrapperRemoveTree(workspace)
    WrapperRemoveTree(project)

    workspace = f"{projectName2}.xcworkspace/"
    project = f"Build/Sandbox/{projectName2}.xcodeproj/"
    WrapperRemoveTree(workspace)
    WrapperRemoveTree(project)

    print(f"Cleaned {workspace}, {project}.")
    return

def WinCleanVS():
    solution = f"Laniakea.sln"
    RemoveFile(solution)

    project = f"Build/Core/{projectName1}.vcxproj"
    user = f"Build/Core/{projectName1}.vcxproj.user"
    filt = f"Build/Core/{projectName1}.vcxproj.filters"
    vsfold = ".vs/"
    RemoveFile(project)
    RemoveFile(user)
    RemoveFile(filt)

    project = f"Build/Sandbox/{projectName2}.vcxproj"
    user = f"Build/Sandbox/{projectName2}.vcxproj.user"
    filt = f"Build/Sandbox/{projectName2}.vcxproj.filters"
    vsfold = ".vs/"
    RemoveFile(project)
    RemoveFile(user)
    RemoveFile(filt)

    WrapperRemoveTree(vsfold)
    return

if __name__ == "__main__":
    system = platform.system()
    projectName1 = "Library"
    projectName2 = "Sandbox"

    if system == "Darwin":
        if os.path.exists("Makefile") or os.path.exists(".vscode/"):
            DarwinCleanMakefile()
        if os.path.exists(f"{projectName1}.xcworkspace") and os.path.exists(f"{projectName2}.xcworkspace"):
            DarwinCleanXCode()

    elif system == "Windows":
        WinCleanVS()

    elif system == "Linux":
        DarwinCleanMakefile()

    else:
        print("Unidentified system. Halting.")
        exit(1)

    WrapperRemoveTree("Binaries/")
    choice = input("Remove documentation folder? [y/N] ")
    if choice == "y" or choice == "Y":
        WrapperRemoveTree("Docs")
    print("Nothing left to clean. Halting")
    exit(0)
