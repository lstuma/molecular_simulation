
colors = {
          "info": "\033[;34m",
          "warning": "\033[;33m",
          "warning2": "\033[;31m",
          "success": "\033[;32m\033[2m",
          "error": "\033[;31m",
          "critical": "\033[1;31m",
          "reset": "\033[0;0m"
}

prefix = {
          "info": "\033[;34mINFO          ",
          "warning": "\033[;33mWARNING       ",
          "warning2": "\033[1;31mCRIT \033[0;33mWARNING  ",
          "success": "\033[;32mSUCCESS       ",
          "error": "\033[;31mERROR         ",
          "critical": "\033[1;31mCRITICAL      "
}


def log(msg, level="info", end="\n", begin=""):
    print(begin+colors[level]+prefix[level]+"\t\033[;37m|\033[0;0m "+msg+"\033[0;0m", end=end)

def log_up(lines=1):
    print('\x1B['+str(2+lines)+'A')

def log_down(lines=1):
    print('\n'*lines)