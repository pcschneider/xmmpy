#!/usr/bin/python3
import subprocess

python_call = "/home/majestix/hdd/python/bin/python3.12"
xmmpy_root = "/home/majestix/hdd/tools/xmmpy/"

if __name__ == "__main__":
    import sys
    print(sys.argv)
    subprocess.run([python_call, xmmpy_root+"/bin/xmm_image_pdf.py", *sys.argv[1:]])
    