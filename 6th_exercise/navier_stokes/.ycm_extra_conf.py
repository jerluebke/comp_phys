# -*- coding: utf-8 -*-

import os

FLAGS = [
        '-Wall',
        '-Wextra',
        #  '-Werror',
        #  '-fexceptions',
        #  '-ferror-limit=10000',
        #  '-DNDEBUG',
        '-std=c++17',
        '-xc++',
        #  '-isystem', 'C:\\tools\\Anaconda3\\include',
        #  '-isystem', 'C:\\Program Files\\LLVM\\include',
        #  '-isystem', 'D:\\source\\gsl\\build-dir',
        #  '-isystem', 'D:\\source\\opencv\\build\\include',
        #  '-isystem', 'C:\\MinGW64\\mingw64\\include',
        #  '-isystem',
        #      'D:\\source\\qt\\qt\\5.11.1\\mingw53_32\\include\\QtWidgets',
        '-I', '.\\include',
]


def readVtkPaths(vtkRoot):
    global FLAGS
    pathDict = {path : (direc, file) for path, direc, file in os.walk(vtkRoot)}
    for direc in pathDict[vtkRoot][0]:
        tmpPath = os.path.join(vtkRoot, direc)
        for subDir in pathDict[tmpPath][0]:
            FLAGS.extend(['-isystem', os.path.join(tmpPath, subDir)])


def Settings(**kwds):
    readVtkPaths('D:\\source\\Libs\\VTK-8.1.2')
    readVtkPaths('D:\\source\\Libs\\VTK-8.1.2\\build')
    return {'flags' : FLAGS}
