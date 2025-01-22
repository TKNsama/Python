#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   农艺性状调查.py
@Time    :   2022/04/14 15:17:04
@Author  :   Tan kenan 
@Version :   1.0
@Contact :   513060174@qq.com
'''

# here put the import lib
import myTTest
import os
import pandas as pd
# 参数-------------------------------------------------------
fileList = [ f for f in os.listdir("E:/Data/python/数据/农艺性状") if '.xlsx' in f]
file3smpList =  [f for f in os.listdir("E:/Data/python/数据/农艺性状/3smp") if '.xlsx' in f]
theDataNum = 3 # 从第几列开始计算
# -----------------------------------------------------------

for i in fileList:
    theDir = "E:/code/python/数据处理/农艺性状统计（t检验）/" + os.path.splitext(i)[0]
    file = "E:/Data/python/数据/农艺性状/" + i
    fileName = "E:/code/python/数据处理/农艺性状统计（t检验）/" + os.path.splitext(i)[0] + "/filtRes.xlsx"
    os.mkdir(os.path.splitext(i)[0])
    os.chdir(theDir)
    myTTest.meanToOutlier(theFile = file,theNum = theDataNum,method="sgm")
    myTTest.basicStats(theFile  = fileName,theNum = theDataNum + 1)
    myTTest.tTestInd(theFile = fileName, theNum = theDataNum + 1)
    filename = [ f for f in os.listdir(theDir) if '.xlsx' in f]
    index = 0
    dfs = []
    dfs = pd.DataFrame(dfs)
    for name in filename:
        print(index)
        dfs = pd.concat([dfs,pd.read_excel(os.path.join(theDir,name),header=None)],ignore_index=True,join="outer")
        index += 1
    print("合并了",index,"个表格")
    dfs.to_excel("results.xlsx")
    os.remove("filtRes.xlsx")
    os.remove("T-testResults.xlsx")
    os.remove("描述统计结果.xlsx")
    os.chdir("..")
'''
for i in file3smpList:
    theDir = "E:/code/python/数据处理/农艺性状统计（t检验）/" + os.path.splitext(i)[0]
    file = "E:/Data/python/数据/农艺性状/3smp/" + i
    fileName = "E:/code/python/数据处理/农艺性状统计（t检验）/" + os.path.splitext(i)[0] + "/filtRes.xlsx"
    os.mkdir(os.path.splitext(i)[0])
    os.chdir(theDir)
    myTTest.meanToOutlier(theFile = file,theNum = theDataNum,method="sgm")
    myTTest.basicStats(theFile  = fileName,theNum = theDataNum + 1)
    myTTest.tTestInd3samp(theFile = fileName, theNum = theDataNum +1)
    filename = [ f for f in os.listdir(theDir) if '.xlsx' in f]
    index = 0
    dfs = []
    dfs = pd.DataFrame(dfs)
    for name in filename:
        print(index)
        dfs = pd.concat([dfs,pd.read_excel(os.path.join(theDir,name),header=None)],ignore_index=True,join="outer")
        index += 1
    print("合并了",index,"个表格")
    dfs.to_excel("results.xlsx")
    os.remove("filtRes.xlsx")
    os.remove("T-testResults.xlsx")
    os.remove("描述统计结果.xlsx")
    os.chdir("..")
'''
