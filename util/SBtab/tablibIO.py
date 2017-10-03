#!/usr/bin/python

import mimetypes
import sys
import tablib
import tablib.core
import csv
import os
from util.SBtab import misc
import xlrd


def sheets(self):  # Added to excess sheets of Databook
    return self._datasets
try:
    tablib.Databook.sheets
except:
    tablib.Databook.sheets = sheets

def importSetNew(sbtabfile,filename,separator=None):
    mimetypes.init()
    file_mimetype = mimetypes.guess_type(filename)[0]
    
    if separator:
        return haveTSV(sbtabfile,separator)
    elif file_mimetype == 'application/vnd.ms-excel' or file_mimetype == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' or file_mimetype == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
        return haveXLS(sbtabfile, True, True)        
    else:
        separator = misc.getDelimiter(sbtabfile)
        return haveTSV(sbtabfile, separator)            

    '''
    if file_mimetype == 'text/csv':
        return haveTSV(sbtabfile,'c')
    elif file_mimetype == 'text/tab-separated-values':
        return haveTSV(sbtabfile, 't')
    elif file_mimetype == 'text/tsv':
        return haveTSV(sbtabfile, 't')    
    elif file_mimetype == 'application/vnd.ms-excel' or file_mimetype == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
        return haveXLS(sbtabfile, True, True)
    else:
        return None
        #raise TypeError("%s is not in a supported format" % fpath)
    '''

def importSet(fpath):
    if not os.path.isfile(fpath):
        raise NameError("%s is not a valid filename" % fpath)
    mimetypes.init()
    file_mimetype = mimetypes.guess_type(fpath)[0]
    if file_mimetype == 'text/csv':
        return loadCSV(fpath, False)  # True says it has headers
    elif file_mimetype == 'text/tab-separated-values':
        return loadTSV(fpath, False)
    elif file_mimetype == 'text/tsv':
        return loadTSV(fpath, False)    
    elif file_mimetype == 'application/vnd.ms-excel' or file_mimetype == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
        return loadXLS(fpath, False, True)
    else:
        raise TypeError("%s is not in a supported format" % fpath)


def importBook(fpath):
    if not os.path.isfile(fpath):
        raise NameError("%s is not a valid filename" % fpath)
    mimetypes.init()
    file_mimetype = mimetypes.guess_type(fpath)[0]
    if file_mimetype == 'application/vnd.oasis.opendocument.spreadsheet':
        return loadODS(fpath, True, False)  # Second flag is for set import only
    elif file_mimetype == 'application/vnd.ms-excel' or file_mimetype == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
        return loadXLS(fpath, True, False)
    else:
        raise TypeError("%s is not in a supported format" % fpath)


def loadCSV(fpath, headers):
    csvfile = open(fpath, 'r')
    in_stream = csvfile.read()
    csvfile.close()
    dset = tablib.Dataset()
    rows = csv.reader(in_stream.splitlines(), quotechar='"')
    for i, row in enumerate(rows):
        if (i == 0) and (headers):
            dset.headers = row
        else:
            dset.append(row)
    return dset


def haveCSV(csvfile,headers):
    '''
    not needed anymore?
    '''
    #csvfile = open(fpath, 'r')
    in_stream = csvfile
    #csvfile.close()
    dset = tablib.Dataset()
    rows = csv.reader(in_stream.splitlines(), delimiter=',',quotechar='"')
    try:
        longest = max([len(x) for x in rows])
    except:
        raise TypeError("File is empty!")
    for i, row in enumerate(rows):
        # Skip empty rows
        if not row:
            continue
        if len(row) < longest:
            for i in range(longest - len(row)):
                row.append('')  # ro
        if (i == 0) and (headers):
            dset.headers = row
        else:
            dset.append(row)

    return dset

def haveTSV(tsvfile,separator):
    
    in_stream = tsvfile     #.read()

    dset = tablib.Dataset()
    rows = list(csv.reader(in_stream.splitlines(), delimiter=separator, quotechar='"'))    

    try:
        longest = max([len(x) for x in rows])
    except:
        raise TypeError("File is empty!")
    for i, row in enumerate(rows):
        # Skip empty rows
        if not row:
            continue
        if len(row) < longest:
            for i in range(longest - len(row)):
                row.append('')  # ro
        #if (i == 0) and (headers):
        #    dset.headers = row
        #else:
        dset.append(row)

    return dset 

def loadTSV(fpath, headers):
    tsvfile = open(fpath, 'r')
    in_stream = tsvfile.read()
    tsvfile.close()
    dset = tablib.Dataset()
    rows = list(csv.reader(in_stream.splitlines(), delimiter='\t', quotechar='"'))
    try:
        longest = max([len(x) for x in rows])
    except:
        raise TypeError("File is empty!")
    for i, row in enumerate(rows):
        # Skip empty rows
        if not row:
            continue
        if len(row) < longest:
            for i in range(longest - len(row)):
                row.append('')  # ro
        if (i == 0) and (headers):
            dset.headers = row
        else:
            dset.append(row)
    return dset

def haveXLS(file,headers,set_only):

    dbook = tablib.Databook()
    xl = xlrd.open_workbook(file_contents=file)

    for sheetname in xl.sheet_names():
        dset = tablib.Dataset()
        dset.title = sheetname
        sheet = xl.sheet_by_name(sheetname)
        for row in range(sheet.nrows):
            if (row == 0) and (headers):
                dset.headers = sheet.row_values(row)
            dset.append(sheet.row_values(row))
        dbook.add_sheet(dset)

    if set_only:
        return dbook.sheets()[0]
    else:
        return dbook   

def loadXLS(fpath, headers, set_only):
    dbook = tablib.Databook()
    f = open(fpath, 'rb')
    xl = xlrd.open_workbook(file_contents=f.read())
    for sheetname in xl.sheet_names():
        dset = tablib.Dataset()
        dset.title = sheetname
        sheet = xl.sheet_by_name(sheetname)
        for row in range(sheet.nrows):
            if (row == 0) and (headers):
                dset.headers = sheet.row_values(row)
            else:
                dset.append(sheet.row_values(row))
        dbook.add_sheet(dset)
    if set_only:
        return dbook.sheets()[0]
    else:
        return dbook


def writeCSV(data, fpath):
    outputfile = open(fpath + '.csv', 'wb')
    print('c')
    try:
        outputfile.write(data.csv)
        outputfile.close()
    except:
        print("It's a Databook! Writing single csv's for each sheet!")
        for sheet in data.sheets():
            outputfile = open(fpath + '_' + sheet.title + '.csv', 'w')
            outputfile.write(sheet.csv)
            outputfile.close()

def writeTSV(data, fpath):
    outputfile = open(fpath + '.tsv', 'wb')
    print('t')
    try:
        outputfile.write(data.tsv)
        outputfile.close()
    except:
        print("It's a Databook! Writing single tsv's for each sheet!")
        for sheet in data.sheets():
            outputfile = open(fpath + '_' + sheet.title + '.tsv', 'w')
            outputfile.write(sheet.tsv)
            outputfile.close()

def writeXLS(data, fpath):
    outputfile = open(fpath + '.xls', 'wb')
    outputfile.write(data.xls)
    outputfile.close()


def writeXLSX(data, fpath):
    outputfile = open(fpath + '.xlsx', 'wb')
    outputfile.write(data.xlsx)
    outputfile.close()


def writeODS(data, fpath):
    outputfile = open(fpath + '.ods', 'wb')
    outputfile.write(data.ods)
    outputfile.close()

'''example
    import tablib
    from tablibIO import *
    dataset = importSet('test.csv') # or TSV, ODS, XLS, XLSX
    dataset.headers  # will show the headers
    dataset.get_col(0)  # gets first column
    dataset.csv  # a csv representation of the table
    writeCSV(dataset, 'myoutput')  # writes your table to 'myoutput.cvs'
    writeXLS(dataset, 'myoutput')  # writes your table to 'myoutput.xls'
'''

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise NameError("no filename entered")
    elif not os.path.isfile(sys.argv[1]):
        raise NameError("%s is not a valid filename" % sys.argv[1])
    filepath = sys.argv[1]
    data = importBook(filepath)
    writeTSV(data, 'output')
