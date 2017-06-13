#!/usr/bin/python
###############################################################################
## Arplib to Origen Binary Library (.obl) converter
## Nicholas C Sly
## Mar 10, 2015
###############################################################################
## Designed to take a .arplib file type and convert it to the Origen Binary
## Library (.obl) file type complete with appropriate ID and Interpolable Tags.
## Initial intention is for the path to the .arplib files to be hard-coded
## and to use the arpdata.txt file to tag as appropriate.  Conversion will use
## the obconvert application to convert the library type and tagging will use
## the obtagmod application, both of which can be found in 
## ${SCALE}/src/Origen/Tools/.
##
## This script will be stealing heavily from the mpihpcshiftrun.py script
## on OIC, but is coded for use on jupiter.
###############################################################################

import os, sys, fnmatch, shutil, re, string, time
from subprocess import *
#from datetime import date

##---------------------------------------------------------------------------##
## PRESET VARIABLES
##---------------------------------------------------------------------------##

## Change if location of obiwan executable is different.
## Alternatively, use the --obiwan option to specify the location
exec_line = 'obiwan '

##---------------------------------------------------------------------------##
## Usage info if requested.
##---------------------------------------------------------------------------##

if ('--help' in sys.argv):
    print ''
    print 'arp2obl.py -- Script for the transformation of ORIGEN libraries to'
    print '              the .obl filetype while using the information available'
    print '              in the arpdata.txt file to add relevant id tag and'
    print '              interp tag information for use with obiwan interpolation.'
    print ''
    print 'Usage: /usr/bin/python arp2obl.py [--arpdata /path/to/arpdata.txt]'
    print '       [--libdir /path/to/origen/library/directory/]'
    print '       [--obiwan /path/to/obiwan/executable]'
    print '       [--help]'
    print ''
    print 'Note: If no --arpdata option is used, the program will search for the'
    print '      arpdata.txt file in the current working directory.'
    print ''
    print 'Note: This script produces significant output.  If not desired,'
    print '      redirect output by appending &>/dev/null to the end of the'
    print '      execution call.'
    print ''
    sys.exit(0)

##---------------------------------------------------------------------------##
## Checking system arguments for non-default values.
##---------------------------------------------------------------------------##

libdir = '.'

if ('--arpdata' in sys.argv):
    arpdata_loc = sys.argv[sys.argv.index('--arpdata') + 1]

if ('--libdir' in sys.argv):
    libdir = sys.argv[sys.argv.index('--libdir') + 1]
print 'ORIGEN Library directory used: ', libdir

if ('--obiwan' in sys.argv):
    exec_line = sys.argv[sys.argv.index('--obiwan') + 1]
print 'Obiwan executable location used: ', exec_line

##---------------------------------------------------------------------------##
## Parsing arpdata.txt
##---------------------------------------------------------------------------##

if 'arpdata_loc' in globals():
    if 'arpdata.txt' not in arpdata_loc:
        if ('/'!=arpdata_loc[len(arpdata_loc)-1]): arpdata_loc+='/'
        arpdata_loc+='arpdata.txt'
else: arpdata_loc=os.getcwd() + '/arpdata.txt'
print 'Arpdata.txt location used: ', arpdata_loc

try:
    arpdata = open(arpdata_loc)
except IOError:
    print 'Could not locate arpdata.txt file.  Try using `arp2obl.py --arpdata /path/to/arpdata.txt`'
num_files = 0
burn_list = []
num_burnsteps = ''
for line in arpdata:
    ## Each section of the arpdata.txt file starts with an '!' followed by the
    ## name for the assembly type.
    if (re.search("^!",line)):
        line_list=line.split('!')
        name=line_list[1].rstrip()
        print "Found Assembly Type: {0}".format(name)
        line_num=1
        file_names = []
        burnup_list = []
        type = 0
    ## The next line of the section tells you how many items to expect in the 
    ## following sections.
    elif line_num == 2:
        line_tmp=line.split('\n')[0].lstrip(' ')
        line_list=line_tmp.split(' ')
        i=0
        while i<len(line_list):
            if line_list[i] == '':
                del line_list[i]
            i+=1
        ## If this line contains 3 numbers, it means that the fuel is Uranium
        ## based.  The numbers tell you the number of enrichments and the number
        ## of moderator densities for which libraries have been created.  The last
        ## number tells you how many burnups are on each library.  The number of
        ## files to expect can be calculated as the number of enrichments multiplied
        ## by the number of moderator densities.
        if len(line_list) == 3:
            num_enrich = line_list[0]
            num_modden = line_list[1]
            num_files = int(num_enrich)*int(num_modden)
            num_burnsteps = line_list[2]
            print "These libraries are uranium fueled.  There are {0} libraries\
                   for the {1} assembly type.".format(num_files,name)
            type = 1
        ## If this line contains 5 numbers, it means that the fuel is MOX based.
        ## The first number tells you the number of plutonium levels characterized.
        ## The second number tells you how much of that plutonium is Pu-239.
        ## The third number is not used and is always 1.
        ## The fourth number tells you the moderator densities and the last number
        ## tells you the number of burnup steps on each library.  The number of
        ## files to expect can be calculated as 
        ## (the number of plutonium levels)*(Pu-239 content)*(moderator densities).
        elif len(line_list) == 5:
            num_pu_content = line_list[0]
            num_pu239_vals = line_list[1]
            unused_var = line_list[2]
            num_modden = line_list[3]
            num_files = int(num_pu_content)*int(num_pu239_vals)*int(num_modden)
            num_burnsteps = line_list[4]
            print "These libraries are MOX fueled.  There are {0} libraries for\
                   the {1} assembly type.".format(num_files,name)
            print "The number of Plutonium contents is {0}.".format(num_pu_content)
            print "The number of Pu-239 values is {0}.".format(num_pu239_vals)
            print "The number of moderator densities is {0}.".format(num_modden)
            type = 2
    ## For uranium fuel libraries, this line contains the enrichments for which
    ## libraries have been generated.
    elif type == 1 and line_num == 3:
        line_tmp=line.split('\n')[0]
        line_list=line_tmp.split(' ')
        enrich = []
        for item in line_list:
            if item != '': enrich.append(item)
    ## For MOX fuel libraries, this line contains the percentages of the fuel
    ## that are plutonium.
    elif type == 2 and line_num == 3:
        line_tmp=line.split('\n')[0]
        line_list=line_tmp.split(' ')
        pu_content = []
        for item in line_list:
            if item != '': pu_content.append(item)
        print 'Plutonium content is: {0}'.format(pu_content)
    ## For uranium fuel libraries, this line contains the moderator densities
    ## for which libraries have been generated.
    elif type == 1 and line_num == 4:
        line_tmp=line.split('\n')[0]
        line_list=line_tmp.split(' ')
        modden = []
        for item in line_list:
            if item != '': modden.append(item)
    ## For MOX fuel libraries, this line contains the percentages of the plutonium
    ## contents that are Pu-239.
    elif type == 2 and line_num == 4:
        line_tmp=line.split('\n')[0]
        line_list=line_tmp.split(' ')
        pu239_vals = []
        for item in line_list:
            if item != '': pu239_vals.append(item)
        print 'Pu-239 content is: {0}'.format(pu239_vals)
    ## For MOX fuel libraries, this number is not used and is always 1.0
    elif type == 2 and line_num == 5: 
        line_tmp=line.split('\n')[0]
        line_list=line_tmp.split(' ')
        if line_list[0] == 1:
            print "Check"
    ## For MOX fuel libraries, this line contains the moderator densities for which
    ## libraries have been generated.
    elif type == 2 and line_num == 6:
        line_tmp=line.split('\n')[0]
        line_list=line_tmp.split(' ')
        modden = []
        for item in line_list:
            if item != '': modden.append(item)
        print 'Pu moderator densities are: {0}'.format(modden)
    ## After the above sections, the arpdata.txt file lists the filenames
    ## of all of the libraries for this assembly type in a certain order.
    elif ( ( (type == 1 and line_num >= 4) or (type == 2 and line_num >= 7)) 
             and (len(file_names) != num_files)):
        line_tmp=line.split('\n')[0]
        line_list = line_tmp.split("'")
        for item in line_list:
            if len(item) > 5:
                print 'Name of library: {0}'.format(item)
                file_names.append(item)
        if (len(file_names) == num_files and type == 1):
            print 'Enrichments: {0}'.format(enrich)
            print 'Moderator Densities: {0}'.format(modden)
        if (len(file_names) == num_files and type == 2):
            print 'Plutonium Contents: {0}'.format(pu_content)
            print 'Plutonium-239 Contents: {0}'.format(pu239_vals)
            print 'Moderator Densities: {0}'.format(modden)
    ## Finally, after the filenames, the values for the burnups are given.
    ## These are (I believe) in MWD/MTHM.
    elif (len(file_names) == num_files):
        line_tmp=line.split('\n')[0]
        line_list=line_tmp.split(' ')
        for item in line_list:
            if item != '': burnup_list.append(item)
        print 'Length of burnup_list = {0} and num_burnsteps = {1}' \
              .format(len(burnup_list),num_burnsteps)
        ## With all of that information collected, it's time to call old Obiwan.
        if (len(burnup_list) == int(num_burnsteps)):
            print 'Entering execution section of the program.'
            burn_string = ''
            for burns in burnup_list:
                burn_string+=burns
                if burns != burnup_list[len(burnup_list) - 1]: burn_string += ' '
            if type == 1:
                for val in enrich:
                    for den in modden:
                        file_num = enrich.index(val)*len(modden)+modden.index(den)
                        os.system(("{0} tag -idtags='Assembly Type={1}," + \
                                   "Fuel Type=Uranium,Burnup Times={2}' " + \
                                   "-interptags='Enrichment={3}," + \
                                   "Moderator Density={4}' {5}") \
                                   .format(exec_line,name,burn_string,val,\
                                           den,file_names[file_num]))
                        print "Modified {0}".format(file_names[file_num])
            if type == 2:
                for den in modden:
                    for pu239 in pu239_vals:
                        for pu in pu_content:
                            file_num = (modden.index(den)*len(pu239_vals) + \
                                        pu239_vals.index(pu239))*len(pu_content) + \
                                        pu_content.index(pu)
                            os.system(("{0} tag -idtags='Assembly Type={1}," + \
                                       "Fuel Type=MOX,Burnup Times={2}' " + \
                                       "-interptags='Moderator Density={3}," + \
                                       "Plutonium Content={4}," + \
                                       "Plutonium-239 Content={5}' {6}") \
                                       .format(exec_line,name,burn_string,den,\
                                               pu,pu239,file_names[file_num]))
                            print "Modified {0}".format(file_names[file_num])
    line_num+=1

sys.exit(0)
