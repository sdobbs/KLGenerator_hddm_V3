

import os
import sys
import subprocess
import glob

# Get HALLD_HOME environment variable, verifying it is set
halld_home = os.getenv('HALLD_RECON_HOME')
if(halld_home == None):
        print('HALLD_RECON_HOME environment variable not set!')
        exit(-1)

# Get HALLD_MY if it exists. Otherwise use HALLD_HOME
halld_my = os.getenv('HALLD_MY', halld_home)

# Get executable name
executablename = os.path.split(os.getcwd())[1]

# Add SBMS directory to PYTHONPATH
sbmsdir = "%s/src/SBMS" % halld_home
sys.path.append(sbmsdir)

import sbms

# Get command-line options
SHOWBUILD = ARGUMENTS.get('SHOWBUILD', 0)

# Get platform-specific name
osname = os.getenv('BMS_OSNAME', 'build')

# Get architecture name
arch = subprocess.Popen(["uname"], stdout=subprocess.PIPE).communicate()[0].strip()

# Setup initial environment
executableinstalldir = halld_my
gluexinstalldir = "%s/%s" %(halld_home, osname)
include = "%s/include" % (gluexinstalldir)
bin = "%s/bin" % (executableinstalldir)
lib = "%s/lib" % (executableinstalldir)
executables = "%s/%s/executables" % (executableinstalldir,osname)
env = Environment(        ENV = os.environ,  # Bring in full environement, including PATH
                      CPPPATH = [include],
                      LIBPATH = ["%s/%s/lib" %(halld_home, osname)],  # n.b. add HALLD_HOME here and prepend HALLD_MY below
                  variant_dir = ".%s" % (osname))

# Only add HALLD_MY library search path if it already exists
# since we'll get a warning otherwise
if (os.path.exists(lib)): env.PrependUnique(lib)

# These are SBMS-specific variables (i.e. not default scons ones)
env.Replace(INSTALLDIR    = executableinstalldir,
                                OSNAME        = osname,
                                INCDIR        = include,
                                BINDIR        = bin,
                                LIBDIR        = lib,
#                                PLUGINSDIR    = plugins,
                                ALL_SOURCES   = [],        # used so we can add generated sources
                                SHOWBUILD     = SHOWBUILD,
                 COMMAND_LINE_TARGETS = COMMAND_LINE_TARGETS)
# Use terse output unless otherwise specified
if SHOWBUILD==0:
        env.Replace(   CCCOMSTR       = "Compiling  [$SOURCE]",
                                  CXXCOMSTR       = "Compiling  [$SOURCE]",
                                  FORTRANPPCOMSTR = "Compiling  [$SOURCE]",
                                  FORTRANCOMSTR   = "Compiling  [$SOURCE]",
                                  SHCCCOMSTR      = "Compiling  [$SOURCE]",
                                  SHCXXCOMSTR     = "Compiling  [$SOURCE]",
                                  LINKCOMSTR      = "Linking    [$TARGET]",
                                  SHLINKCOMSTR    = "Linking    [$TARGET]",
                                  INSTALLSTR      = "Installing [$TARGET]",
                                  ARCOMSTR        = "Archiving  [$TARGET]",
                                  RANLIBCOMSTR    = "Ranlib     [$TARGET]")

# Get compiler from environment variables (if set)
env.Replace( CXX = os.getenv('CXX', 'g++'),
             CC  = os.getenv('CC' , 'gcc'),
             FC  = os.getenv('FC' , 'gfortran') )

# Use C++11
env.PrependUnique(    CXXFLAGS = ['-std=c++11'])

# Add local directory, directories from HALLD_MY and HALLD_HOME to include search path
#env.PrependUnique(CPPPATH = ['#'])
env.PrependUnique(CPPPATH = ['%s/src' % halld_my, '%s/src/libraries' % halld_my, '%s/src/libraries/include' % halld_my])
env.PrependUnique(CPPPATH = ['%s/src' % halld_home, '%s/src/libraries' % halld_home, '%s/src/libraries/include' % halld_home])

# Turn on debug symbols and warnings
env.PrependUnique(      CFLAGS = ['-g', '-fPIC', '-Wall'])
env.PrependUnique(    CXXFLAGS = ['-g', '-fPIC', '-Wall'])
env.PrependUnique(FORTRANFLAGS = ['-g', '-fPIC', '-Wall'])

# Apply any platform/architecture specific settings
#sbms.ApplyPlatformSpecificSettings(env, arch)
#sbms.ApplyPlatformSpecificSettings(env, osname)


# Make executable from source in this directory
sbms.AddDANA(env)
sbms.AddROOT(env)
sbms.executable(env,".%s/%s" % (osname,executablename))

# Make install target
env.Alias('install', env.Install(bin, ".%s/%s" % (osname,executablename)))

