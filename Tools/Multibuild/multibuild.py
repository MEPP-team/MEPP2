
import subprocess
import sys
import time
import datetime
import collections
import os
import shutil
import re
import zlib
import ConfigParser
#DBG import pdb



BuildConfig = collections.namedtuple('BuildConfig', 'buildtype, use_cgal, use_openmesh, use_aif, use_gui, use_qt5, use_pcl, use_vtk, build_examples')

configurations = []
#
# define build & test configurations
#
# all datastructures, Debug then Release
configurations.append(BuildConfig(buildtype='Debug', use_cgal='ON', use_openmesh='ON', use_aif='ON', use_gui='ON', use_qt5='OFF', use_pcl='OFF', use_vtk='OFF', build_examples='ON'))
configurations.append(BuildConfig(buildtype='Release', use_cgal='ON', use_openmesh='ON', use_aif='ON', use_gui='ON', use_qt5='OFF', use_pcl='OFF', use_vtk='OFF', build_examples='ON'))
# Openmesh only, Debug then release
configurations.append(BuildConfig(buildtype='Debug', use_cgal='OFF', use_openmesh='ON', use_aif='OFF', use_gui='ON', use_qt5='OFF', use_pcl='OFF', use_vtk='OFF', build_examples='ON'))
configurations.append(BuildConfig(buildtype='Release', use_cgal='OFF', use_openmesh='ON', use_aif='OFF', use_gui='ON', use_qt5='OFF', use_pcl='OFF', use_vtk='OFF', build_examples='ON'))
# AIF only, Debug then release
configurations.append(BuildConfig(buildtype='Debug', use_cgal='OFF', use_openmesh='OFF', use_aif='ON', use_gui='ON', use_qt5='OFF', use_pcl='OFF', use_vtk='OFF', build_examples='ON'))
configurations.append(BuildConfig(buildtype='Release', use_cgal='OFF', use_openmesh='OFF', use_aif='ON', use_gui='ON', use_qt5='OFF', use_pcl='OFF', use_vtk='OFF', build_examples='ON'))
# CGAL only, Debug then Release
configurations.append(BuildConfig(buildtype='Debug', use_cgal='ON', use_openmesh='OFF', use_aif='OFF', use_gui='ON', use_qt5='OFF', use_pcl='OFF', use_vtk='OFF', build_examples='ON'))
configurations.append(BuildConfig(buildtype='Release', use_cgal='ON', use_openmesh='OFF', use_aif='OFF', use_gui='ON', use_qt5='OFF', use_pcl='OFF', use_vtk='OFF', build_examples='ON'))

# global variables
CGAL_DIR = ''
OPENMESH_DIR = ''
GIT_PATH = ''
KIT_ROOT = ''
MSBUILD_PATH = ''
CMAKE_PATH = ''
CTEST_PATH = ''

# load parameters
def init(argv):
	global CGAL_DIR
	global OPENMESH_DIR
	global GIT_PATH
	global KIT_ROOT
	global MSBUILD_PATH
	global CMAKE_PATH
	global CTEST_PATH

	# test command line arguments
	if len(argv) != 2:
		print('Launch compiling and testing of a pre-defined set of MEPP2 configurations.')
		print('Usage:  ' + argv[0] + '  path_to_multibuild.conf')
		print('Example:  ' + argv[0] + '  Tools/Multibuild/multibuild.conf')
		sys.exit(1)

	# read configuration file
	check_path(argv[1], 'configuration file')
	config = ConfigParser.RawConfigParser()
	config.read(argv[1])

	if sys.platform.startswith('linux'):
		# linux
		print('Linux detected')
		CGAL_DIR = config.get('Linux', 'CGAL_DIR').strip('"')
		OPENMESH_DIR = config.get('Linux', 'OPENMESH_DIR').strip('"')
	elif sys.platform == "win32":
		# Windows
		print('Windows detected')
		#DBG pdb.set_trace() #breakpoint
		GIT_PATH = config.get('Windows', 'GIT_PATH').strip('"')
		check_path(GIT_PATH, 'git executable')
		KIT_ROOT = config.get('Windows', 'KIT_ROOT').strip('"')
		check_path(KIT_ROOT, 'KIT_ROOT dir')
		MSBUILD_PATH = config.get('Windows', 'MSBUILD_PATH').strip('"')
		check_path(MSBUILD_PATH, 'msbuild executable')
		CMAKE_PATH = KIT_ROOT + '\\_utils_\\cmake-3.4.3-win32-x86\\bin\\cmake.exe'
		check_path(CMAKE_PATH , 'cmake executable')
		CTEST_PATH = KIT_ROOT + '\\_utils_\\cmake-3.4.3-win32-x86\\bin\\ctest.exe'
		check_path(CTEST_PATH , 'ctest executable')
	elif sys.platform == "darwin":
		# OS X
		print('Mac OSX detected, not yet supported')
		sys.exit(1)


def check_path(path, message):
	"""Check that path exists, terminate if the path is invalid"""
	if not os.path.exists(path):
		print('Could not find ' + message + ': ' + path)
		sys.exit(1)


def run_command(title, command, wdir, logfile):
	"""Run a command, write output to log file"""
	"""Return True if the command was successfull"""
	logfile.write('\n')
	logfile.write('*************\n')
	logfile.write('RUNNING ' + title + '\n')
	logfile.write(time.asctime() + '\n')
	logfile.write('*************\n')
	logfile.write('\n')
	logfile.write('directory:\n')
	logfile.write('  ' + wdir + '\n')
	logfile.write('command:\n')
	logfile.write('  ' + command + '\n')
	logfile.write('\n')
	logfile.flush()

	try:
		if sys.platform.startswith('linux'):
			subprocess.call(command, stdout=logfile, stderr=logfile, cwd=wdir, shell=True)
		elif sys.platform == "win32":
			subprocess.call(command, stdout=logfile, stderr=logfile, cwd=wdir)
		elif sys.platform == "darwin":
			subprocess.call(command, stdout=logfile, stderr=logfile, cwd=wdir, shell=True)

		logfile.write('\n')
		logfile.write('**********************\n')
		logfile.write(title + ' RUN SUCCESSFULLY\n')
		logfile.write(time.asctime() + '\n')
		logfile.write('**********************\n')
		logfile.flush()
		return True #success

	except subprocess.CalledProcessError as e:
		logfile.write('the command:\n')
		logfile.write('\n')
		logfile.write(e.cmd + '\n')
		logfile.write('\n')
		logfile.write('failed with output:\n')
		logfile.write('\n')
		logfile.write(e.output + '\n')
		logfile.write('\n')
		logfile.write('***************\n')
		logfile.write(title + ' FAILED!!!\n')
		logfile.write(time.asctime() + '\n')
		logfile.write('***************\n')
		logfile.flush()
		return False #success


def build_one_config(config, builddir, logfilename):
	"""config is a BuildConfig object"""
	start_time = time.time()

	# init log file
	logfile = open(logfilename, 'a')

	# write configuration to logfile
	logfile.write(str(config))
	logfile.write('\n')

	# run cmake
	if sys.platform.startswith('linux'):
		cmakecmd = ''
		if CGAL_DIR:
			cmakecmd += ' CGAL_DIR=' + CGAL_DIR
		if OPENMESH_DIR:
			cmakecmd += ' OPENMESH_DIR=' + OPENMESH_DIR
		cmakecmd += ' cmake'
		cmakecmd += ' -DCMAKE_BUILD_TYPE=' + config.buildtype
	elif sys.platform == "win32":
		cmakecmd = CMAKE_PATH
		cmakecmd += ' -G"Visual Studio 14 2015 Win64"'
		cmakecmd += ' -DMSVC_KIT_ROOT="' + KIT_ROOT + '"'
	elif sys.platform == "darwin":
		pass
		#TODO

	cmakecmd += ' -DBUILD_USE_CGAL=' + config.use_cgal
	cmakecmd += ' -DBUILD_USE_OPENMESH=' + config.use_openmesh
	cmakecmd += ' -DBUILD_USE_AIF=' + config.use_aif
	cmakecmd += ' -DBUILD_USE_GUI=' + config.use_gui
	cmakecmd += ' -DBUILD_USE_QT5=' + config.use_qt5
	cmakecmd += ' -DBUILD_USE_PCL=' + config.use_pcl
	cmakecmd += ' -DBUILD_USE_VTK=' + config.use_vtk
	cmakecmd += ' -DBUILD_EXAMPLES=' + config.build_examples
	#cmakecmd += ' -DBUILD_ADVANCED_TESTS_COMPR_VALENCE=ON'
	cmakecmd += ' ..'

	success = run_command('CMAKE', cmakecmd, builddir, logfile)

	# run compiler
	if success:
		if sys.platform.startswith('linux'):
			buildcmd = 'make -j2'
		elif sys.platform == "win32":
			buildcmd = MSBUILD_PATH
			buildcmd += ' /m'
			buildcmd += ' /property:Configuration=' + config.buildtype
			buildcmd += ' MEPP2.sln'
			# note: "/maxcpucount" option make msbuild fail with some sort of syntax error
		elif sys.platform == "darwin":
			pass
			#TODO

		success = run_command('BUILD', buildcmd, builddir, logfile)

	# run ctest
	if success:
		if sys.platform.startswith('linux'):
			ctestcmd = 'ctest  -j 2  --output-on-failure'
			ctestcmd += ' -E mepp-gui'
		elif sys.platform == "win32":
			ctestcmd = CTEST_PATH + '  -j 2  --output-on-failure'
		elif sys.platform == "darwin":
			pass
			#TODO

		success = run_command('CTEST', ctestcmd, builddir, logfile)

	# display build time
	build_time = time.time() - start_time
	logfile.write('\n')
	logfile.write('total build time: ' + str(datetime.timedelta(seconds=int(build_time))))
	logfile.write('\n')
	logfile.close()


def update_report(report, buildtimestamp, buildstatus, outputfilename):
	file = open(outputfilename, 'w')
	file.write('<!DOCTYPE html>\n')
	file.write('<html>\n')
	file.write('  <head>\n')
	file.write('    <meta charset="utf-8" />\n')
	file.write('    <meta http-equiv="refresh" content="30">\n')
	file.write('    <title>Multibuild report</title>\n')
	file.write('  </head>\n')
	file.write('\n')
	file.write('  <body>\n')
	file.write('    <h1 style="text-align:center">Build report ' + buildtimestamp + '</h1>\n')
	file.write('    <h2 style="text-align:center;color:' + buildstatus[0] + '">' + buildstatus[1] + '</h2>\n')

	for rep in report:
			file.write('    <hr>\n')
			file.write('    <h2>build ' + str(rep['configid']) + '</h2>\n')
			file.write('    ' + rep['configstr'] + '<br>\n')
			file.write('    started at  ' + rep['starttime'] + '<br>\n')
			file.write('    <a href="file:///' + os.getcwd() + '/' + rep['logfilename'] + '">log file</a>\n')
			if rep['durationstr']:
				color = 'green'
				if int(rep['warningsnbr']) > 0:
					color = 'orange'
				file.write('&emsp;<span style="color:' + color + '">' + rep['warningsstr'] + '</span>\n')
				color = 'green'
				if int(rep['errorsnbr']) > 0:
					color = 'red'
				file.write('&emsp;<span style="color:' + color + '">' + rep['errorsstr'] + '</span>\n')
				color = 'green'
				if int(rep['testssnbr']) < 100:
					color = 'red'
				file.write('&emsp;<span style="color:' + color + '">' + rep['testsstr'] + '</span><br>\n')
				file.write('    ' + rep['durationstr'] + '<br>\n')

	file.write('  </body>\n')
	file.write('</html>\n')
	file.close()


def display_file_in_browser(filename):
	if sys.platform.startswith('linux'):
		cmd = ['xdg-open', filename]
		subprocess.Popen(cmd)
	elif sys.platform == "win32":
		cmd = ['start', filename]
		subprocess.Popen(cmd, shell=True)
	elif sys.platform == "darwin":
		cmd = ['open', filename]
		subprocess.Popen(cmd)


def get_commit_info():
	if sys.platform.startswith('linux'):
		gitcmd = ['git', 'log', '-1']
	elif sys.platform == "win32":
		gitcmd = [GIT_PATH, 'log', '-1']
	elif sys.platform == "darwin":
		gitcmd = ['git', 'log', '-1']

	#DBG pdb.set_trace() #breakpoint
	commit_info = subprocess.check_output(gitcmd)
	return commit_info


def build_all_configs():
	commit_info = get_commit_info()
	buildtimestamp = time.strftime("%Y%m%d-%H%M%S")
	buildname = 'build.' + buildtimestamp
	os.mkdir(buildname)
	reportfilename = buildname + '/multibuild_report.html'
	build_compile_error_flag = False
	build_compile_warning_flag = False
	build_test_error_flag = False
	buildstatus = ('blue', 'in progress...')
	#dbg print('report in ' + reportfilename)
	report = [] # a list of dictionnaries
	update_report(report, buildtimestamp, buildstatus, reportfilename)
	display_file_in_browser(reportfilename)

	configid = 0
	for c in configurations:
		configid += 1
		configname = 'build_' + str(configid)
		logfilename = buildname + '/' + configname + '.log'
		print('start ' + configname)
		# init report info
		rep = {}
		rep['configid'] = configid
		rep['configname'] = configname
		rep['logfilename'] = logfilename
		rep['configstr'] = str(c)
		rep['starttime'] = time.strftime("%Y-%m-%d  %H:%M:%S")
		rep['durationstr'] = None
		report.append(rep)
		update_report(report, buildtimestamp, buildstatus, reportfilename)
		# clean or create build dir
		builddir = 'build.' + str(configid)
		shutil.rmtree(builddir, ignore_errors=True)
		os.mkdir(builddir)
		# create log file and write commit info
		with open(logfilename, 'w') as f:
			f.write(commit_info)
			f.write('\n')
		# do build
		build_one_config(c, builddir, logfilename)
		print('completed ' + configname)
		# get build summary
		summary = ''
		logfile = open(logfilename, 'r')
		for line in logfile:
			if sys.platform.startswith('linux'):
				if re.search(r'BuildConfig|warning:|error:|tests failed out of|total build time', line):
					summary += line
			elif sys.platform == "win32":
				if re.search(r'BuildConfig|Avertissement\(s\)|Warning\(s\)|Erreur\(s\)|Error\(s\)|tests failed out of|total build time', line):
					summary += line
			elif sys.platform == "darwin":
				if re.search(r'BuildConfig|warning:|error:|tests failed out of|total build time', line):
					summary += line
		logfile.close()
		# store summary in report
		report[-1]['configstr'] = re.search(r'BuildConfig.*', summary).group(0)

		if sys.platform.startswith('linux'):
			try:
				report[-1]['warningsnbr'] = str(len(re.findall('warning:', summary)))
				report[-1]['warningsstr'] = report[-1]['warningsnbr'] + ' warning(s)'
			except:
				report[-1]['warningsstr'] = '999 warnings info not available, look at the log file !!!'

			try:
				report[-1]['errorsnbr'] = str(len(re.findall('error:', summary)))
				report[-1]['errorsstr'] = report[-1]['errorsnbr'] + ' error(s)'
			except:
				report[-1]['errorsstr'] = '999 errors info not available, look at the log file !!!'

			try:
				report[-1]['testsstr'] = re.search(r'[0-9]+% tests passed.*', summary).group(0)
			except:
				report[-1]['testsstr'] = '0% tests info not available, look at the log file !!!'

			report[-1]['durationstr'] = re.search(r'total build time.*', summary).group(0)
			report[-1]['testssnbr'] = re.search(r'^[0-9]+', rep['testsstr']).group(0)
			report[-1]['duration'] = re.search(r'[0-9]+:[0-9]+:[0-9]+$', rep['durationstr']).group(0)
		elif sys.platform == "win32":
			try:
				report[-1]['warningsstr'] = re.findall(r'[0-9]+ Avertissement\(s\)', summary)
				# returns a list
			except:
				try:
					report[-1]['warningsstr'] = re.findall(r'[0-9]+ Warning\(s\)', summary)
					# returns a list
				except:
					report[-1]['warningsstr'] = ['999 warnings info not available, look at the log file !!!']

			try:
				report[-1]['errorsstr'] = re.findall(r'[0-9]+ Erreur\(s\)', summary)
				# returns a list
			except:
				try:
					report[-1]['errorsstr'] = re.findall(r'[0-9]+ Error\(s\)', summary)
					# returns a list
				except:
					report[-1]['errorsstr'] = ['999 errors info not available, look at the log file !!!']

			try:
				report[-1]['testsstr'] = re.search(r'[0-9]+% tests passed.*', summary).group(0)
			except:
				report[-1]['testsstr'] = '0% tests info not available, look at the log file !!!'

			report[-1]['durationstr'] = re.search(r'total build time.*', summary).group(0)

			report[-1]['warningsnbr'] = 0
			for warnstr in rep['warningsstr']:
				warnnbr = int(re.search(r'^[0-9]+', warnstr).group(0))
				report[-1]['warningsnbr'] += warnnbr
			report[-1]['warningsnbr'] = str(report[-1]['warningsnbr'])
			report[-1]['warningsstr'] = str(report[-1]['warningsstr'])

			report[-1]['errorsnbr'] = 0
			for errstr in rep['errorsstr']:
				errnbr = int(re.search(r'^[0-9]+', errstr).group(0))
				report[-1]['errorsnbr'] += errnbr
			report[-1]['errorsnbr'] = str(report[-1]['errorsnbr'])
			report[-1]['errorsstr'] = str(report[-1]['errorsstr'])

			report[-1]['testssnbr'] = re.search(r'^[0-9]+', rep['testsstr']).group(0)
			report[-1]['duration'] = re.search(r'[0-9]+:[0-9]+:[0-9]+$', rep['durationstr']).group(0)
		elif sys.platform == "darwin":
			pass
			#TODO

		update_report(report, buildtimestamp, buildstatus, reportfilename)

		if int(report[-1]['errorsnbr']) > 0:
			build_compile_error_flag = True
		if int(report[-1]['warningsnbr']) > 0:
			build_compile_warning_flag = True
		if int(report[-1]['testssnbr']) < 100:
			build_test_error_flag = True

		# fix spaces at the beginning of the lines
		#summary = re.sub(r'(BuildConfig)', r'    \1', summary)
		#summary = re.sub(r'([0-9]+% tests passed)', r'    \1', summary)
		#summary = re.sub(r'(total build time)', r'    \1', summary)
		#print(configname + ' summary:')
		#print(summary)
		#sys.stdout.flush()

	# update global build status
	if build_compile_error_flag:
		buildstatus = ('red', 'failed with compilation error(s)')
	elif build_test_error_flag:
		buildstatus = ('red', 'failed with test error(s)')
	elif build_compile_warning_flag:
		buildstatus = ('orange', 'passed with compilation warning(s)')
	else:
		buildstatus = ('green', 'passed')
	update_report(report, buildtimestamp, buildstatus, reportfilename)


if __name__ == '__main__':
	init(sys.argv)
	#build_one_config(config='Debug', use_cgal='OFF', use_aif='OFF')
	build_all_configs()
