#! /usr/bin/env python
# encoding: utf-8
# Harald Klimach 2011
#
# Compilation will create the executable and execute all unit tests,
# if you only to compile, without execution, use --notests
import os
import glob

APPNAME = 'ateles'
VERSION = '1'

top = '.'
out = 'build'

def options(opt):
    opt.load('eclipse')
    opt.recurse('treelm')
    opt.add_option('--no_harvesting', action='store_true',
                   default=False,
                   help = 'Do not include harvesting files for compilation.',
                   dest = 'no_harvesting')

def configure(conf):
    from waflib import Logs
    if not conf.options.no_harvesting:
      conf.env.build_hvs = True
    else:
      conf.env.build_hvs = False
    conf.recurse('treelm')
    conf.setenv('')
    Logs.warn('Ateles specific configuration:')
    conf.recurse('polynomials')

    # Avoid some warnings in gfortran:
    if not conf.options.nowarn:
        for key, fenv in conf.all_envs.items():
            if fenv.FC_NAME == 'GFORTRAN':
                fenv.FCFLAGS.append('-Wno-unused-dummy-argument')
                fenv.FCFLAGS.append('-Wno-uninitialized')
                fenv.FCFLAGS.append('-fmax-stack-var-size=131072')

    Logs.warn('Ateles modified flags:')
    Logs.warn('Default flags: {0}'.format(' '.join(conf.all_envs[''].FCFLAGS)))
    Logs.warn('Debug   flags: {0}'.format(' '.join(conf.all_envs['debug'].FCFLAGS)))

    conf.setenv('')
    conf.setenv('ford', conf.env)
    conf.env.ford_mainpage = 'mainpage.md'

def build(bld):

    bld(rule='cp ${SRC} ${TGT}', source=bld.env.COCOSET, target='coco.set')

    bld.add_group()

    if bld.cmd != 'gendoxy':
        # Build Treelm and aotus recursive
        bld.recurse('treelm')

    else:
        bld(rule='cp ${SRC} ${TGT}',
            source = bld.path.find_node(['treelm', 'source', 'arrayMacros.inc']),
            target = bld.path.find_or_declare('arrayMacros.inc'))

    bld.recurse('polynomials')

    # Build the Ateles objects
    build_atl_objs(bld)



def build_atl_objs(bld):

    ateles_sources = bld.path.ant_glob('source/*.f90',
                                       excl='source/ateles.f90')
    ateles_sources += bld.path.ant_glob('source/equation/*.f90')
    ateles_sources += bld.path.ant_glob('source/filter/*.f90')
    ateles_sources += bld.path.ant_glob('source/flux/*.f90')
    ateles_sources += bld.path.ant_glob('source/scheme/*.f90')
    ateles_sources += bld.path.ant_glob('source/scheme/modg/*.f90')
    ateles_sources += bld.path.ant_glob('source/scheme/modg_2d/*.f90')
    ateles_sources += bld.path.ant_glob('source/scheme/modg_1d/*.f90')
    ateles_sources += bld.path.ant_glob('source/scheme/modg/fpt/*.f90')
    ateles_sources += bld.path.ant_glob('source/timestepping/*.f90')
    ateles_sources += bld.path.ant_glob('source/material/*.f90')
    ateles_sources += bld.path.ant_glob('source/exchange/*.f90')

    atl_ppsources = bld.path.ant_glob('source/*.fpp')
    atl_ppsources += bld.path.ant_glob('source/equation/*.fpp')
    atl_ppsources += bld.path.ant_glob('source/filter/*.fpp')
    atl_ppsources += bld.path.ant_glob('source/material/*.fpp')
    atl_ppsources += bld.path.ant_glob('source/scheme/*.fpp')
    atl_ppsources += bld.path.ant_glob('source/scheme/modg/*.fpp')
    atl_ppsources += bld.path.ant_glob('source/scheme/modg_2d/*.fpp')
    atl_ppsources += bld.path.ant_glob('source/scheme/modg_1d/*.fpp')

    ateles_sources += atl_ppsources

    if bld.cmd != 'gendoxy':

        compile_atl(bld, ateles_sources)
        bld(
            features = 'fc fcprogram',
            source = 'source/ateles.f90',
            use     = ['atl_objs', 'tem_objs', 'ply_objs', 'aotus', bld.env.distcrc,
                       bld.env.mpi_mem_c_obj, 'fftw_mod_obj', 'NAG',
                       'fxtp_wrap_obj', 'fxtp_obj', 'fxtp_wrapper',
                       'PRECICE','MPICXX', 'PYLIB', 'STDCXX', 'RT', 'ZLIB', 'PETSC',
                       'BOOST_system', 'BOOST_filesystem'],
            target = 'ateles')
        bld(
            features = 'fc fcprogram',
            source = 'peons/solve_euler_riemann.f90',
            use     = ['atl_objs', 'tem_objs', 'ply_objs', 'aotus',
                       bld.env.mpi_mem_c_obj, 'fftw_mod_obj', 'NAG', bld.env.distcrc,
                       'fxtp_wrap_obj', 'fxtp_obj', 'fxtp_wrapper',
                       'PRECICE', 'MPICXX', 'PYLIB', 'STDCXX', 'RT', 'ZLIB', 'PETSC', 
                       'BOOST_system', 'BOOST_filesystem'],
            target = 'solve_euler_riemann')
        if bld.env.build_hvs and not bld.options.no_harvesting:
            bld(
                features = 'fc fcprogram',
                source = 'source/atl_harvesting/atl_harvesting.f90',
                use     = ['atl_objs', 'tem_objs', 'ply_objs', 'aotus', bld.env.distcrc,
                           bld.env.mpi_mem_c_obj, 'fftw_mod_obj', 'NAG', 'base64',
                           'fxtp_wrap_obj', 'fxtp_obj', 'fxtp_wrapper',
                           'PRECICE','MPICXX', 'PYLIB', 'STDCXX', 'RT', 'ZLIB', 'PETSC', 
                           'BOOST_system', 'BOOST_filesystem'],
                target = 'atl_harvesting')
    else:

       bld.recurse('treelm', 'post_doxy')

       bld(
           features = 'coco',
           source = atl_ppsources)

def compile_atl(bld, ateles_sources):
    from waflib import Logs
    from waflib.extras.utest_results import utests

    bld(
        features = 'coco fc',
        source = ateles_sources,
        use = ['FFTW3', 'NAG'],
        target = 'atl_objs')

    test_dep = ['aotus', 'atl_objs', 'tem_objs', 'ply_objs', bld.env.distcrc,
                bld.env.mpi_mem_c_obj, 'fftw_mod_obj', 'NAG',
                'fxtp_wrap_obj', 'fxtp_obj', 'fxtp_wrapper',
                'PRECICE', 'MPICXX', 'PYLIB', 'STDCXX', 'RT', 'ZLIB']
    utests(bld = bld, use = test_dep, coco=True)

    if bld.env.LIB_FFTW3:
       utests(bld = bld, use = test_dep, path = 'utests/with_fftw')


# clean output files 
# add different extension format to remove
# by extending outputfiles list using outfiles.extend(glob.glob('*.yourfileext')
def cleanoutput(ctx):
    outputfiles=[]
    outputfiles.extend(glob.glob('*.vtk'))
    for output in outputfiles:
        os.remove(output)

#clean build directory and coco completely to create the build from scratch
def cleanall(ctx):
    from waflib import Options
    Options.commands = ['distclean'] + Options.commands
    ctx.exec_command('rm coco')
    
