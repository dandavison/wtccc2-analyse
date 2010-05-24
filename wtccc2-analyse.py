#!/usr/bin/env python

import sys, os, time
from cmdline.cmdline import CommandLineApp
import subprocess
from subprocess import PIPE
from dedpy.ded import *

__width__ = 30
__home__ = subprocess.Popen('echo $HOME', shell=True, stdout=PIPE).communicate()[0].strip()
__datadir__ = '/data/oak/project/wtccc2'
__dry_run__ = False
__insectdir__ = 'insect_out'
__progname__ = 'wtccc2-analyse'

class App(CommandLineApp):
    def __init__(self):
        CommandLineApp.__init__(self)
        
        op = self.option_parser
        usage = []
        usage.append("%s --snptest --controls 58C --cases NBS --chrom '21 22' --local"
                     % __progname__)
        usage.append('%s --cohorts "58C NBS POBI" --snpfile ~/wtccc2/celine-snps-ill'% __progname__)
        op.set_usage('\n' + '\n'.join(usage))

        op.add_option('--pca', dest='pca', default=False, action='store_true',
                      help='Perform PCA of selected cohorts')
        op.add_option('', '--make-geno', dest='make_geno', default=False, action='store_true',
                      help="Convert to .geno input data format used by Dan's PCA software " + \
                          "(one line per SNP, *no spaces*, {0,1,2,9})")
        op.add_option('--snptest', dest='snptest', default=False, action='store_true',
                      help='Perform association tests of selected cases & controls')
        op.add_option('--sstat', dest='sstat', default=False, action='store_true',
                      help='Compute summary statistics for selected cohorts')
        op.add_option('--cohorts', dest='cohorts', type='string', default='',
                      help='Space-separated list of cohorts')
        op.add_option('--cases', dest='cases', type='string', default='',
                      help='Space-separated list of case cohorts')
        op.add_option('--controls', dest='controls', type='string', default='',
                      help='Space-separated list of control cohorts') 
        op.add_option('--outfile', dest='outfile', type='string', default='results',
                      help='Name of output file')
        op.add_option('--local', dest='local', default=False, action='store_true',
                      help='Run PCA / snptest locally, not on the cluster')
        op.add_option('--chrom', dest='chroms', type='string', default=None,
                      help='Chromosomes to include (space-separated; default is 1-22)')
        op.add_option('--snpfile', dest='snpfile', type='string', default=None,
                      help='File containing SNPs to use (.map format)')
        op.add_option('--samplefile', dest='samplefile', type='string', default=None,
                      help='sample file to use (chiamo .sample format; defaults to WTCCC2 cohort sample file)')
        op.add_option('--snptest-opts', dest='snptest_opts', type='string', default='',
                      help='Options to pass to snptest')
        op.add_option('--factorfile', dest='factor_file', type='string', default=None,
                      help='File containing integer factor classifying individuals (minimum factor level is 1)')
        op.add_option('--exclude', dest='excludefile', type='string', default=None,
                      help='File containing IDs of individuals to exclude' + \
                          'in addition to project exclusions (one ID per line)')
        op.add_option('-n', dest='dry_run', default=False, action='store_true',
                      help="Don't actually do anything, just print commands")
        op.add_option('--platform', dest='platform', type='string', default=None)

    def main(self):
        if self.options.chroms:
            self.chroms = map(int, self.options.chroms.split())
        else:
            self.chroms = [c+1 for c in range(22)]
        self.snpfile = self.options.snpfile
        self.samplefile = self.options.samplefile
        self.platform = self.options.platform
        self.options.snptest_opts = self.options.snptest_opts.replace('*', ' ')
        self.options.combine_cohorts = self.options.pca or self.options.make_geno

        if self.options.pca or self.options.make_geno:
            self.analysis = 'PCA'
            self.cohorts = self.options.cohorts.split()
            self.say_hello()
            self.sanity_check()
            self.create_data_set()
            if self.options.pca:
                self.pca()
        elif self.options.snptest:
            self.analysis = 'snptest'
            self.cases = self.options.cases.split()
            self.controls = self.options.controls.split()
            self.cohorts = self.cases + self.controls
            self.say_hello()
            self.sanity_check()
            self.create_data_set()
            self.snptest()
        elif self.options.sstat:
            self.analysis = 'sstat'
            self.cohorts = self.options.cohorts.split()
            self.sanity_check()
            self.sstat()

    def sanity_check(self):
        actions = ['pca','snptest','sstat', 'make_geno']
        requested_actions = [getattr(self.options, action) for action in actions]
        if len(filter(None, requested_actions)) != 1:
            raise Exception('Use either %s' % ' or '.join(actions))
        if self.options.pca:
            if not self.options.cohorts:
                raise Exception('Select cohorts using --cohorts')
        elif self.options.snptest:
            if not self.options.cases:
                raise Exception('Select case cohorts using --cases')
            if not self.options.controls:
                raise Exception('Select control cohorts using --controls')
        elif self.options.sstat:
            if not self.options.cohorts:
                raise Exception('Select cohorts using --cohorts')
            if len(self.cohorts) != 1:
                print(self.cohorts)
                raise Exception('Select a single cohort with --sstat')
        if self.platform not in ['illumina','affymetrix']:
            raise Exception('Select platform using --platform illumina or --platform affymetrix')
        if not self.options.snpfile:
            raise Exception('Specify SNPs to include with --snpfile')

    def say_hello(self):
        print(time.ctime())
        print('Analysis'.ljust(__width__) + '%s' % self.analysis)
        print('Cohorts'.ljust(__width__) + '%s' % self.cohorts)
        print('Chromosomes'.ljust(__width__) + '%s' % self.chroms)
        print('SNP file'.ljust(__width__) + '%s' % self.snpfile)
        if self.options.dry_run:
            print('Dry run')

    def create_data_set(self):
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('Intersecting chromosome files\n')
        fnames = ['%s/%s-%02d.tmp' % (__insectdir__, coh, chrom) \
                      for coh in self.cohorts \
                      for chrom in self.chroms]
        if not all(map(os.path.exists, fnames)):
            self.insect_chromosome_files()

        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('Concatenating chromosomes\n')
        fnames = [coh + '.gen' for coh in self.cohorts]
        if not all(map(os.path.exists, fnames)):
            self.concatenate_chromosomes()
            system('rm %s/*' % __insectdir__)
            system('rmdir %s' % __insectdir__)

        def files_exist(bnames):
            format = 'geno' if self.options.pca else 'gen'
            geno = [b + '.' + format for b in bnames]
            sample = [b + '.sample' for b in bnames]
            maps = [b + '.map' for b in bnames] if self.options.pca else []
            return all(map(os.path.exists, flatten([geno, sample, maps])))

        rfiles = [restricted_genofile(coh, self.snpfile) for coh in self.cohorts]
        xfiles = [excluded_genofile(coh, self.snpfile) for coh in self.cohorts]

        if not (files_exist(rfiles) or files_exist(xfiles)):
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print('Restricting to selected SNPs\n')
            self.restrict_to_selected_SNPs()
        rmapfiles = [rfile + '.map' for rfile in rfiles]
        assert_files_identical(rmapfiles)

        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('Excluding individuals\n')
        if not files_exist(xfiles):
            self.exclude_individuals()

        if self.options.combine_cohorts:
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print('Combining data across cohorts\n')
            if not files_exist([excluded_genofile('all', self.snpfile)]):
                self.combine_cohorts()

    def insect_chromosome_files(self):
        outdir = __insectdir__
        if not os.path.exists(outdir): os.mkdir(outdir)
        for chrom in self.chroms:
            fnames = ['%s-%02d.tmp' % (coh, chrom) for coh in self.cohorts]
            for i in range(len(self.cohorts)):
                coh = self.cohorts[i]
                with open(fnames[i], 'w') as f:
                    Popen(['gunzip', '-vc', gen_gz_file(coh, chrom, self.platform)], stdout=f).communicate()
                    
            cmd = ['insect', '-v', '--unique', "-d ' '", '-f 2', '-o ' + outdir] + fnames
            # subprocess.Popen(cmd, shell=True).communicate()
            system(' '.join(cmd))
            map(os.remove, fnames)

    def concatenate_chromosomes(self):
        for coh in self.cohorts:
            with open(coh + '.gen', 'w') as f:
                cmd = 'cat %s/%s-*' % (__insectdir__, coh)
                Popen([cmd], shell=True, stdout=f).communicate()
            if not(os.path.exists(coh + '.sample')):
                os.symlink(sample_file(coh, self.platform), coh + '.sample')

    def restrict_to_selected_SNPs(self):
        for coh in self.cohorts:
            cmd = 'shellfish --make-%s --file %s %s --out %s' % \
                ('gen' if self.options.snptest else 'geno',
                 coh,
                 '--file2 %s' % self.snpfile if self.snpfile else '',
                 restricted_genofile(coh, self.snpfile) )
            print(cmd)
            system(cmd)
            system('mv %s.sample %s.sample' % (coh, restricted_genofile(coh, self.snpfile)))
            system('rm %s.gen %s.map' % (coh,coh))

    def exclude_individuals(self):
        for coh in self.cohorts:
            
            # Make sorted list of IDs to be excluded
            if self.samplefile is None or \
                    not os.path.exists(user_sample_file(self.samplefile, coh)):
                project_excludeglob =  exclude_dir(coh, self.platform) + '/*.exclude.txt'
                ## TODO: test for non-empty glob expansion
                cmd = 'cat %s %s | sort | uniq > %s.xids' % \
                    (project_excludeglob, self.options.excludefile or "", coh)
            else:
                ## Get IDs to keep
                tempfile = "/tmp/%s.ids" % coh
                cmd = "sed 1,2d %s | cut -d ' ' -f 1 > %s ; " % \
                    (user_sample_file(self.samplefile, coh), tempfile)
                cmd += "sed 1,2d %s | cut -d ' ' -f 1 | grep -vf %s > %s.xids" % \
                    (sample_file(coh, self.platform), tempfile, coh)   
            system(cmd, verbose=True)

            # Get cohort indices of individuals to be excluded
            # These are the (line index in sample file) - 2, because sample file has 2 header lines.
            cmd = "sed 1,2d %s | cut -d ' ' -f 1 | match %s.xids > %s.xidx" % \
                (sample_file(coh, self.platform), coh, coh)
            system(cmd, verbose=True)

            # Check for IDs that did not appear in cohort sample file
            cmd = 'echo "%s: `grep -F NA %s.xidx  | wc -l` excluded individuals not recognised"' % \
                (coh, coh)
            system(cmd)
            cmd = 'grep -vF NA %s.xidx | sort -n > tmp && mv tmp %s.xidx' % \
                (coh, coh)
            system(cmd, verbose=True)

            if self.options.snptest:
                format = 'gen'
                # Compute columns of .gen file to be excluded
                idx = map(int, read_lines('%s.xidx' % coh))
                firstofthree = [6 + (i-1)*3 for i in idx]
                idx = flatten([range(s, s+3) for s in firstofthree])
                write_lines(idx, '%s.xidx' % coh)
            else:
                format = 'geno'

            # Exclude individuals from genotype data
            cmd = 'columns %s -v -f %s.xidx < %s.%s > %s.%s' % (
                '-s' if format == 'gen' else '',
                coh,
                restricted_genofile(coh, self.snpfile), format,
                excluded_genofile(coh, self.snpfile), format)
            system(cmd, verbose=True)
                
            # Get IDs of included individuals
            cmd = "sed 1,2d %s | cut -d ' ' -f 1 | slice -v --line-file %s.xidx > %s.ids" % \
                (sample_file(coh, self.platform), coh, excluded_genofile(coh, self.snpfile))
            system(cmd, verbose=True)

            system('rm %s.%s' % (restricted_genofile(coh, self.snpfile), format), verbose=True)
            if self.options.pca:
                system('mv %s.map %s.map' % (
                        restricted_genofile(coh, self.snpfile),
                        excluded_genofile(coh, self.snpfile)), verbose=True)
            system('mv %s.sample %s.sample' % (
                    restricted_genofile(coh, self.snpfile),
                    excluded_genofile(coh, self.snpfile)), verbose=True)

    def combine_cohorts(self):
        geno_files = [excluded_genofile(coh, self.snpfile) + '.geno' for coh in self.cohorts]
        map_files = [excluded_genofile(coh, self.snpfile) + '.map' for coh in self.cohorts]
        id_files = [excluded_genofile(coh, self.snpfile) + '.ids' for coh in self.cohorts]
        combined_basename = excluded_genofile('all', self.snpfile)

        ## Combine genotype data across cohorts
        cmd = "paste -d '\\0' %s > %s.geno" % (' '.join(geno_files), combined_basename)
        system(cmd, verbose=True)

        ## Combine .ids files
        cmd = 'cat %s > %s.ids' % (' '.join(id_files), combined_basename)
        system(cmd, verbose=True)
        
        ## Propagate .map file
        cmd = 'cp %s %s.map' % (map_files[0], combined_basename)
        system(cmd, verbose=True)

        ## Clean up
        map(os.remove, geno_files)
        map(os.remove, map_files)
        map(os.remove, id_files)

    def snptest(self):
        case_files = [excluded_genofile(coh, self.snpfile) for coh in self.cases]
        control_files = [excluded_genofile(coh, self.snpfile) for coh in self.controls]

        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

        if self.options.local:
            print('Running shellfish on local machine\n')
            cmd = 'shellfish --snptest --maxprocs 1'
            cmd += ' --cases ' + ' '.join(case_files) + ' --controls ' + ' '.join(control_files)
            cmd += ' --outfile %s' % self.options.outfile
            cmd += ' --snptest-chunk 1000 --snptest-opts %s' % self.options.snptest_opts
            system(cmd, verbose=True)

        else:
            print('Running shellfish on remote machine\n')
    
            remote = 'login2-cluster1'
            remote_dir = 'shellfish-%s' % datetimenow()
    
            case_files_string = ' '.join([f + '.gen' for f in case_files])
            control_files_string = ' '.join([f + '.gen' for f in control_files])
    
            cmd = "ssh %s 'mkdir -p %s'" % (remote, remote_dir)
            system(cmd, verbose=True)
            
            cmd = 'scp %s %s %s:%s/' % (case_files_string, control_files_string, remote, remote_dir)
            system(cmd, verbose=True)
            
            remote_cmd = 'shellfish --snptest --sge --sge-level 2 --maxprocs 100'
            remote_cmd += ' --cases ' + ' '.join(case_files) + ' --controls ' + ' '.join(control_files)
            remote_cmd += ' --outfile %s/%s' % (remote_dir, self.options.outfile)
            remote_cmd += ' --snptest-chunk 1000 --snptest-opts %s' % self.options.snptest_opts
            remote_cmd = "'nohup %s < /dev/null > %s/log 2>&1 &'" % (remote_cmd, remote_dir)
            
            cmd = 'ssh %s %s' % (remote, remote_cmd)
            system(cmd, verbose=True)

    def pca(self):
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        remote = 'login2-cluster1'
        remote_dir = 'shellfish-%s' % datetimenow()
        print('Running shellfish on %s in %s\n' % (remote, remote_dir))

        if not os.path.exists(excluded_genofile('all', self.snpfile) + '.evecs'):

            cmd = "ssh %s 'mkdir -p %s'" % (remote, remote_dir)
            system(cmd)
            
            tup = ((excluded_genofile('all', self.snpfile),) * 3) + (remote, remote_dir)
            cmd = 'scp %s.geno %s.map %s.ids %s:%s/' % tup
            system(cmd)
            
            remote_cmd = "shellfish --pca --sge --sge-level 2 --numpcs 10 --maxprocs 500 "
            remote_cmd += "--file %s --out %s" % ((excluded_genofile('all', self.snpfile),) * 2)
            remote_cmd = "'cd %s && nohup %s < /dev/null > log 2>&1'" % (remote_dir, remote_cmd)

            cmd = 'ssh %s %s &' % (remote, remote_cmd)
            system(cmd)

    def sstat(self):
        """Run sstat on each cohort file for each chromosome"""
        coh = self.cohorts[0]
        nsample = count_lines(sample_file(coh, self.platform)) - 2 
        nfac = count_lines(self.options.factor_file)
        if nsample != nfac:
            raise Exception('Number of individuals in sample file (%d) does not match number if factor file (%d)' % (
                    (nsample, nfac)))
        for chrom in self.chroms:
            system('gunzip -c %s | sstat -n %d -p -f %s > %s-%02d.sstat' % (
                    gen_gz_file(coh, chrom, self.platform), nsample, self.options.factor_file, coh, chrom), verbose=True)

def gen_gz_file(coh, chrom, platform):
    return '%s/%s/%s/calls/%s_%02d_%s.gen.gz' % \
        (__datadir__, coh, platform, coh, chrom, platform)

def sample_file(coh, platform):
    return '%s/%s/%s/calls/%s_%s.sample' % \
        (__datadir__, coh, platform, coh, platform)

def user_sample_file(basename, coh):
    return '%s.%s' % (basename, coh)

def restricted_genofile(coh, snpfile):
    f = coh
    if snpfile:
        f += '-' + os.path.basename(snpfile)
    return f

def exclude_dir(coh, platform):
    return '%s/%s/%s/exclusions' % (__datadir__, coh, platform)

def excluded_genofile(coh, snpfile):
    f = coh + 'x'
    if snpfile:
        f += '-' + os.path.basename(snpfile)
    return f

def Popen(cmd, shell=False, stdout=None):
    print(' '.join(cmd) + (' > ' + stdout.name if stdout else ''))
    if app.options.dry_run:
        return subprocess.Popen('', shell=True)
    else:
        return subprocess.Popen(cmd, shell=shell, stdout=stdout)

if __name__ == '__main__':
      app = App()
      # app.options, main_args = app.option_parser.parse_args()      
      app.run()
