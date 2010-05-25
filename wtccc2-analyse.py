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
        op.add_option('', '--make-gen', dest='make_gen', default=False, action='store_true',
                      help="Convert to .gen input data format used by chiamo and related software")
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
        op.add_option('--outfile', dest='outfile', type='string', default=None,
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
        global opts
        opts = self.options
        if opts.chroms:
            opts.chroms = map(int, opts.chroms.split())
        else:
            opts.chroms = [c+1 for c in range(22)]
        opts.snptest_opts = opts.snptest_opts.replace('*', ' ')
        opts.combine_cohorts = opts.pca or opts.make_gen
        self.format = 'geno' if opts.pca else 'gen'
        self.insect_dir = opts.outfile + '-' + 'insect_out'

        if opts.pca or opts.make_gen:
            self.analysis = 'PCA'
            self.cohorts = opts.cohorts.split()
            self.say_hello()
            self.sanity_check()
            self.create_data_set()
            if opts.pca:
                self.pca()
        elif opts.snptest:
            self.analysis = 'snptest'
            self.cases = opts.cases.split()
            self.controls = opts.controls.split()
            self.cohorts = self.cases + self.controls
            self.say_hello()
            self.sanity_check()
            self.create_data_set()
            self.snptest()
        elif opts.sstat:
            self.analysis = 'sstat'
            self.cohorts = opts.cohorts.split()
            self.sanity_check()
            self.sstat()

    def sanity_check(self):
        actions = ['pca','snptest','sstat', 'make_gen']
        requested_actions = [getattr(opts, action) for action in actions]
        if len(filter(None, requested_actions)) != 1:
            raise Exception('Use either %s' % ' or '.join(actions))
        if opts.pca:
            if not opts.cohorts:
                raise Exception('Select cohorts using --cohorts')
        elif opts.snptest:
            if not opts.cases:
                raise Exception('Select case cohorts using --cases')
            if not opts.controls:
                raise Exception('Select control cohorts using --controls')
        elif opts.sstat:
            if not opts.cohorts:
                raise Exception('Select cohorts using --cohorts')
            if len(self.cohorts) != 1:
                print(self.cohorts)
                raise Exception('Select a single cohort with --sstat')
        if opts.platform not in ['illumina','affymetrix']:
            raise Exception('Select platform using --platform illumina or --platform affymetrix')
        if opts.outfile is None:
            raise Exception('Supply output filename prefix with --outfile')

    def say_hello(self):
        print(time.ctime())
        print('Analysis'.ljust(__width__) + '%s' % self.analysis)
        print('Cohorts'.ljust(__width__) + '%s' % self.cohorts)
        print('Chromosomes'.ljust(__width__) + '%s' % opts.chroms)
        print('SNP file'.ljust(__width__) + '%s' % opts.snpfile)
        print('Output file/prefix'.ljust(__width__) + '%s' % opts.outfile)
        if opts.dry_run:
            print('Dry run')

    def create_data_set(self):
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('Intersecting chromosome files\n')
        fnames = ['%s/%s-%02d.tmp' % (self.insect_dir, coh, chrom) \
                      for coh in self.cohorts \
                      for chrom in opts.chroms]
        if not all(map(os.path.exists, fnames)):
            self.insect_chromosome_files()

        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('Concatenating chromosomes\n')
        fnames = [coh + '.gen' for coh in self.cohorts]
        if not all(map(os.path.exists, fnames)):
            self.concatenate_chromosomes()
            system('rm %s/*' % self.insect_dir)
            system('rmdir %s' % self.insect_dir)

        def files_exist(bnames):
            geno = [b + '.' + self.format for b in bnames]
            sample = [b + '.sample' for b in bnames]
            maps = [b + '.map' for b in bnames] if self.format == 'geno' else []
            return all(map(os.path.exists, flatten([geno, sample, maps])))

        rfiles = [self.restricted_genofile(coh) for coh in self.cohorts]
        xfiles = [self.excluded_genofile(coh) for coh in self.cohorts]

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

        if opts.combine_cohorts:
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print('Combining data across cohorts\n')
            if not files_exist([self.excluded_genofile('all')]):
                self.combine_cohorts()

    def insect_chromosome_files(self):
        outdir = self.insect_dir
        if not os.path.exists(outdir): os.mkdir(outdir)
        for chrom in opts.chroms:
            fnames = ['%s-%s-%02d.tmp' % (opts.outfile, coh, chrom) for coh in self.cohorts]
            for i in range(len(self.cohorts)):
                coh = self.cohorts[i]
                with open(fnames[i], 'w') as f:
                    Popen(['gunzip', '-vc', gen_gz_file(coh, chrom, opts.platform)], stdout=f).communicate()
                    
            cmd = ['insect', '-v', '--unique', "-d ' '", '-f 2', '-o ' + outdir] + fnames
            # subprocess.Popen(cmd, shell=True).communicate()
            system(' '.join(cmd))
            map(os.remove, fnames)

    def concatenate_chromosomes(self):
        for coh in self.cohorts:
            gen_file = opts.outfile + '-' + coh + '.gen'
            sample_file = opts.outfile + '-' + coh + '.sample'
            with open(gen_file, 'w') as f:
                cmd = 'cat %s/%s-%s-*' % (self.insect_dir, opts.outfile, coh)
                Popen([cmd], shell=True, stdout=f).communicate()
            if not(os.path.exists(sample_file)):
                os.symlink(wtccc2_sample_file(coh, opts.platform), sample_file)

    def restrict_to_selected_SNPs(self):
        for coh in self.cohorts:
            cmd = 'shellfish --make-%s --file %s %s --out %s' % \
                (self.format,
                 coh,
                 '--file2 %s' % opts.snpfile if opts.snpfile else '',
                 self.restricted_genofile(coh) )
            system(cmd, verbose=True)
            system('mv %s.sample %s.sample' % (coh, self.restricted_genofile(coh)))
            if opts.snpfile or self.format == 'geno': # In which case a new genotype file has been created
                system('rm %s.gen' % coh)

    def exclude_individuals(self):
        for coh in self.cohorts:
            
            # Make sorted list of IDs to be excluded
            if opts.samplefile is None or \
                    not os.path.exists(user_sample_file(opts.samplefile, coh)):
                project_excludeglob =  exclude_dir(coh, opts.platform) + '/*.exclude.txt'
                ## TODO: test for non-empty glob expansion
                cmd = 'cat %s %s | sort | uniq > %s.xids' % \
                    (project_excludeglob, opts.excludefile or "", coh)
            else:
                ## Get IDs to keep
                tempfile = "/tmp/%s-%s.ids" % (opts.outfile, coh)
                cmd = "sed 1,2d %s | cut -d ' ' -f 1 > %s ; " % \
                    (user_sample_file(opts.samplefile, coh), tempfile)
                cmd += "sed 1,2d %s | cut -d ' ' -f 1 | grep -vf %s > %s.xids" % \
                    (wtccc2_sample_file(coh, opts.platform), tempfile, coh)   
            system(cmd, verbose=True)

            # Get cohort indices of individuals to be excluded
            # These are the (line index in sample file) - 2, because sample file has 2 header lines.
            cmd = "sed 1,2d %s | cut -d ' ' -f 1 | match %s.xids > %s.xidx" % \
                (wtccc2_sample_file(coh, opts.platform), coh, coh)
            system(cmd, verbose=True)

            # Check for IDs that did not appear in cohort sample file
            cmd = 'echo "%s: `grep -F NA %s.xidx  | wc -l` excluded individuals not recognised"' % \
                (coh, coh)
            system(cmd)
            cmd = 'grep -vF NA %s.xidx | sort -n > %s-tmp && mv %s-tmp %s.xidx' % \
                (coh, coh, opts.outfile, opts.outfile)
            system(cmd, verbose=True)

            if self.format == 'gen':
                # Compute columns of .gen file to be excluded
                idx = map(int, read_lines('%s.xidx' % coh))
                firstofthree = [6 + (i-1)*3 for i in idx]
                idx = flatten([range(s, s+3) for s in firstofthree])
                write_lines(map(str, idx), '%s.xidx' % coh)

            # Exclude individuals from genotype data
            cmd = 'columns %s -v -f %s.xidx < %s.%s > %s.%s' % (
                '-s' if self.format == 'gen' else '',
                coh,
                self.restricted_genofile(coh), self.format,
                self.excluded_genofile(coh), self.format)
            system(cmd, verbose=True)
                
            # Get IDs of included individuals
            cmd = "sed 1,2d %s | cut -d ' ' -f 1 | slice -v --line-file %s.xidx > %s.ids" % \
                (wtccc2_sample_file(coh, opts.platform), coh, self.excluded_genofile(coh))
            system(cmd, verbose=True)

            # clean up
            system('rm %s.%s' % (self.restricted_genofile(coh), self.format), verbose=True)
            system('rm %s.xids %s.xidx' % (coh, coh))

            if self.format == 'geno':
                system('mv %s.map %s.map' % (
                        self.restricted_genofile(coh),
                        self.excluded_genofile(coh)), verbose=True)
            system('mv %s.sample %s.sample' % (
                    self.restricted_genofile(coh),
                    self.excluded_genofile(coh)), verbose=True)

    def combine_cohorts(self):
        gen_files = [self.excluded_genofile(coh) + '.' + self.format for coh in self.cohorts]
        map_files = [self.excluded_genofile(coh) + '.map' for coh in self.cohorts]
        id_files = [self.excluded_genofile(coh) + '.ids' for coh in self.cohorts]
        sample_files = [self.excluded_genofile(coh) + '.sample' for coh in self.cohorts]
        combined_basename = self.excluded_genofile('all')

        ## Combine genotype data across cohorts
        if self.format == 'geno':
            cmd = "paste -d '\\0' %s > %s.geno" % (' '.join(gen_files), combined_basename)
            system(cmd, verbose=True)
    
            ## Propagate .map file
            cmd = 'cp %s %s.map' % (map_files[0], combined_basename)
            system(cmd, verbose=True)

        else:
            gen_only_files = [self.excluded_genofile(coh) + '.gen_only' for coh in self.cohorts]
            for (gen_file, map_file, gen_only_file) in zip(gen_files, map_files, gen_only_files):
                system("cut -d ' ' -f 1-5 < %s > %s" % (gen_file, map_file), verbose=True)
                system("cut -d ' ' -f 6- < %s > %s" % (gen_file, gen_only_file), verbose=True)
            assert_files_identical(map_files)
            cmd = "paste -d ' ' %s %s > %s.gen" % (
                gen_files[0], ' '.join(gen_only_files[1:]), combined_basename)
            system(cmd, verbose=True)
            map(os.remove, gen_only_files)

        ## Combine .ids files
        cmd = 'cat %s > %s.ids' % (' '.join(id_files), combined_basename)
        system(cmd, verbose=True)
            
        ## Clean up
        map(os.remove, gen_files)
        map(os.remove, map_files)
        map(os.remove, id_files)
        map(os.remove, sample_files)

    def snptest(self):
        case_files = [self.excluded_genofile(coh) for coh in self.cases]
        control_files = [self.excluded_genofile(coh) for coh in self.controls]

        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

        if opts.local:
            print('Running shellfish on local machine\n')
            cmd = 'shellfish --snptest --maxprocs 1'
            cmd += ' --cases ' + ' '.join(case_files) + ' --controls ' + ' '.join(control_files)
            cmd += ' --outfile %s' % opts.outfile
            cmd += ' --snptest-chunk 1000 --snptest-opts %s' % opts.snptest_opts
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
            remote_cmd += ' --outfile %s/%s' % (remote_dir, opts.outfile)
            remote_cmd += ' --snptest-chunk 1000 --snptest-opts %s' % opts.snptest_opts
            remote_cmd = "'nohup %s < /dev/null > %s/log 2>&1 &'" % (remote_cmd, remote_dir)
            
            cmd = 'ssh %s %s' % (remote, remote_cmd)
            system(cmd, verbose=True)

    def pca(self):
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        remote = 'login2-cluster1'
        remote_dir = 'shellfish-%s' % datetimenow()
        print('Running shellfish on %s in %s\n' % (remote, remote_dir))

        if not os.path.exists(self.excluded_genofile('all') + '.evecs'):

            cmd = "ssh %s 'mkdir -p %s'" % (remote, remote_dir)
            system(cmd)
            
            tup = ((self.excluded_genofile('all'),) * 3) + (remote, remote_dir)
            cmd = 'scp %s.geno %s.map %s.ids %s:%s/' % tup
            system(cmd)
            
            remote_cmd = "shellfish --pca --sge --sge-level 2 --numpcs 10 --maxprocs 500 "
            remote_cmd += "--file %s --out %s" % ((self.excluded_genofile('all'),) * 2)
            remote_cmd = "'cd %s && nohup %s < /dev/null > log 2>&1'" % (remote_dir, remote_cmd)

            cmd = 'ssh %s %s &' % (remote, remote_cmd)
            system(cmd)

    def sstat(self):
        """Run sstat on each cohort file for each chromosome"""
        coh = self.cohorts[0]
        nsample = count_lines(wtccc2_sample_file(coh, opts.platform)) - 2 
        nfac = count_lines(opts.factor_file)
        if nsample != nfac:
            raise Exception('Number of individuals in sample file (%d) does not match number if factor file (%d)' % (
                    (nsample, nfac)))
        for chrom in opts.chroms:
            system('gunzip -c %s | sstat -n %d -p -f %s > %s-%02d.sstat' % (
                    gen_gz_file(coh, chrom, opts.platform), nsample, opts.factor_file, coh, chrom), verbose=True)

    def restricted_genofile(self, coh):
        f = opts.outfile + '-' + coh + 'r'
        if opts.snpfile:
            f += '-' + os.path.basename(opts.snpfile)
        return f
    
    def excluded_genofile(self, coh):
        f = opts.outfile + '-' + coh + 'x'
        if opts.snpfile:
            f += '-' + os.path.basename(opts.snpfile)
        return f

def gen_gz_file(coh, chrom, platform):
    return '%s/%s/%s/calls/%s_%02d_%s.gen.gz' % \
        (__datadir__, coh, platform, coh, chrom, platform)

def wtccc2_sample_file(coh, platform):
    return '%s/%s/%s/calls/%s_%s.sample' % \
        (__datadir__, coh, platform, coh, platform)

def user_sample_file(basename, coh):
    return '%s.%s' % (basename, coh)

def exclude_dir(coh, platform):
    return '%s/%s/%s/exclusions' % (__datadir__, coh, platform)

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
