import os
import sys
import uuid
import numpy as np

from glob import glob
from collections import defaultdict

PROG_NAME   = 'ms_fraction_merge'
DESCRIPTION = 'This software merges fractionated mass spectrometry abundance profiles from replicate experiments where fraction data is not immediately comparable between replicates and must first be aligned to combine equivalent fractions'

DEFAULT_MIN_REPS  = 3
DEFAULT_OUT_FRACS = 25
DEFAULT_OUT_FILE  = 'merged_replicates.csv'

QUIET   = False
LOGGING = False
LOG_FILE_PATH = 'ms-frac-merge-out-%s.log' % uuid.uuid4()
LOG_FILE_OBJ = None # Created when needed

DEFAULT_OFFSET_MAX  = 5.0
DEFAULT_OFFSET_STEP = 0.5
DEFAULT_SCALE_MAX   = 2.5
DEFAULT_SCALE_STEP  = 0.2


def report(msg):
 
  if LOGGING:
    if not LOG_FILE_OBJ:
      LOG_FILE_OBJ = open(LOG_FILE_PATH, 'w')
      
    LOG_FILE_OBJ.write(msg)
  
  if not QUIET:
    print(msg)


def warn(msg, prefix='WARNING'):

  report('%s: %s' % (prefix, msg))

 
def critical(msg, prefix='FAILURE'):

  report('%s: %s' % (prefix, msg))
  sys.exit(0)


def info(msg, prefix='INFO'):

  report('%s: %s' % (prefix, msg))


def check_regular_file(file_path):

  msg = ''
  
  if not os.path.exists(file_path):
    msg = 'File "%s" does not exist'
    critical(msg % file_path)
  
  if not os.path.isfile(file_path):
    msg = 'Location "%s" is not a regular file'
    critical(msg % file_path)
  
  if os.stat(file_path).st_size == 0:
    msg = 'File "%s" is of zero size '
    critical(msg % file_path)
    
  if not os.access(file_path, os.R_OK):
    msg = 'File "%s" is not readable'
    critical(msg % file_path)


def load_spectral_counts(file_paths, min_reps, clip_fracs):
  """
  Load replicate spectral count data from CSV files into NumPy arrays.
  Data will be first row normalised (per protein) by adjusting the minimum to zero and dividing by the
  row sum. Columns are scaled by dividing by their median.
  """
  
  data = []

  col_max = 0
  col_min = sys.maxint
  rep_counts = defaultdict(int)
  
  for r, rep_file in enumerate(file_paths):
    rep_data = {}
    data.append(rep_data)
    clip_frac = clip_fracs[r]
    
    with open(rep_file, 'rU') as file_obj:
      head = file_obj.readline().strip()
      
      if head.count(',') > head.count('\t'):
        sep = ','
      else:
        sep = '\t'  
      
      head = head.split(sep)
      n_cols = len(head)
      
      cols = [int(x) for x in head[1:]]
      col_min = min(col_min, min(cols))
      col_max = max(col_max, max(cols))
 
      for i, line in enumerate(file_obj):
        vals = line.strip().split(sep)
        
        if len(vals) != n_cols:
          msg = 'Input data format problem at %s line %d (and maybe others thereafter) ' % rep_file
          msg += 'Expected %d fields (where first is protein ID) but found %d' % (n_cols, len(vals))
          critical(msg)
        
        for i, val in enumerate(vals):
          if val[0] == '"':
            val = val[1:-1] 
          
          val = val.replace('""', '"')
          val = val.strip()
          vals[i] = val
        
        pid = vals[0]
        vals = [float(x) for x in vals[1:]]

        prot_dict = {}
        for i, col in enumerate(cols):
          
          if (clip_frac is not None) and (col >= clip_frac):
            break
          
          prot_dict[col] = vals[i]
         
        rep_data[pid] = prot_dict
        rep_counts[pid] += 1
  
  # Get common proteins with representation in the stated number of replicates
  common_pids = {pid for pid in rep_counts if rep_counts[pid] >= min_reps}
  common_pids = sorted(common_pids)
  
  cols = np.arange(col_min, col_max+1)
  n_cols = len(cols)
  n_prot = len(common_pids)

  info("  loaded data for %d components over a total of %d fraction columns" % (n_prot, n_cols))
  
  comp_matrices = []
  comp_matrices_orig = []
  
  for r, rep_data in enumerate(data):
    matrix = np.zeros((n_cols, n_prot), float)
 
    for j, pid in enumerate(common_pids):
      if pid in rep_data:
        for col in  rep_data[pid]:
          i = col - col_min
          matrix[i, j] = rep_data[pid].get(col, 0.0)
               
    comp_matrices_orig.append(matrix.copy())
    
    for i in range(n_prot):
      matrix[:,i] -= matrix[:,i].min()
      matrix[:,i] /= matrix[:,i].sum() or 1.0

    for i in range(n_cols):
      col = matrix[i]
      col = col[col.nonzero()]
      
      if len(col):
        f = np.median(col)
        matrix[i] /= float(f) or 1.0
      
    comp_matrices.append(matrix)

  return comp_matrices, comp_matrices_orig, common_pids, rep_counts


def get_comp_matrix_regions(comp_matrix):
  """
  Small function to get the starting fraction region boundaries for
  a given composition matrix
  """
  
  n_col, n_prot = comp_matrix.shape
  regions = np.array(range(0, n_col+1), float)

  return regions


def interpolate_comp_matrix_regions(comp_matrix):
  """
  Function to add interpolated fractional composition data to missing, intermediate fractions
  """
  
  n_col, n_prot = comp_matrix.shape
  col_sums = comp_matrix.sum(axis=1)
  comp_matrix_2 = comp_matrix.copy()
  
  msg = "   interpolating composition matrix over %d columns and %d components (%d empty columns)"
  
  is_zero = col_sums == 0.0
  is_zero_idx = is_zero.nonzero()[0]
  
  info(msg % (n_col, n_prot, len(is_zero_idx)))
  
  for i in range(n_col):
    if i > 0:
      if col_sums[i] == 0:
        j = i
        k = i
        while j >= 0 and col_sums[j] == 0:
          j -= 1
        
        if col_sums[j] == 0:
          continue

        while k < n_col-1 and col_sums[k] == 0:
          k += 1
        
        if col_sums[k] == 0:
          continue
        
        f = float(i-j)
        g = float(k-i)
        t = f+g
        
        f /= t
        g /= t
        
        comp_matrix_2[i] = g * comp_matrix[j] + f * comp_matrix[k]
        
  col_sums_2 = comp_matrix_2.sum(axis=1)
 
  i = 0
  j = n_col - 1
  
  while (col_sums_2[i] == 0) and (i < j):
    i += 1
  
  while (col_sums_2[j] == 0) and (j > i):
    j -= 1
  
  return comp_matrix_2[i:j+1]
  
  
def get_sim_score(abun_a, abun_b, weight, default=-1.0):
  """
  Function to get the width-weighted similarity score (Pearson's correlation) between abundance profiles
  """
  
  idx = (abun_a * abun_b).nonzero()[0]  
  
  if len(idx):
    score = np.corrcoef(abun_a[idx], abun_b[idx])[0,1]
    score *= weight # Scale score according to overlap region width
  
  else:
    score = 0.0
 
  return score


def overlap_region_comps(regions_a, comps_a, regions_b, comps_b,
                         contrib_a=1.0, contrib_b=1.0, insignificant=1e-4):
  """
  Function to calculate a score for a trial overlap between two fraction region
  specifications with corresponding component/protein compositions
  """
  
  x_min = min(regions_a[ 0], regions_b[ 0])
  x_max = max(regions_a[-1], regions_b[-1])
  
  p = 0
  
  n = len(regions_a)
  m = len(regions_b)
  
  zeros_a = set([i for i in range(n-1) if comps_a[i].sum() == 0])
  zeros_b = set([j for j in range(m-1) if comps_b[j].sum() == 0])
  
  # Some proteins are missing in some replicates
  n_prot = len(comps_a[0]) 
  missing_a = [ i for i in range(n_prot) if comps_a[:,i].sum() == 0] # Missing across all fractions
  missing_b = [ j for j in range(n_prot) if comps_b[:,j].sum() == 0]
    
  merge_regions = [x_min, ]
  
  sim_score = 0.0
  score_width = 0.0  
  sim_scores = []
  
  # First regions
  
  if regions_b[0] < regions_a[0]:
    i = 0
    j = 1
    merge_comps = [comps_b[0],]
    sim_scores.append(0.0)
    
  elif regions_a[0] < regions_b[0]:
    i = 1
    j = 0
    merge_comps = [comps_a[0],]
    sim_scores.append(0.0)
  
  else:
    i = 1
    j = 1    
    
    comp = (contrib_a * comps_a[0] + contrib_b * comps_b[0])/(contrib_a+contrib_b)

    if missing_a: # Missing proteins are not averaged with zero
      comp[missing_a] = comps_b[0][missing_a] 

    if missing_b:
      comp[missing_b] = comps_a[0][missing_b]
    
    merge_comps = [comp]
    delta = min(regions_a[1], regions_b[1]) - x_min
    s = get_sim_score(comps_a[0], comps_b[0], delta)
    sim_score += s
    sim_scores.append(s/delta)
    score_width += delta
    p += 1
  
  x = x_min
  
  # Intermediate regions
  
  while (i < n-1) or (j < m-1):
    #info(i, j, sim_score, sim_score/(score_width or 1.0))
    
    if j >= m-1:
      merge_regions.append(regions_a[i])
      merge_comps.append(comps_a[i])
      sim_scores.append(0.0)
      i += 1
   
    elif i >= n-1:
      merge_regions.append(regions_b[j])
      merge_comps.append(comps_b[j])
      sim_scores.append(0.0)
      j += 1
         
    elif regions_a[i] < regions_b[j]:
      merge_regions.append(regions_a[i])
      if j == 0:
        sim_scores.append(0.0)
        merge_comps.append(comps_a[i])
      else:      
        if i in zeros_a:
          sim_scores.append(0.0)
          merge_comps.append(comps_b[j])
 
        elif j in zeros_b:
          sim_scores.append(0.0)
          merge_comps.append(comps_a[i])
 
        else:
          delta = merge_regions[-1] - merge_regions[-2]
          
          if delta > insignificant: # Ignore, insignificant small merged fraction widths 
            
            # Contribution weighted compositions
            comp = (contrib_a*comps_a[i] + contrib_b*comps_b[j])/(contrib_a+contrib_b)
            
            # Some proteins are missing in some replicates
            if missing_a: # Missing proteins are not averaged with zero
              comp[missing_a] = comps_b[j][missing_a]

            if missing_b:
              comp[missing_b] = comps_a[i][missing_b]
            
            merge_comps.append(comp)
           
            s = get_sim_score(comps_a[i], comps_b[j], delta)
            sim_score += s
            sim_scores.append(s/delta)
            score_width += delta
            p += 1
            
          else:
            merge_regions.pop()
 
        
      i += 1
      
    elif regions_b[j] < regions_a[i]:
      merge_regions.append(regions_b[j])
       
      if i == 0:
        merge_comps.append(comps_b[j])
        sim_scores.append(0.0)
      
      else:
        if i in zeros_a:
          merge_comps.append(comps_b[j])
          sim_scores.append(0.0)

        elif j in zeros_b:
          merge_comps.append(comps_a[i])
          sim_scores.append(0.0)
 
        else:
          delta = merge_regions[-1] - merge_regions[-2]          
          
          if delta > 1e-4:
            weight_a = contrib_a #* delta/(regions_a[i]-regions_a[i-1])
            weight_b = contrib_b #* delta/(regions_b[j]-regions_b[j-1])
            comp = (weight_a*comps_a[i] + weight_b*comps_b[j])/(weight_a+weight_b)
            
            if missing_a: # Missing proteins are not averaged with zero
              comp[missing_a] = comps_b[j][missing_a]

            if missing_b:
              comp[missing_b] = comps_a[i][missing_b]
             
            merge_comps.append(comp)

            s = get_sim_score(comps_a[i], comps_b[j], delta)
            sim_score += s
            sim_scores.append(s/delta)
            score_width += delta
            p += 1
          
          else:
            merge_regions.pop()
          
      j += 1
    
    else:
      merge_regions.append(regions_a[i])
      
      if i in zeros_a:
        merge_comps.append(comps_b[j])
        sim_scores.append(0.0)
        
      elif j in zeros_b:
        merge_comps.append(comps_a[i])
        sim_scores.append(0.0)
      
      else:
        delta = merge_regions[-1] - merge_regions[-2]
        
        if delta > 1e-4:
          comp = (contrib_a * comps_a[i] + contrib_b * comps_b[j])/(contrib_a+contrib_b)
            
          if missing_a: # Missing proteins are not averaged with zero
            comp[missing_a] = comps_b[j][missing_a]

          if missing_b:
            comp[missing_b] = comps_a[i][missing_b]
         
          merge_comps.append(comp)
 
          s = get_sim_score(comps_a[i], comps_b[j], delta)
            
          sim_score += s
          sim_scores.append(s/delta)
          score_width += delta
          p += 1
        else:
          merge_regions.pop()
        
      i += 1
      j += 1
 
  merge_regions.append(x_max)
  
  # Last region
  
  if (regions_a[-1] == x_max) and (regions_b[-1] == x_max):
    merge_comps.append((contrib_a * comps_a[-1] + contrib_b * comps_b[-1])/(contrib_a+contrib_b))
    delta = merge_regions[-1] - merge_regions[-2]
    
    s = get_sim_score(comps_a[-1], comps_b[-1], delta)
    sim_score += s
    sim_scores.append(s/delta)
    score_width += delta
    p += 1
    
  elif regions_a[-1] == x_max:
    merge_comps.append(comps_a[-1])
    sim_scores.append(0.0)
  
  else:
    merge_comps.append(comps_b[-1])
    sim_scores.append(0.0)

  sim_score /= score_width or 1.0
    
  sim_scores = np.array(sim_scores)
  merge_regions = np.array(merge_regions)
  merge_comps =  np.array(merge_comps)
  
  return sim_scores, merge_regions, merge_comps, sim_score, score_width, p
  

def merge_replicates(regions_a, comps_a, regions_b, comps_b, contrib_a=1.0, contrib_b=1.0,
                     offset_max=4.0, offset_step=0.5, scale_max=2.0, scale_step=0.1):
   """
   Function to do the exhaustive offset and scale parameter search to find the best
   alignment of input fractions.
   """
   
   best_score = -1.0
   best_offset = 0.0
   best_scale_start = -1.0
   best_scale_end = -1.0
   best_comps = None
   best_scores = None
   best_regions = None

   n_b = len(regions_b)
   points_b = np.arange(n_b)
   
   min_scale = 1.0/scale_max
   
   for offset in np.arange(-offset_max, offset_max, offset_step):

     for scale_start in np.arange(min_scale, scale_max, scale_step):
       
       for scale_end in np.arange(min_scale, scale_max, scale_step):
                   
         stretch = np.interp(points_b, [0, n_b], [scale_start, scale_end])
 
         stretch_regions = (regions_b * stretch) + offset   
 
         d = overlap_region_comps(regions_a, comps_a,
                                  stretch_regions, comps_b,
                                  contrib_a, contrib_b)
                                  
         sim_scores, merge_regions, merge_comps, sim_score, score_width, p = d
         
         if sim_score > best_score:
           best_scores = sim_scores
           best_score = sim_score
           best_scale_start = scale_start
           best_scale_end = scale_end
           best_offset = offset
           best_comps = merge_comps
           best_regions = merge_regions

   info('    best score: %5.2f at offset: %5.2f start scale: %5.2f end scale: %5.2f' % (best_score, best_offset, best_scale_start, best_scale_end))
   
   return best_regions, best_comps, best_score


def sample_region_comps(regions, comps, n_samples):
  """
  Function to generate a merged pseudo-fraction data set from aligned
  composition data by sampling regularly spaced averages.
  """
  
  n_cols, n_prot = comps.shape
  x_min = regions[0]
  x_max = regions[-1]
  dx = (x_max - x_min)/n_samples
  sampled_comp = np.zeros((n_samples, n_prot), float)
  weights = np.zeros(n_samples, float)
  
  n = len(regions)-1
  
  for i in range(n):
    x1 = regions[i]
    x2 = regions[i+1]
    
    p1 = int((x1-x_min)/dx)
    p2 = int((x2-x_min)/dx)
    
    if p1 == p2:
      sampled_comp[p1] += comps[i]
      weights[p1] += 1.0
    
    else:
      p = p1 + 1
      x = p * dx
      sampled_comp[p1] += comps[i]
      weights[p1] += 1.0
      x0 = x
      
      while (p <= p2) and (p < n_samples):
        x = min(x2, x0+dx)
        sampled_comp[p] += comps[i]
        weights[p] += 1.0
        p += 1   
        
  for i in range(n_samples):
    sampled_comp[i] /= weights[i]
      
  return sampled_comp   

  
def ms_fraction_merge(rep_csv_paths, out_csv_path, marker_csv_paths,
                      num_pseudo_fracs=DEFAULT_OUT_FRACS, min_num_reps=DEFAULT_MIN_REPS, clip_fracs=None,
                      offset_max=DEFAULT_OFFSET_MAX, offset_step=DEFAULT_OFFSET_STEP,
                      scale_max=DEFAULT_SCALE_MAX, scale_step=DEFAULT_SCALE_STEP):
  
  if len(rep_csv_paths) < 2:
    critical('Must have at least two input replicate CSV files')
  
  for file_path in rep_csv_paths:
    check_regular_file(file_path)
  
  if offset_max < 0.0:
    critical('Maximum offset must not be negative')

  if offset_step < 0.01:
    critical('Offset step must be at least 0.01')

  if offset_step >= offset_max:
    critical('Offset step must smaller than the maximum offset')
  
  if scale_max < 1.0:
    critical('Maximum scale must be at least 1.0')
  
  if scale_step < 0.01:
    critical('Scale step must be at least 0.01')

  if scale_step >= scale_max:
    critical('Scale step must be smaller than the maximum scale')
  
  num_reps = len(rep_csv_paths)
  
  if min_num_reps > num_reps:
    min_num_reps = num_reps
    warn('Minimum number of replicates set to %d' % num_reps)
  
  if not clip_fracs:
    clip_fracs = []
  
  while len(clip_fracs) < num_reps:
    clip_fracs.append(None)

  info('Loading and normalising spectra count data')

  comp_matrices, comp_matrices_orig, pids, rep_counts = load_spectral_counts(rep_csv_paths, min_num_reps, clip_fracs)
  
  
  info('Interpolating missing intermediate column/fraction data')  
  
  comp_matrices = [interpolate_comp_matrix_regions(mat) for mat in comp_matrices]
  
  
  info('Merging replicate data')  
  
  regions = [get_comp_matrix_regions(mat) for mat in comp_matrices]

  categ_dict = {}

  for marker_path in marker_csv_paths:
    with open(marker_path) as file_obj:
      
      
      for line in file_obj:
        if '\t' in line:
          data = line.strip().split()
        else:
          data = line.strip().split(',')
          
        pid = data[0]
        pid = pid.split('.')[0]
        
        compartment = data[1]
        categ_dict[pid] = compartment
        
        # Tweak for protein isoform IDs
        for j in range(1,5):
          categ_dict[pid + '.%d' % j] = compartment

  info('Course comparison of replicate pairs')  
  
  pair_scores = []
  for i in range(num_reps-1):
    for j in range(i+1, num_reps):
       # Initial scoring alignment is course for speed
        
       info('  comparing %d:%d' % (i+1, j+1))
       regs, comps, score = merge_replicates(regions[i], comp_matrices[i], regions[j], comp_matrices[j], 1.0, 1.0,
                                             offset_max, 1.0, scale_max, 0.5)
       pair_scores.append((score, regs, comps, i, j,))
  
  
  info('Optimising replicate alignment')  
  
  merged = set()
  merge_regions = None
  merge_comps = None
  pair_scores.sort()
  
  while pair_scores: # Merge in score order, adding only to the previous merged data
    score, regs, comps, i, j = pair_scores.pop() # (Next) highest score
    
    if merged: # Merge to previous
      if (i in merged) and (j in merged):
        continue # Already added in some other order
       
      if i in merged:
        k = j
      
      elif j in merged:
        k = i
      
      else: 
        continue # Not a merge to previous
      
      weight = float(len(merged))
      merged.add(k)
      
      info('  merging %s' % '+'.join([str(x+1) for x in sorted(merged)]))
      
      merge_regions, merge_comps, score = merge_replicates(merge_regions, merge_comps, regions[k], comp_matrices[k],
                                                           weight, 1.0, offset_max, offset_step, scale_max, scale_step)
       
    else: # First starting pair, accept merge
      merged = {i,j}
      info('  merging %s' % '+'.join([str(x+1) for x in sorted(merged)]))
      
      merge_regions, merge_comps, score = merge_replicates(regions[i], comp_matrices[i], regions[j], comp_matrices[j],
                                                           1.0, 1.0, offset_max, offset_step, scale_max, scale_step)
  
  info('Resampling merged data over %d pseudo-fractions' % num_pseudo_fracs)

  resampled_comp = sample_region_comps(merge_regions, merge_comps, num_pseudo_fracs)

  n_cols = len(resampled_comp)

  with open(out_csv_path, 'w') as file_obj:
    head = ['component', 'category', 'num_reps'] + [str(i) for i in range(n_cols)]
    line = ','.join(head) + '\n'
    file_obj.write(line)
 
    for i, pid in enumerate(pids):
      abun = resampled_comp[:,i]
      category = categ_dict.get(pid, 'unknown')
      n_rep = rep_counts[pid]

      data = [pid, category, '%d' % n_rep] + ['%.7f' % x for x in abun]
      line = ','.join(data) + '\n'
      file_obj.write(line)


  info('Merged data written to %s' % out_csv_path)
  

if __name__ == '__main__':

  from argparse import ArgumentParser
   
  epilog = 'For further help on running this program please email tjs23@cam.ac.uk.\n\n'
  epilog += 'Example use:\n\n'
  epilog += 'ms_fraction_merge example_data/test_spec_count_data_*.csv -o merged_replicates.csv -c example_data/test_categories.tsv -m 2 -clip 50 47 37'

  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('rep_csvs', nargs='+', metavar='CSV_FILES',
                         help='Input file paths of CSV files containing spectral count data, one for each replicate (may contain wildcards)') 

  arg_parse.add_argument('-o', metavar='OUT_CSV_FILE', default=DEFAULT_OUT_FILE,
                         help='Optional putput file path for CSV file containing merged spectral count data; Default %s' % DEFAULT_OUT_FILE) 

  arg_parse.add_argument('-c', nargs='+', metavar='CATEGORY_CSV_FILE', default=None,
                         help='One or more optional CSV files describing categories/classes, e.g. for marker proteins, relating to the input data (may contain wildcards)') 

  arg_parse.add_argument('-n', metavar='NUM_FRACTIONS', default=DEFAULT_OUT_FRACS, type=int,
                         help='Number of output (pseudo-)fraction columns to output for the merged data. Default: %d' % DEFAULT_OUT_FRACS) 

  arg_parse.add_argument('-m', metavar='MIN_REPS', default=DEFAULT_MIN_REPS, type=int,
                         help='Minimum number of replicate occurrences required for a protein to be accepted. Default: %d' % DEFAULT_MIN_REPS) 

  arg_parse.add_argument('-omax', metavar='MAX_OFFSET', default=DEFAULT_OFFSET_MAX, type=float,
                         help='The maximum number of whole columns that could offset the first values between different replicate data; defines the width of the offset parameter search space. Default: %.2f' % DEFAULT_OFFSET_MAX) 
  
  arg_parse.add_argument('-ostep', metavar='STEP', default=DEFAULT_OFFSET_STEP, type=float,
                         help='The fractional increment in column widths to use when aligning replicate data; defines the granularty of the offset parameter search space. Default: %.2f' % DEFAULT_OFFSET_STEP) 
  
  arg_parse.add_argument('-smax', metavar='MAX_SCALE', default=DEFAULT_SCALE_MAX, type=float,
                         help='The maximum stretch/compression scale factor between different replicate data; defines the width of the scale parameter search space. Default: %.2f' % DEFAULT_SCALE_MAX) 
  
  arg_parse.add_argument('-sstep', metavar='STEP', default=DEFAULT_SCALE_STEP, type=float,
                         help='The increment in stretch scale factor use when aligning replicate data; defines the granularty of the scale parameter search space. Default: %.2f' % DEFAULT_SCALE_STEP) 
   
  arg_parse.add_argument('-clip', nargs='+', metavar='FRACTION_NUMS', default=None, type=int,
                         help='Optional fraction numbers to clip input data at (one for each replicate in input order), i.e. beyond which there is no useful data') 
  
  arg_parse.add_argument('-q', default=False, action='store_true',
                         help='Sets quiet mode to supress on-screen reporting.')
  
  arg_parse.add_argument('-log', default=False, action='store_true',
                         help='Log all reported output to a file.')
  
  args = vars(arg_parse.parse_args())

  rep_csvs   = args['rep_csvs']
  out_csv    = args['o']
  mark_csvs  = args['c'] or []
  out_fracs  = args['n']
  min_reps   = args['m']
  clip_fracs = args['clip']
  off_max    = args['omax']
  off_step   = args['ostep']
  scale_max  = args['smax']
  scale_step = args['sstep']
  
  # Globals
  QUIET   = args['q']
  LOGGING = args['log']  
  
  ms_fraction_merge(rep_csvs, out_csv, mark_csvs, out_fracs, min_reps, clip_fracs, off_max, off_step, scale_max, scale_step)

