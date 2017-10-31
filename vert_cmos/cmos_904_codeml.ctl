      seqfile = sate_904_forcodeml.fasta
     treefile = RAxML_bestTree.constrained_904
      outfile = cmos-904

        noisy = 9   * how much rubbish on the screen
      verbose = 1   * detailed output
      runmode = 0   * user tree 

      seqtype = 1   * codons
    CodonFreq = 2   * F3X4
        clock = 0   * no clock
        model = 0   * no variation in dN/dS among branches

      NSsites = 3   * discrete among-site dN/dS distribution w/ # of categories
                    * specified by ncatG (see below)
        icode = 0   * standard genetic code

    fix_kappa = 0   * kappa to be estimated
        kappa = 2   * initial kappa
    fix_omega = 0   * dN/dS to be estimated
        omega = .4  * initial dN/dS

    fix_alpha = 1   * no gamma-distributed rate heterogeneity
        alpha = 0   * no gamma-distributed rate heterogeneity
       Malpha = 0   * different alphas for genes
        ncatG = 3   * default # of dN/dS categories (as in Yang et al. 2000)

        getSE = 0   * no standard errors of estimates
 RateAncestor = 1   * force the estimation of ancestral states

   Small_Diff = .5e-6
    cleandata = 0   * do not remove ambiguous sites
  fix_blength = 1   * use the RAxML branch lengths as initial values
       method = 1   * update branch lengths one by one: faster when model = 0