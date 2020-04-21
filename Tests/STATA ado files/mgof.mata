*! version 1.0.0  Ben Jann  12jan2008

* all functions are part of the -moremata- package (using different names)
* see -ssc describe moremata-

version 10.0
mata:

struct SJmgof_subsetinfo {
    real scalar    n, k, i, j, counter, algorithm
    real colvector x, y
}
real matrix SJmgof(
    real colvector f,    // observed count
  | real colvector h0,   // expected distribution
    string scalar m0,    // method: "approx" (default), "mc", or "ee"
    string vector s0,    // test-stats: "x2", "lr", "cr", "mlnp", "ksmirnov"
    real scalar lambda,  // lambda for Cressie-Read statistic
    real scalar arg1,    // nfit for "approx"; reps for "mc"
    real scalar dots,    // display dots with mc/ee
    real scalar arg2)    // will be replaced; used by mgof_ee to return counter
{
    real scalar     i, n, hsum
    string scalar   m, s
    string vector   stats
    real vector     p
    pointer(real colvector) scalar       h
    pointer(real scalar function) vector s1

// parsing and defaults
    if (args()<7) dots = 0
    m = SJmgof_strexpand(m0, ("approx", "mc", "ee"), "approx", 0)
    stats = ("x2", "lr", "cr")
    if (m!="approx") stats = (stats, "mlnp", "ksmirnov")
    if (length(s0)<1) s = "x2"
    else {
        s = J(length(s0),1,"")
        for (i=1;i<=length(s0);i++) {
            s[i] = SJmgof_strexpand(s0[i], stats, "x2", 1)
        }
    }
    if (m!="approx" & lambda<0) _error(3498, "lambda<0 not allowed with "+m)

// check data / prepare null distribution
    if (missing(f)) _error(3351, "f has missing values")
    if (any(f:<0)) _error(3498, "f has negative values")
    if (any(f:<0)) _error(3498, "f has negative values")
    if (m=="ee" | (m=="mc" & anyof(s, "mlnp"))) {
        if (f!=trunc(f)) _error(3498, "non-integer f not allowed in this context")
    }
    if (anyof(stats,"cr") & lambda<0 & anyof(f,0))
     _error(3498, "some f are 0: lambda<0 not allowed")
    if (args()<2) h0 = 1
    if (rows(f)!=rows(h0) & rows(h0)!=1) _error(3200)
    if (missing(h0)) _error(3351, "h has missing values")
    if (any(h0:<=0)) _error(3498, "some h are negative or 0")
    n = colsum(f)
    hsum = colsum(h0)
    if (n!=hsum) h = &(h0*n/hsum)
    else h = &h0

// return (0,1) if n. of obs is 0 or n. of categories < 2
    if ( n==0 | rows(f)<2 ) return(J(rows(s),1,0),J(rows(s),1,1))

// gather statistics functions
    s1 = J(rows(s),1,NULL)
    for (i=1;i<=length(s);i++) {
        if      (s[i]=="x2")        s1[i] = &SJmgof_x2()
        else if (s[i]=="lr")        s1[i] = &SJmgof_lr()
        else if (s[i]=="cr")        s1[i] = &SJmgof_cr()
        else if (s[i]=="mlnp")      s1[i] = &SJmgof_mlnp()
        else if (s[i]=="ksmirnov")  s1[i] = &SJmgof_ksmirnov()
        else _error(3498,s[i]+ " invalid")
    }

// approximate chi2-test
    if (m=="approx") return(SJmgof_approx(s1, f, *h, lambda, arg1))

// Monte Carlo test
    if (m=="mc") return(SJmgof_mc(s1, f, *h, lambda, arg1, dots))

// exhaustive enumeration / random composition test
    if (s=="mlnp") {
        return(SJmgof_ee(NULL, f, *h))
    }
    p = SJmgof_which(s:=="mlnp")
    if (rows(p)>0) {
        p = p[1]
        p = p \ (p>1 ? (1::p-1) : J(0,1,.)) \ (p<rows(s) ? (p+1::rows(s)) : J(0,1,.))
        s1 = s1[p[| 2 \ .|]]
        p = invorder(p)
    }
    else p = (2::length(s)+1)
    return(SJmgof_ee(s1, f, *h, lambda, anyof(s,"ksmirnov"), dots, arg2)[p,])
}

// large sample chi2-approximation test for multinomial distributions
real matrix SJmgof_approx(
 pointer(real scalar function) vector s, // statistics to be tested
 real colvector f,     // observed count
 real colvector h0,    // expected count
 | real scalar lambda, // lambda for Cressie-Read statistic
   real scalar nfit)   // n. of fitted parameters (default 0)
{
    real scalar                     n, i, k, nstat
    real colvector                  stat
    pointer(real colvector) scalar  h

// compute hypothetical distribution if missing
    n = colsum(f)
    k = rows(f)
    if (rows(h0)<=1) h = &(n :* J(k,1,1/k))
    else             h = &h0
// compute statistics and p-values
    nstat = length(s)
    stat = J(nstat,1,.)
    for (i=1; i<=nstat; i++) stat[i] = (*s[i])(f,*h,lambda)
    return(stat, chi2tail(rows(f) - (nfit<. ? nfit : 0) - 1,stat))
}

// Monte Carlo exact test for discrete distributions
// (based on sampling from the hypothetical distribution)
real matrix SJmgof_mc(
 pointer(real scalar function) vector s, // statistics to be tested
 real colvector f,     // observed count
 real colvector h0,    // expected count
 | real scalar lambda, // lambda for Cressie-Read statistic
   real scalar reps0,  // replications (default=10000)
   real scalar dots)   // display dots
{
    real scalar                     n, j, i, k, nstat, reps, uni
    real scalar                     doti, ndots
    real colvector                  fj, H, stat, p, hp
    pointer(real colvector) scalar  h

// set defaults and compute hypothetical distribution
    if (args()<5 | reps0>=.) reps = 10000
    else                     reps = reps0
    if (args()<6) dots = 0
    n = round(colsum(f)) // so that, e.g., 10 are sampled if n=9.9999...
    k = rows(f)
    if (uni = (rows(h0)==1)) {
        uni = 1
        h  = &(n :* J(k,1,1/k))
        H  = rangen(1/k,1,k)
        hp = J(k,1,1/k)
    }
    else {
        uni = SJmgof_isconstant(h0)
        H  = runningsum(h0)
        H  = H / H[rows(H)]
        hp = (h0) / colsum(h0)
        h  = &h0
    }
// compute statistics and p-values
    nstat = length(s)
    stat = p = J(nstat,1,0)
    for (i=1; i<=nstat; i++) stat[i] = (*s[i])(f,*h, lambda, n, hp, H)
    if (reps<=0) return(stat, J(nstat,1,.))
    if (dots) {
        printf("\n{txt}Percent completed ({res}%g{txt} replications)\n",reps)
        display("{txt}0 {hline 5} 20 {hline 6} 40 {hline 6} 60 {hline 6} 80 {hline 5} 100")
        ndots = (reps>=50 ? 1 : (reps>=25 ? 2 : (reps>=10 ? 5 :
                (reps>=5 ? 10 : (reps>=2 ? 25 : 50)))))
        doti = 1
    }
    for (j=1; j<=reps; j++) {
        if (uni) fj = SJmgof_mc_usmpl(n, k) // = mm_srswr(n,rows(f),1)
        else     fj = SJmgof_mc_smpl(n, H) // = mm_upswr(n, h, 1)  (slower)
        for (i=1; i<=nstat; i++)  p[i] = p[i] +
          (stat[i] <= (*s[i])(fj,*h, lambda, n, hp, H))
        if (dots) {
            if (j*50 >= doti*ndots*reps) {
                printf(ndots*".")
                displayflush()
                doti++
            }
        }
    }
    if (dots) display("")
    return(stat, p/reps)
}
// Monte Carlo subroutine: sample from uniform multinomial
real colvector SJmgof_mc_usmpl(
 real scalar n, // sample size
 real scalar k) // number of categories
{
    real scalar     i
    real colvector  u, res

    u = ceil(uniform(n,1)*k)
    res = J(k,1,0)
    if (n>1.5^(k-2)) { // faster for small k relative to n
        for (i=1; i<=k; i++) res[i] = colsum(u:==i)
    }
    else {
        for (i=1;i<=n;i++) res[u[i]] = res[u[i]] + 1
    }
    return(res)
}
// Monte Carlo subroutine: sample from non-uniform multinomial
real colvector SJmgof_mc_smpl(
 real scalar n,    // sample size
 real colvector F) // distribution function (cumulative)
{
    real scalar     i, j
    real colvector  u, res

    if (n>2.3^(rows(F)-3)) {  // faster for small k relative to n
        u = uniform(n,1)
        res = J(rows(F),1,0)
        res[1] = sum(u:<=F[1])
        for (j=2;j<=rows(F);j++) {
            res[j] = sum(u:>F[j-1] :& u:<=F[j])
        }
        return(res)
    }
    u = sort(uniform(n,1),1)
    j=1
    res = J(rows(F),1,0)
    for (i=1;i<=n;i++) {
        while (u[i]>F[j]) j++
        res[j] = res[j]+1
    }
    return(res)
}

// exhaustive enumeration exact test (performs the exact
// multinomial g.o.f. test plus, optionally, exact tests
// for specified additional statistics)
// WARNING: use this test for very small samples only since
// computation time is linear to the number of compositions
// =comb(n+k-1,k-1), where n is the sample size and k is
// the number of categories.
// IMPORTANT EXCEPTION: when the null distribution is uniform,
// compution time is determined by the number of partitions
// which is magnitudes smaller than the number of compositions
real matrix SJmgof_ee(
 pointer(real scalar function) vector s0, // statistics to be tested
 real colvector f,     // observed count
 real colvector h0,    // expected count
 | real scalar lambda, // lambda for Cressie-Read statistic
   real scalar force,  // do not use partitions (for ksmirnov)
   real scalar dots,   // display dots
   real scalar counter) // will be replaced
{
    real scalar                           n, i, k, nstat, uni, w
    real scalar                           doti, ndots, reps
    real colvector                        fj, j, H, stat, statj, p, hp
    pointer(real colvector) scalar        h
    pointer(real scalar function) vector  s
    struct SJmgof_subsetinfo scalar           info

// set defaults and compute hypothetical distribution
    if (args()<6) dots = 0
    if (s0==NULL) s = &SJmgof_mlnp()
    else {
        if (rows(s0)!=1) s = &SJmgof_mlnp() \ s0
        else             s = &SJmgof_mlnp() , s0
    }
    n = colsum(f)
    k = rows(f)
    uni = (force==1 ? 0 : 1)
    if (rows(h0)==1) {
        h = &(n :* J(k,1,1/k))
        H = rangen(1/k,1,k)
        hp = J(k,1,1/k)
    }
    else {
        h = &h0
        if (SJmgof_isconstant(*h)==0) uni = 0 // reset to 0 if h not constant
        H = runningsum(*h)
        H = H / H[rows(H)]
        hp = (*h) / colsum(*h)
    }
// compute statistics and p-values
    nstat = length(s)
    stat = statj = p = j = J(nstat,1,0)
    for (i=1; i<=nstat; i++)  stat[i] = (*s[i])(f,*h, lambda, n, hp, H)
    if (uni) {         // uniform null distribution: work with partitions
        if (dots) {
            reps = SJmgof_npartitions(n,k)
            printf("\n{txt}Percent completed ({res}%g{txt} partitions)\n",reps)
            display("{txt}0 {hline 5} 20 {hline 6} 40 {hline 6} 60 {hline 6} 80 {hline 5} 100")
            ndots = (reps>=50 ? 1 : (reps>=25 ? 2 : (reps>=10 ? 5 :
                    (reps>=5 ? 10 : (reps>=2 ? 25 : 50)))))
            doti = 1
        }
        info = SJmgof_partitionsetup(n, k)
        while ((fj = SJmgof_partition(info,1)) != J(0,1,.)) {
            w = exp(lnfactorial(info.k) - quadcolsum(lnfactorial(SJmgof_panels(fj))))
            for (i=1; i<=nstat; i++) {
                statj[i] = (*s[i])(fj,*h, lambda, n, hp, H)
                if (stat[i] <= statj[i]) {
                    j[i] = j[i] + w // statj[1] = -ln(p) => p = exp(-statj[1])
                    p[i] = p[i] + (exp(-statj[1])-p[i])*w/j[i]  // mean updating
                }
            }
            if (dots) {
                if (info.counter*50 >= doti*ndots*reps) {
                    printf(ndots*".")
                    displayflush()
                    doti++
                }
            }
        }
        if (dots) display("")
        counter = info.counter
        return(stat, p:*j)
    }
    if (dots) {
        reps = SJmgof_ncompositions(n,k)
        printf("\n{txt}Percent completed ({res}%g{txt} compositions)\n",reps)
        display("{txt}0 {hline 5} 20 {hline 6} 40 {hline 6} 60 {hline 6} 80 {hline 5} 100")
        ndots = (reps>=50 ? 1 : (reps>=25 ? 2 : (reps>=10 ? 5 :
                (reps>=5 ? 10 : (reps>=2 ? 25 : 50)))))
        doti = 1
    }
    info = SJmgof_compositionsetup(n, k)
    while ((fj = SJmgof_composition(info)) != J(0,1,.)) {
        for (i=1; i<=nstat; i++) {
            statj[i] = (*s[i])(fj,*h, lambda, n, hp, H)
            if (stat[i] <= statj[i]) {
                j[i] = j[i] + 1  // statj[1] = -ln(p) => p = exp(-statj[1])
                p[i] = p[i] + (exp(-statj[1])-p[i])/j[i]  // mean updating
            }
        }
        if (dots) {
            if (info.counter*50 >= doti*ndots*reps) {
                printf(ndots*".")
                displayflush()
                doti++
            }
        }
    }
    if (dots) display("")
    counter = info.counter
    return(stat, p:*j)
}

// Pearson's chi2 statistic subroutine
real scalar SJmgof_x2(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | transmorphic matrix opt1, // not used
      transmorphic matrix opt2, // not used
      transmorphic matrix opt3, // not used
      transmorphic matrix opt4) // not used
{
    return(quadcolsum((f-h):^2:/h))
}

// log likelihood ratio statistic subroutine
real scalar SJmgof_lr(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | transmorphic matrix opt1, // not used
      transmorphic matrix opt2, // not used
      transmorphic matrix opt3, // not used
      transmorphic matrix opt4) // not used
{
    return(2*quadcolsum(f:*(ln(f)-ln(h))))
}

// Cressie-Read statistic
real scalar SJmgof_cr(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | real scalar lambda,       // lambda parameter; default is 2/3
      transmorphic matrix opt2, // not used
      transmorphic matrix opt3, // not used
      transmorphic matrix opt4) // not used
{
    real scalar L
    L = (lambda<. ? lambda : 2/3)
    if (abs(L)<1e-6)        return( 2 * quadcolsum(f:*(ln(f)-ln(h))) )
    else if (abs(L+1)<1e-6) return( 2 * quadcolsum(h:*(ln(h)-ln(f))) )
    else if (L>0) return( 2/(L*(L+1)) * quadcolsum(f:*((f:/h):^L :- 1)) )
                  return( 2/(L*(L+1)) * quadcolsum(f:*((h:/f):^-L :- 1)) )
}

// -ln(p) statistic subroutine (multinomial test)
real scalar SJmgof_mlnp(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | transmorphic matrix opt1, // not used
      real scalar n,            // n. of obs [ =colsum(f) ]
      real colvector p,         // hypothical propotions [ = h/colsum(h) ]
      transmorphic matrix opt4) // not used
{
    if (args()>=5) return(-(lnfactorial(n) -
     quadcolsum(lnfactorial(f)) + quadcolsum(f:*ln(p))))
    return(-(lnfactorial(quadcolsum(f)) -
     quadcolsum(lnfactorial(f)) + quadcolsum(f:*ln(h/quadcolsum(h)))))
}

// kolmogorov-smirnov statistic subroutine
real scalar SJmgof_ksmirnov(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | transmorphic matrix opt1, // not used
      real scalar n,            // n. of obs [ =colsum(f) ]
      transmorphic matrix opt3, // not used
      real colvector H)         // hypothetical cumulative [ = runsum(h)/sum(h) ]
{
    if (args()>=6) return(max(abs(H-runningsum(f)/n)))
    return(max(abs(runningsum(h)/colsum(h)-runningsum(f)/colsum(f))))
}


// combinatorial algorithms
real scalar SJmgof_ncompositions(
 real scalar n,
 | real scalar k0)
{
    real scalar k

    k = trunc(k0<. ? k0 :  n)
    return(comb(trunc(n)+k-1,k-1))
}
struct SJmgof_subsetinfo scalar SJmgof_compositionsetup(
 real scalar n,
 | real scalar k)
{
    struct SJmgof_subsetinfo scalar s

    s.n = trunc(n)
    s.k = (k<. ? trunc(k) :  s.n)
    if (s.n>=.) _error(3351)
    if (s.k<1) {
        s.i = 0
        s.x = J(0,1,.)
    }
    else if (s.n<1) {
        s.i = 0
        s.x = J(s.k,1,0)
    }
    else if (s.k<2) {
        s.i = 0
        s.x = J(s.k,1,s.n)
    }
    else {
        s.x = s.n \ J(s.k-1,1,0)
        s.i = 1
    }
    s.j = 2
    s.counter = 0
    return(s)
}
real colvector SJmgof_composition(struct SJmgof_subsetinfo scalar s)
{
    real colvector subset

    subset = s.x
    if (s.i<1) {
        if (subset!=J(0,1,.)) s.counter = s.counter + 1
        s.x = J(0,1,.)
        return(subset) /*done*/
    }
    while (s.x[s.i]==0 & s.i>1) { // move back
        s.x[s.i] = s.x[s.j]
        s.x[s.j] = 0
        s.j = s.i--
    }
    if (s.x[s.i]>0) {
        s.x[s.i] = s.x[s.i]-1
        s.x[s.j] = s.x[s.j]+1
        if (s.j<s.k) s.i = s.j++
    }
    else {
        s.i = 0
        s.x = J(0,1,.)
    }
    s.counter = s.counter + 1
    return(subset)
}
real scalar SJmgof_npartitions(
 real scalar n,
 | real scalar k)
{
    real scalar      u, i, j, m
    real colvector   p, a

    m = min((n,k))
    if (m==n) return(_SJmgof_npartitions(n))
    p = J(n+1,1,1)
    for (u=2;u<=m;u++) {
        a = p
        p = J(n+1,1,0)
        for (i=0;i<=n;i=i+u) {
            for (j=i+1;j<=n+1;j++) {
                p[j] = p[j] + a[j-i]
            }
        }
    }
    return(p[n+1])
}
real scalar _SJmgof_npartitions(
 real scalar n)
{
    real scalar      i, j, s, k
    real colvector   p

    p = J(n+1,1,1)
    for (i=1;i<=n;i++) {
        j = 1 ; k = 1 ; s = 0
        while (j>0) {
            j = i - (3*k*k + k)/2
            if (j>=0) s = s - (-1)^k * p[j+1]
            j = i - (3*k*k - k)/2
            if (j>=0) s = s - (-1)^k * p[j+1]
            k = k + 1
        }
        p[i+1] = s
    }
    return(p[n+1])
}
struct SJmgof_subsetinfo scalar SJmgof_partitionsetup(
 real scalar n,
 | real scalar k)
{
    struct SJmgof_subsetinfo scalar s

    s.n = trunc(n)
    s.k = trunc(k)
    if (s.n>=.) _error(3351)
    if (s.k>=.) s.k = s.n
//  else if (s.k>s.n) _error(3300,"k may not be larger than n")
    if (s.n<1 | s.k<1 ) s.i = -1
    else {
        s.x = s.n \ J(s.k-1, 1, 1)
        s.i = 1
    }
    s.j = 1
    s.counter = 0
    return(s)
}
real colvector SJmgof_partition(     // original ZS1 Algorithm does not
 struct SJmgof_subsetinfo scalar s,  // support k-way partitions; changes
 real scalar pad)                    // are indicated
{
    real scalar    r, t, l
    real colvector subset

    if (s.i==-1) return(J(0,1,.)) /*done*/
    subset = s.x[|1 \ s.j|]
    if (pad) subset = subset \ J(s.k-rows(subset),1,0)
    l = (s.k<s.n ? s.k : s.n)     // added
    if (s.j>=l) {                 // added; original stopping rule is (s.x[1]==1)
        if (s.x[1]<=s.x[l]+1) {      // changed
            s.i = -1
            s.counter = s.counter + 1
            return(subset)
        }
        while (s.x[s.i]<=s.x[l]+1) {  // step back loop added
            s.j = s.j + s.x[s.i] - 1
            s.i = s.i - 1
        }
    }
    if (s.x[s.i]==2) {
        s.j      = s.j + 1
        s.x[s.i] = 1
        s.x[s.j] = 1             // added
        s.i      = s.i - 1
    }
    else {
        r = s.x[s.i] - 1
        t = s.j - s.i + 1
        s.x[s.i] = r
        while (t>=r) {
            s.i = s.i + 1
            s.x[s.i] = r
            t = t - r
        }
        if (t==0) {
            s.j = s.i
        }
        else {
            s.j = s.i + 1
            if (t==1) s.x[s.i+1] = 1    // added
            if (t>1) {
                s.i    = s.i + 1
                s.x[s.i] = t
            }
        }
    }
    s.counter = s.counter + 1
    return(subset)
}

// other functions
real scalar SJmgof_isconstant(X)
{
    if (length(X)<2) return(1)
    return(all(X:==X[1,1]))
}
real colvector SJmgof_freq(transmorphic matrix x,
 | real colvector w, transmorphic matrix levels)
{
    real colvector p

    if (args()<2) w = 1
    if (args()<3) levels = .
    if (cols(x)==0) return(_SJmgof_freq(x, w, levels))
    if (rows(w)==1) return(_SJmgof_freq(sort(x,1..cols(x)), w, levels))
    p = order(x,1..cols(x))
    return(_SJmgof_freq(x[p,], w[p,], levels))
}

real colvector _SJmgof_freq(transmorphic matrix x,
 | real colvector w, transmorphic matrix levels)
{
    real scalar    i, j, l
    real colvector result

    if (args()<2) w = 1
    if (args()<3) levels = .
    if (rows(w)!=1 & rows(w)!=rows(x)) _error(3200)
    if (levels==.) levels = _SJmgof_uniqrows(x)
    if (rows(x)==0) return(J(0,1, .))
    l = rows(levels)
    result = J(l,1,0)
    j = 1
    for (i=1; i<=rows(x); i++) {
        for (;j<=l;j++) {
            if (x[i,]==levels[j,]) break
        }
        if (j>l) break
        result[j] = result[j] + (rows(w)!=1 ? w[i] : w)
    }
    return(result)
}
transmorphic matrix _SJmgof_uniqrows( // uniqrows() for sorted X
 transmorphic matrix x)
{
        real scalar             i, j, n, ns
        transmorphic matrix     res

        if (rows(x)==0) return(J(0,cols(x), missingof(x)))
        if (cols(x)==0) return(J(1,0, missingof(x)))

        ns = 1
        n = rows(x)
        for (i=2;i<=n;i++) {
                if (x[i-1,]!=x[i,]) ns++
        }
        res = J(ns, cols(x), x[1,1])
        res[1,] = x[1,]
        for (i=j=2;i<=n;i++) {
                if (x[i-1,]!=x[i,]) res[j++,] = x[i,]
        }
        return(res)
}
string scalar SJmgof_strexpand(string scalar s, string vector slist,
 | string scalar def, real scalar unique, string scalar errtxt)
{
    real scalar   err
    string scalar res

    if (args()<5) errtxt = `"""' + s + `"" invalid"'
    if (args()<4) unique = 0
    err = _SJmgof_strexpand(res, s, slist, def, unique)
    if (err) _error(err, errtxt)
    return(res)
}

real scalar _SJmgof_strexpand(res, string scalar s,
 string vector slist, | string scalar def, real scalar unique)
{
    real scalar i, l, match

    if (s=="") {
        res = def
        return(0)
    }
    if (args()<5) unique = 0
    l = strlen(s)
    if (unique) {
        match = 0
        for (i=1; i<=length(slist); i++) {
            if (s==substr(slist[i], 1, l)) {
                if (match) return(3498)
                match = i
            }
        }
        if (match) {
            res = slist[match]
            return(0)
        }
    }
    else {
        for (i=1; i<=length(slist); i++) {
            if (s==substr(slist[i], 1, l)) {
                res = slist[i]
                return(0)
            }
        }
    }
    return(3499)
}
real matrix SJmgof_which(real vector I)
{
    if (cols(I)!=1) return(select(1..cols(I), I))
    else return(select(1::rows(I), I))
}
real colvector SJmgof_panels(transmorphic vector X, | real scalar np)
{
    real scalar i, j, r
    real colvector res

    r = length(X)
    if (r<1) return(J(0,1,.))
    if (args()<2) np = r
    res = J(np, 1, 1)
    j = 1
    for (i=2; i<=r; i++) {
        if (X[i]==X[i-1]) res[j] = res[j] + 1
        else j++
    }
    if (j==r) return(res)
    return(res[|1 \ j|])
}

end
