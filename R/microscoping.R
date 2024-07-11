# This is package microscoping 

".onLoad" <-
function( libname, pkgname){}


"ckmr_laugh_test" <-
function(
  nata,
  wata,
  pmatata,
  samp_pata,
  samp_ny,
  MAX_JUVE_AGE,
  MIN_CATCH_CURVE_AGE=0,
  beqy_POP=FALSE,
  etriv = 0.01,
  extraplus= FALSE,
  want_gory_details=FALSE
){
## How many POPs, HSPs, and GGPs might we expect ?
## Deliberately ruffff

  # Assume steady-state: Z_a = catch-curve slope at age a...
  # ... but may not be equal across ages
  # Try to give it a bit of flexibility

  # Replace actual n[a] with smoothed version (since we assume steady state);
  # expand (possible) plus-group up to "ancient" age-class beyond which lies trivia

  # Fit catch curve, omitting final age-class in case it's a plus-group

  MJA <- MAX_JUVE_AGE
  stopifnot( my.all.equal( length( nata), length( wata), length( pmatata), length( samp_pata)))

  # I want everything to start at age 0
  atathings <- cq( nata, wata, pmatata, samp_pata)
  # Check they're all the same length, with some "ingenious" code...
stopifnot( length( unique( lengths( mget( atathings)))) == 1)

  for( ithing in atathings) {
    thing <- get( ithing)
    if( thing %is.not.an% 'offarray') {

      assign( ithing, offarray( thing, first=c( AGE=0), dim=length( thing)))
    }
  }

  Aplus <- dimrange( nata)[2]
  Apreplus <- Aplus - 1
  anonplus <- 0:Apreplus

  direct_TRO <- sum( wata * pmatata * nata) # Get this directly, not post-smoothing
  nata[ VECSUB= which(nata==0)] <- 1e-3 # one thousandth of an animal !
  smoon_a <- offarray( -1, dimseq=list( 0:Aplus))
  smoon_a[ anonplus] <- smooth_n( nata[ MIN_CATCH_CURVE_AGE %upto% Apreplus])

  # Plus-group--- for now, all stuck together
  # NB SLICE to avoid single-element offarrays and dimseq mismatches...
  if( extraplus) {
    # "Real" plus-group ignored, and just extrapolated from preplus
    ps <- smoon_a[ SLICE=Apreplus] / smoon_a[ SLICE=Apreplus-1]
    smoon_a[ Aplus] < smoon_a[ SLICE=Apreplus] * ps / (1-ps)
  } else {
    # Figure out (pre)plus-group-surv that would spread "real" plus-group to line up with preplus
    ps <- nata[ SLICE=Aplus] / (smoon_a[ SLICE=Apreplus] + nata[ SLICE=Aplus])
    smoon_a[ Aplus] <- nata[ SLICE=Aplus] # will get smeared out across expanded ages
  }

  # Need Age-extrap range above which fec contrib is neglig (<= etriv*total)

  totfec_preplus <- sum( wata[ anonplus] * pmatata[ anonplus] * smoon_a[ anonplus])
  fec_plus <- c( wata[ Aplus] * pmatata[ Aplus] * smoon_a[ Aplus])

  ancient <- log( etriv * (1 + totfec_preplus/fec_plus)) / log( ps)
  ancient <- floor( 1+max( ancient, 0)) # indeed, ancient *should* not be < 0, but...

  AMAX <- Aplus + ancient
  ages <- 0:AMAX

  # Extending an offarray does take extra effort--- tedious
  template <- offarray( -1, first=0L, last=AMAX)
  template[ anonplus] <- wata[ anonplus]
  template[ Aplus + 0:ancient] <- tail( wata, 1)
  wata <- template

  template[ anonplus] <- pmatata[ anonplus]
  template[ Aplus + 0:ancient] <- tail( pmatata, 1)
  pmatata <- template

  template[ anonplus] <- smoon_a[ anonplus]
  template[ Aplus + 0:ancient] <-
      smoon_a[ SLICE=Aplus] * (1-ps) * (ps ^ (0:ancient))
  nata <- template

  ny <- length( samp_ny)
  years <- 1:ny

## Actual CKMR stuff begins here...

  # Just do female RO and female-female comps for now; fiddle the numbers at the end to allow for males
  nfata <- nata / 2 # females
  fecata <- wata * pmatata
  TROF <- nfata %**% fecata
  inv_TROF <- 1/TROF

  fec_fn <- function( a) {
      # return 0 unless 0 <= a < AMAX
      # must be vectorized, for autoloop, and so must clamp index range
      # NOOFF since 'a' won't be consec ints
      # fecata is zero-based
      (a > 0) * fecata[ NOOFF= a |> clamp( 0, AMAX-1)] 
    }

  ## POPs
  extract.named( autoloop( A1=ages, Y1=years, A2=ages, Y2=years, {
    B1 <- Y1-A1
    B2 <- Y2-A2

    # Fish 1 the parent (Mother) of Fish 2 (Daughter)?
    Pr_M1D2 <- (B1 < B2) * (Y1 >= B2+!beqy_POP) * fec_fn( A1 - (Y1-B2)) * inv_TROF
    # or vice versa
    Pr_D1M2 <- (B2 < B1) * (Y2 >= B1+!beqy_POP) * fec_fn( A2 - (Y2-B1)) * inv_TROF

    Pr_MDP <- Pr_M1D2 + Pr_D1M2 # Mother/daughter (or vice versa)
    Pr_POP <- Pr_MDP   # This is *only* for female-female comp, but I'll call it POP for now
  returnList( Pr_POP,Pr_MDP, Pr_M1D2, Pr_D1M2)
  }))

  # HSPs and GGPs: do it conditional on age of "unobserved parent"
  # Then weighted sum over prob of AgeOfUnobsPar

  Pr_Mum1Age <- offarray( (nfata * fecata) *  inv_TROF, first=0, last=AMAX) # unobserved parent

  psurv_ay <- get_psurv( nata, ny+2*AMAX)

  # ALL-FEMALE-ONLY cases (ie boths sibs and their mum are Female, or Grandma/Mum/Grand-daughter)
  extract.named( autoloop( A1=ages, Y1=years, A2=ages, Y2=years, PARAGE=ages, {
    B1 <- Y1-A1
    B2 <- Y2-A2
    Bearly <- pmin( B1, B2)
    Blate <- pmax( B2, B1)
    Bdiff <- Blate - Bearly
    Yearly <- ifelse( B1<B2, Y1, Y2)
    Ylate <- Y1 + Y2 - Yearly
    Aearly <- ifelse( B1<B2, A1, A2)
    Alate <- A1 + A2 - Aearly

    Pr_HSP_cond <- psurv_ay[ PARAGE, Bdiff] * fec_fn( PARAGE+Bdiff) * inv_TROF

    # GGP where missing middle is Mum
    Bmid <- Blate - PARAGE
    Pr_GGP_cond <- (Yearly > Bmid) * fec_fn( Aearly - (Yearly-Bmid)) * inv_TROF

    # Anti-double-count goes into ncomps...
    # but anti-same-cohort goes here:
    Pr_GGP_cond <- Pr_GGP_cond * (Bearly != Blate)
    Pr_HSP_cond <- Pr_HSP_cond * (Bearly != Blate) # surely 0 anyway !

    # Kin-probs are *conditional* on Mum's age
    # Make it joint prob
    Pr_HSP_and_Mum1Age <- Pr_HSP_cond * Pr_Mum1Age[ PARAGE]
    Pr_GGP_and_Mum1Age <- Pr_GGP_cond * Pr_Mum1Age[ PARAGE]
  returnList( Pr_GGP_and_Mum1Age, Pr_HSP_and_Mum1Age) # ie these two are extracted
  }))

  # DOUBLE because the unobserved "mother" could be a "father" !
  # Not applicable for POPs
  Pr_HSP <- 2*sumover( Pr_HSP_and_Mum1Age, 'PARAGE')
  Pr_GGP <- 2*sumover( Pr_GGP_and_Mum1Age, 'PARAGE')


  # Omit when adult catch year == juve birth year; this could be changed, eg if adult catches are
  # ... mainly seasonal and *after* the spawning season within the calendar year
  # ... But, if you are worrying deeply about such details, you should anyway be doing Detailed Design instead ;)
  # Avoid double-counts
  # Disallow same-cohort HSPs
  # Omit possible GGPs

  # *Female* samples per year and age
  # fsamp_ya[ y, a] := 0.5 * samp_ny[ y] * samp_pata[ a]

  # samp_pata needs to expand thru plus-group
  # split the plus-group equally
  nextra <- length( ages) - length( samp_pata)
  # Next bit used to be wrong...
  samp_pata <- offarray( c( head( samp_pata, -1),
      rep( tail( samp_pata, 1) / (nextra+1), (nextra+1))),
      first=0, last=AMAX)
  samp_pata <- samp_pata / sum( samp_pata)

  fsamp_ya <- autoloop( y=years, a=ages, {
    0.5 * samp_ny[ y] * samp_pata[ a]
  })


  # Sequence fish by year then age, and enforce earlier-compared-with-later
  # For same-year same-age comps, that means only half the number (to avoid double-count)
  ncomps <- autoloop( A1=ages, Y1=years, A2=ages, Y2=years, {

    # Avoid double-counting; Y1<Y2, or if equal then A1<A2, or
    # if both equal, then only half the comps
    ppn_to_use <- ifelse(
      (Y1 < Y2) | ((Y1==Y2) & (A1 < A2)),
        1, ifelse(
      (A1==A2) & (Y1==Y2),
        0.5,
        0
    ))
    fsamp_ya[ Y1, A1] * fsamp_ya[ Y2, A2] * ppn_to_use
  })

  #    ( ((Y1 < Y2) | (A1 < A2)) + 0.5 * (Y1==Y2) * (A1==A2))

  # We're not just sampling Females; we get to do F/F, M/F, F/M, and M/M comparisons
  # So...
  ncomps <- ncomps * 4

  nkins <- autoloop( A1=ages, Y1=years, A2=ages, Y2=years, {
    # Distinguish POPs according as *offspring* is above MJA
    POP_offJ <- Pr_M1D2[A1,Y1,A2,Y2] * ncomps[A1,Y1,A2,Y2] * (A2 <= MJA) +
        Pr_D1M2[A1,Y1,A2,Y2] * ncomps[A1,Y1,A2,Y2] * (A1 <= MJA)
    POP_offA <- Pr_M1D2[A1,Y1,A2,Y2] * ncomps[A1,Y1,A2,Y2] * (A2 > MJA) +
        Pr_D1M2[A1,Y1,A2,Y2] * ncomps[A1,Y1,A2,Y2] * (A1 > MJA)

    #POP_AJ <- Pr_POP * ncomps * (A1 <= MJA) * (A2 > MJA)
    #POP_AA <- Pr_POP * ncomps * (A1 > MJA) * (A2 > MJA)
    HSP_JJ <- Pr_HSP[A1,Y1,A2,Y2] * ncomps[A1,Y1,A2,Y2] * (A1 <= MJA) * (A2 <= MJA)
    GGP_JJ <- Pr_GGP[A1,Y1,A2,Y2] * ncomps[A1,Y1,A2,Y2] * (A1 <= MJA) * (A2 <= MJA)
  returnList( POP_offJ, POP_offA, HSP_JJ, GGP_JJ)
  })

  results <- lapply( nkins, sum)
  # extract.named( results) # if doing further processing

  # Total number of samples collected
  results$nsamp_J <- sum( 2*fsamp_ya[,0 %upto% MJA])
  results$nsamp_A <- sum( 2*fsamp_ya[,(MJA+1) %upto% AMAX])
  mode( results) <- 'numeric' # as.numeric() kills names

  results@call <- sys.call()
  params <- as.list( mget( names( formals( sys.function()))))
  class( params) <- 'nullprint'
  results@params <- params

  if( want_gory_details) {
    results@gory <- environment()
  }

return( results)
}


"clt2" <-
function(
  n_sa,
  w_sa,
  pmat_sa,
  nsamp_sya= NULL, # if null, construct from next two
  samp_p_sa,
  nsamp_y,
  MAX_JUVE_AGE,
  MIN_CATCH_CURVE_AGE= 0,
  beqy_POP= FALSE,
  etriv= 0.01,
  extraplus= FALSE, # ie use "real" input plus-group
  want_gory_details= FALSE,
  .BLOCKSIZE= 1e4,
  .REPORTSEC= 5
){
## How many POPs, HSPs, and GGPs might we expect ?
## Deliberately ruffff

  # Convert all the blah_sa args to offarray
  # and tedious sanity checks on dimensions
  # Too tedious to be thorough, so it's basiscally up to the user
  sa_things <- cq( n, w, pmat) %&% '_sa'

  if( !missing( samp_p_sa)) {
    sa_things <- c( sa_things, 'samp_p_sa')
  }

stopifnot( {
      x <- FOR( sa_things, dim( get(.)))
      all( do.on( x, all( .==x[[1]])))
    },
    nrow( n_sa)==2)

  # Check dimnames[[1]], ie Sex, consistent
  # and ensure names( dimnames) is (SEX, AGE)
  sdn1 <- FOR( sa_things, dimnames( get( .))[[1]])
  nnsdn1 <- sdn1 %SUCH.THAT% !is.null( .)

  SEXES <- cq( F, M)
  if( length( nnsdn1)) { # must all be the same
    s1 <- toupper( substring( nnsdn1[[1]], 1, 1)) # F,M or M,F
stopifnot( sort( s1) == SEXES,
      all( do.on( nnsdn1, all( .==nnsdn1[[ 1]]))) )
  } else {
    s1 <- SEXES # none are set; will just *assume* F=1, M=2
    # ... doesn't matter for these calcs
  }

  for( sa_thing in sa_things) {
    i <- get( sa_thing)
    if( !all( dimnames( i)[[1]] == s1)) {
      dimnames( i) <- list( SEX= s1, dimnames( i)[-1])
      i <- i[ SEXES,] # ensure F first; might not be, in original data
    }
    if( i %is.not.an% 'offarray') {
      i <- offarray( c( i), first=c( 1, 0), last=c( 2, Aplus),
          dimnames=list( SEX=SEXES, AGE=NULL))
    }

    names( dimnames( i)) <- cq( SEX, AGE) # AGE is used in glm formula
    assign( sa_thing, i)
  }

  Aplus <- dimrange( n_sa)['AGE',2]
  Apreplus <- Aplus-1

  ## At one point I had code for sex-specific MAX_JUVE_AGE, but complications...
  ## ... so now it's one MJA for both
  if( FALSE) {
    # In the one case of MJA, we'll permit length-1 (rather than sex-specific) arg...
    if( length( MAX_JUVE_AGE)==1) {
      MAX_JUVE_AGE <- rep( unname( MAX_JUVE_AGE), 2)
    }
    if( is.null( names( MAX_JUVE_AGE))) { # then assume same sex-order as n_sa, w_sa, etc
      names( MAX_JUVE_AGE) <- s1
    }
    MJA <- MAX_JUVE_AGE[ SEXES]  # correct order
  }
  MJA <- MAX_JUVE_AGE
stopifnot( length( MJA)==1)

  direct_TRO <- rowSums( w_sa * pmat_sa * n_sa) # Get this directly, not post-smoothing

  n_sa[ VECSUB= (n_sa==0)] <- 1e-3 # one thousandth of an animal !
  smoon_sa <- offarray( -1, first=c( 1, 0), last=c( 2, Aplus), dimnames=list( SEX=SEXES, AGE=NULL))

  plusgroup_surv <- ancient <- structure( c( -1, -1), names=SEXES) # by sex

  anonplus <- 0:Apreplus
  for( sex in SEXES) {
    smoon_sa[ sex, anonplus] <- smooth_n( n_sa[ SLICE=sex, MIN_CATCH_CURVE_AGE %upto% Apreplus])

    # Plus-group--- for now, all stuck together
    if( extraplus) {
      # "Real" plus-group ignored, and just extrapolated from preplus
      ps <- plusgroup_surv[ sex] <- c( smoon_sa[ sex, Apreplus] / smoon_sa[ sex, Apreplus-1])
      smoon_sa[ sex, Aplus] < smoon_sa[ sex, Apreplus] * ps / (1-ps)
    } else {
      # Figure out (pre)plus-group-surv that would spread "real" plus-group to line up with preplus
      ps <- plusgroup_surv[ sex] <- n_sa[ sex, Aplus] / (smoon_sa[ sex, Apreplus] + n_sa[ sex, Aplus])
      smoon_sa[ sex, Aplus] <- n_sa[ sex, Aplus] # will get smeared out across expanded ages
    }

    # Need Age-extrap range by sex, then pick the bigger
    # IE ge (post-plus) above which tot fec contrib is neglig

    # Should use SLICE= here but too lazy
    totfec_preplus <- sum( w_sa[ sex, anonplus] * pmat_sa[ sex, anonplus] * smoon_sa[ sex, anonplus])
    fec_plus <- c( w_sa[ sex, Aplus] * pmat_sa[ sex, Aplus] * smoon_sa[ sex, Aplus])


    ancient[ sex] <- log( etriv * (1 + totfec_preplus/fec_plus)) / log( ps)
    ancient[ sex] <- floor( 1+max( ancient[ sex], 0)) # indeed, ancient *should* not be < 0, but...
  } # for sex

  max_extra <- max( ancient) # break up plus-group into how many age-classes?
  AMAX <- Aplus + max_extra
  AGES <- 0:AMAX
  xthing <- offarray( 0, first=c( SEX=1, AGE=0), last=c( 2, AMAX),
      dimnames=list( SEX=SEXES, AGE=NULL))
  for( thing in cq( w_sa, pmat_sa)) {
    xthing[,0:Aplus] <- get(thing)
    xthing[,Aplus + (1 %upto% max_extra)] <- xthing[, rep( Aplus, max_extra), NOOFF=TRUE]
    assign( thing, xthing)
  }
  fec_sa <- w_sa * pmat_sa
  n_sa <- 0 * pmat_sa

  # Spread out plus-group across extra ages
  # Actual number exactly at Aplus is plus-group-total*(1-surv)
  for( sex in SEXES) {
    n_sa[ sex,] <- c( smoon_sa[ SLICE=sex, anonplus],
        smoon_sa[ SLICE=sex, SLICE=Aplus] *
        (1-plusgroup_surv[ sex]) *
        (plusgroup_surv[ sex] ^ (0 %upto% max_extra) ))
  }

  ny <- if( is.null( nsamp_sya)){
      length( nsamp_y)
    } else {
      dim( nsamp_sya)[2]
    }
  YEARS <- 1:ny

## Actual CKMR stuff begins here...

  TRO <- sumover( n_sa * fec_sa, 'AGE') # still by sex
  # ... now done directly from input n
  inv_TRO <- 1 / TRO

  fec_fn <- function( s, a) {
      # check 0 <= a < AMAX
      # If not, change a to something that allows fec_sa[ awork] to succeed...
      awork <- pmax( 0, pmin( AMAX-1, a))
      # ... and set fec to 0 if a <= 0
      # fec_sa is zero-based
      # This will be called inside autoloop, so need "autovec" version
      fec_sa[ MATSUB=data.frame( s, awork)] * (a > 0)
    }

  ## POPs
  extract.named( autoloop( S1=SEXES, A1=AGES, Y1=YEARS, S2=SEXES, A2=AGES, Y2=YEARS, {
    B1 <- Y1-A1
    B2 <- Y2-A2

    # Fish 1 the Parent of Fish 2 the Offspring ?
    Pr_P1O2 <- (B1 < B2) * (Y1 >= B2+!beqy_POP) * fec_fn( S1, A1 - (Y1-B2)) * inv_TRO[ S1]
    # or vice versa
    Pr_O1P2 <- (B2 < B1) * (Y2 >= B1+!beqy_POP) * fec_fn( S2, A2 - (Y2-B1)) * inv_TRO[ S2]

    Pr_POP <- Pr_P1O2 + Pr_O1P2
  returnList( Pr_POP, Pr_P1O2, Pr_O1P2)
  }))

  # HSPs and GGPs: do it conditional on age of "unobserved parent"
  # Then weighted sum over prob of AgeOfUnobsPar

  Pr_Par1Age <- autoloop( S=SEXES, A=AGES, {
    n_sa[ S, A] * fec_sa[ S, A] * inv_TRO[ S]
  })

  # Survival: leave room for old parents of HSPs
  psurv_say <- autoloop( S=SEXES, A=AGES, DY=0 %upto% (ny+2*AMAX), 0)
  for( S in SEXES) {
    psurv_say[ S,,] <- get_psurv( n_sa[ SLICE=S,], ny+2*AMAX)
  }

  # extract.named( autoloop( S1=SEXES, A1=AGES, Y1=YEARS, S2=SEXES, A2=AGES, Y2=YEARS, PARSEX=SEXES, PARAGE=AGES, {
  extract.named( autoloop( S1=SEXES, A1=AGES, Y1=YEARS, S2=SEXES, A2=AGES, Y2=YEARS,
      SUMOVER=list( PARSEX=SEXES, PARAGE=AGES),
      .BLOCKSIZE=.BLOCKSIZE,
      .REPORTSEC=.REPORTSEC, {
    B1 <- Y1-A1
    B2 <- Y2-A2
    Bearly <- pmin( B1, B2)
    Blate <- pmax( B2, B1)
    Bdiff <- Blate - Bearly
    Yearly <- ifelse( B1<B2, Y1, Y2)
    Ylate <- Y1 + Y2 - Yearly
    Aearly <- ifelse( B1<B2, A1, A2)
    Alate <- A1 + A2 - Aearly
    Searly <- ifelse( B1<B2, S1, S2) # for GGP only
    # Slate doesn't matter for either

    # For both HSP & POP:
    # Pr[ K & PARAGE | PARSEX] where PAR is parent of Later
    # = Pr[ K | PARAGE, PARSEX] * Pr[ PARAGE | PARSEX]
    Pr_HSP <-
        Pr_Par1Age[ PARSEX, PARAGE] *
        psurv_say[ PARSEX, PARAGE, Bdiff] *
        fec_fn( PARSEX, PARAGE+Bdiff) *
        inv_TRO[ PARSEX]

    # Anti-double-count goes into ncomps...
    # but anti-same-cohort goes HERE:
    Pr_HSP <- Pr_HSP * (Bearly != Blate)

    # GGP --- needs Searly...
    Bmid <- Blate - PARAGE
    Pr_GGP <-
        Pr_Par1Age[ PARSEX, PARAGE] *
        (Yearly > Bmid) *
        fec_fn( Searly, Aearly - (Yearly-Bmid)) *
        inv_TRO[ Searly]

  returnList( Pr_GGP, Pr_HSP) # ie these two are extracted
  }))


  # Omit when adult catch year == juve birth year; this could be changed, eg if adult catches are
  # ... mainly seasonal and *after* the spawning season within the calendar year
  # ... But, if you are worrying deeply about such details, you should anyway be doing Detailed Design instead ;)


  if( is.null( nsamp_sya)) {
    # Build from bits
    samp_p_sa <- samp_p_sa / sum( samp_p_sa) # just in case...
    nsamp_sya <- autoloop( S=SEXES, Y=YEARS, A=dimseq( samp_p_sa)[[2]], {
        nsamp_y[ Y] * samp_p_sa[ S, A]
      })
  }

  # nsamp_sya plus-group needs to be split up
  if( !my.all.equal( dimseq( nsamp_sya)[[3]], AGES)) {
    # xthing already has right shape, thx2 w_sa
    ex_nsamp_sya <- autoloop( 0, S=SEXES, Y=YEARS, A=AGES) # just set dims
    for( Y in YEARS) {
      xthing[] <- 0
      xthing[,0:Apreplus] <- nsamp_sya[,SLICE=Y, 0:Apreplus]
      xthing[,Aplus:AMAX] <- nsamp_sya[,SLICE=Y, rep( Aplus, max_extra+1), NOOFF=TRUE] / (max_extra+1)
      ex_nsamp_sya[,Y,] <- xthing
    }
    nsamp_sya <- ex_nsamp_sya
  }


  # Sequence fish by sex then year then age, and enforce earlier-compared-with-later
  # For same-year same-age comps, that means only half the number (to avoid double-count)
  ncomps <- autoloop( S1=SEXES, A1=AGES, Y1=YEARS, S2=SEXES, A2=AGES, Y2=YEARS, {
    # Avoid double-counting; S1<S2, or if equal then Y1<Y2, or if also equal then A1<A2, or
    # if all 3 vars equal, then only half the comps

    ppn_to_use <- ifelse( (
      (S1 < S2) | ((S1==S2) & (
      (Y1 < Y2) | ((Y1==Y2) &
      (A1 < A2))))),
        1, ifelse(
      (S1==S2) & (A1==A2) & (Y1==Y2),
        0.5,
        0
    ))

    nsamp_sya[ S1, Y1, A1] * nsamp_sya[ S2, Y2, A2] * ppn_to_use
  })

  nkins <- autoloop( S1=SEXES, A1=AGES, Y1=YEARS, S2=SEXES, A2=AGES, Y2=YEARS, {
    # Distinguish POPs according as *offspring* is above MJA
    # NB the parent can always be a "juvenile" ie below MJA
    POP_offJ <- (
        Pr_P1O2[ S1, A1, Y1, S2, A2, Y2] * (A2 <= MJA) +
        Pr_O1P2[ S1, A1, Y1, S2, A2, Y2] * (A1 <= MJA)) *
        ncomps[ S1, A1, Y1, S2, A2, Y2]
    POP_offA <- (
        Pr_P1O2[ S1, A1, Y1, S2, A2, Y2] * (A2 > MJA) +
        Pr_O1P2[ S1, A1, Y1, S2, A2, Y2] *  (A1 > MJA)) *
        ncomps[ S1, A1, Y1, S2, A2, Y2]
    # POP_all <- Pr_POP * ncomps # debugging check...

    HSP_JJ <-
        Pr_HSP[ S1, A1, Y1, S2, A2, Y2] *
        ncomps[ S1, A1, Y1, S2, A2, Y2] *
        (A1 <= MJA) *
        (A2 <= MJA)
    # ... NB you might define MJA > Amat if you're sure of age, so J may not be immature...
    # ... so GGP_JJ may not be negligible despite its name ...
    # ... in fact that's one basis for _setting_ MJA (to a value s.t. GGP is neglig)
    GGP_JJ <-
        Pr_GGP[ S1, A1, Y1, S2, A2, Y2] *
        ncomps[ S1, A1, Y1, S2, A2, Y2] *
        (A1 <= MJA) *
        (A2 <= MJA)

  returnList( POP_offJ, POP_offA, HSP_JJ, GGP_JJ)
  })

  results <- lapply( nkins, sum)

  # Total number of samples collected
  results$nsamp_J <- sum( nsamp_sya[,, 0 %upto% MJA])
  results$nsamp_A <- sum( nsamp_sya[,, (MJA+1) %upto% AMAX])

  mode( results) <- 'numeric' # as.numeric() kills names

  results@call <- sys.call()
  params <- as.list( mget( names( formals( sys.function()))))
  class( params) <- 'nullprint'
  results@params <- params

  oldClass( nsamp_sya) <- c( 'nullprint', oldClass( nsamp_sya))
  results@nsamp_sya <- nsamp_sya

  if( want_gory_details) {
    results@gory <- environment()
  }


return( results)
}


"get_psurv" <-
function( n, DYMAX) {
  # n is offarray presumably
  # Anyway, should start at age 0
  AMAX <- length( n) - 1
  BIGA <- 3*AMAX + DYMAX
  
  lpsurv <- diff( log( unclass( n)))
  lpsurv <- c( 0, lpsurv, rep( lpsurv[ AMAX-1], BIGA - AMAX))
  clpsurv <- offarray( cumsum( lpsurv), first=0, dim=length( lpsurv))
  
  # DYMAX = 3*AMAX+ny
  lpsurv_ay <- autoloop( A=0:AMAX, DY=0:DYMAX,
    clpsurv[ DY+A] - clpsurv[ A] # DY more years starting from A
  )

return( exp( lpsurv_ay))
}


"smooth_n" <-
function( nata) {
## Convert n-at-a into smoothly-decreasing sequence, since it's "in equlibrium"
## and so that age-specific eqm surv can be found
## ie fit a catch-curve
## It might be smooth, but it's still ugly...

  # Caller should exclude final age--- presumably/potentially plus-group
  # Try more-flexible models first, and go simpler if they fail
	# Tweedie(1.7) as compromise between Poisson noise and Gamma rct-var
	
	n_df <- as.data.frame( nata, name='n_a') #
	if( nrow( n_df) > 3) {
		# Tweedie(1.7) as compromise between Poisson noise and Gamma rct-var
		catch_curve <- try( glm( n_a ~ AGE + sqrt( AGE), family=Tweedie( p=1.7, link='log'),
				data=n_df), silent=TRUE)
		no_sqrt_term <- (catch_curve %is.a% 'try-error') || !all( diff( fitted( catch_curve)) < 0)
	} else {
		no_sqrt_term <- TRUE
	}

	if( no_sqrt_term){
		# try again, enforcing constant age-specific mortality
		# this can actually go wrong if we go straight to Tweedie with 2 age classes, so total paranoia version
		# quasipoisson in case of non-integer n_a...
		catch_curve <- glm( n_a ~ AGE, family=quasipoisson( link='log'), data=n_df)
		# ... use as startvals for Tweedie version
		test_catch_curve <- try(
				glm( n_a ~ AGE, family=Tweedie( p=1.7, link='log'), data=n_df, mustart=fitted( catch_curve)),
				silent=TRUE)
		if( test_catch_curve %is.not.a% 'try-error') {
			# sometimes Tweedie misbehaves
			catch_curve <- test_catch_curve
		}
	}

	if( !all( diff( fitted( catch_curve, type='response')) < 0)) {
stop( "Can't fit catch-curve with survival < 1")
	}

	newdf <- data.frame( AGE=0:dimrange( nata)[2], n_a=-1)
	smootho <- unname( predict( catch_curve, newdata=newdf, type='response'))

return( smootho)
}

