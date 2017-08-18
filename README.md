
NolteWeisser\_AluminiumKinetics
===============================

The compartmental model for aluminium kinetics in human body of Nolte et al \[1,2\], as reparameterized by Weisser et al \[3\], implemented under [mrgsolve](https://mrgsolve.github.io/) and [R](https://www.r-project.org/).

1.  Nolte E, Beck E, Winklhofer C, Steinhausen C. Compartmental model for aluminium biokinetics. Hum Exp Toxicol. 2001 Feb;20(2):111-7.
2.  Steinhausen C, Kislinger G, Winklhofer C, Beck E, Hohl C, Nolte E, Ittel TH, Alvarez-Brückmann MJ. Investigation of the aluminium biokinetics in humans: a 26Al tracer study. Food Chem Toxicol. 2004 Mar;42(3):363-71.
3.  Weisser K, Stübler S, Matheis W, Huisinga W. Towards toxicokinetic modelling of aluminium exposure from adjuvants in medicinal products. Regul Toxicol Pharmacol. 2017 Aug;88:310-321. doi: 10.1016/j.yrtph.2017.02.018. Epub 2017 Feb 22.

Remarks
-------

The model file contains the mean parameter set — as presented in the supplementary appendix of \[3\] — so this is used in calculations by default. However, it can be changed in R (without changing the model file) as this is only a specification for defaults, and can be overriden.

I tried to double-check everything, but no guarantee of course. I welcome every suggestion or correction!

Example
-------

Necessary libraries:

``` r
library( mrgsolve )
library( lattice )
```

We first load the model:

``` r
mod <- mread( "NolteWeisser" )
```

    ## Compiling NolteWeisser ... done.

Let's first simulate what happens for oral dosing, using the mean parameter set (replicating Figure 4a of the supplementary appendix of \[3\]):

``` r
e <- ev( time = 0, evid = 1, amt = 1, cmt = "STO" )
out <- mrgsim( ev( mod, e ), end = 10000, delta = 0.1 )

xyplot( PLASMA*100 + U*100 + LIV*100 + MUS*100 + BON*100 ~ time, data = out@data[ -1, ], type = "l",
        scale = list( log = TRUE ),
        auto.key = list( text = c( "Plasma", "Urine", "Liver", "Muscle", "Bone" ),
                         lines = TRUE, points = FALSE, columns = 2 ),
        xlab = "Time [h]", ylab = "Aluminium [% ID]" )
```

<img src="README_files/figure-markdown_github-ascii_identifiers/oraldosing-1.png" style="display: block; margin: auto;" />

Iv dosing, using the mean parameter set (replicating Figure 4b of the supplementary appendix of \[3\]):

``` r
e <- ev( time = 0, evid = 1, amt = 1, cmt = "PC" )
out <- mrgsim( ev( mod, e ), end = 10000, delta = 0.1 )

xyplot( PLASMA*100 + U*100 + LIV*100 + MUS*100 + BON*100 + FAE*100 ~ time, data = out@data[ -1, ],
        type = "l", scale = list( log = TRUE ),
        auto.key = list( text = c( "Plasma", "Urine", "Liver", "Muscle", "Bone", "Faeces" ),
                         lines = TRUE, points = FALSE, columns = 2 ),
        xlab = "Time [h]", ylab = "Aluminium [% ID]" )
```

<img src="README_files/figure-markdown_github-ascii_identifiers/ivdosing-1.png" style="display: block; margin: auto;" />
