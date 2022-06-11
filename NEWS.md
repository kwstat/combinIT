combinIT 1.5.0
===========

Updates
-------
-   `Boik_test` updated to have an argument named 'report' to provide users some
     reports on the result of the test and to describe the pattern of interaction. In        addition, the 'alpha' argument was added to perform test at the alpha level.
    
-   `Piepho_test` updated to have an argument named 'report' to provide users some           reports on the result of the test and to describe the pattern of interaction. In        addition, the 'alpha' argument was added to perform test at the alpha level.

-   `KKM_test` updated to have an argument named 'report' to provide users some reports      on the result of the test and to describe the pattern of interaction. In                addition, the 'alpha' argument was added to perform test at the alpha level.

-   `Malik_test` updated to have an argument named 'report' to provide users some            reports on the result of the test and to describe the pattern of interaction. In        addition, the 'alpha' argument was added to perform test at the alpha level.

-   `Franck_test` updated to have an argument named 'report' to provide users some           reports on the result of the test and to describe the pattern of interaction. In        addition, the 'alpha' argument was added to perform test at the alpha level. An         argument named 'plot' was also added to plot the interaction plot and characterize      the hidden interaction structure. Arguments 'vecolor' and 'linetype' are for            customizing the plotted lines in the interaction plot.
 
-   `KKSA_test` updated to have an argument named 'report' to provide users some             reports on the result of the test and to describe the pattern of interaction. In        addition, the 'alpha' argument was added to perform test at the alpha level. An         argument named 'plot' was also added to plot the interaction plot and characterize      the heterogeneous sub-tables. Arguments 'vecolor' and 'linetype' are for                customizing the plotted lines in the interaction plot.

-   `CPI_test` updated to have an argument named 'report' to provide users some              reports on the result of the test and to describe the pattern of interaction. In        addition, the 'alpha' argument was added to perform test at the alpha level. An         argument named 'opvalue' was also added to combine the p-values of some other tests      (in addition to the six considered tests).

New Additions
-------------

-   `Result.Boik`, `Result.Piepho`, `Result.KKM`, `Result.KKSA`, `Result.Franck`, and       `Result.Malik` were added as internal functions to summarize the results of the six      considered interaction tests.

Minor improvements and bug fixes
--------------------------------

-   Fixed some bugs in `CPI_test`, `Boik_test`, `Piepho_test`, and `Franck_test`. Some      typos were corrected in function documents. In addition, `ITtestClass_methods`,         `combtest_methods`, and `plots` were improved. 
