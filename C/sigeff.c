/*
** SigEff.c
** Developed by Paul Rayson
** Part of Wmatrix and previously Xmatrix and Tmatrix
**
** Last updated 4july16
**
** See http://ucrel.lancs.ac.uk/llwizard.html for a description of the measures
**
*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#define TRUE            1
#define FALSE           0       
#define NIL            -1
#define EXPECT_CUT_OFF  3               /* cut off value for expectation */

typedef struct
	{
	char	tag[50];
	char	desc[255];
	} tagdesc;

tagdesc	taglist[500];

int main(argc, argv)
int     argc;
char    *argv[];
{
FILE	*ifp;
float	x, y, z;
char	c, d, text[1000];
char	s1[100], s2[100];
char	effectSize[100];			/* effect size type from command line */
int		s1len, i,
		df = 1;						/* degrees of freedom, always 1 in pairwise comparison */
int		num_desc = 0;
int		ln, lineno = 0, chi_sq_fail;
float	chi_sq_val, sum_obs, 
		sum_tot,					/* sum of total frequencies of the corpora */
		likelihood;
float	totf1, totf2, totf3,		/* total frequencies */
		relfreq1, relfreq2,			/* relative frequencies */
		normfreq1, normfreq2,		/* normalised frequencies */
		expect1, expect2, expect3,	/* expected frequencies */
		minexpect;					/* minimum expected frequency */
float	chisq1, chisq2, chisq3,
		like1, like2, like3,
		ivalue, iloglikelihood;
float	top, bottom;				/* temporary variables for calculations */
float	diff,						/* %DIFF measure */
		bayes,						/* Bayes Factor BIC */
		ell,						/* Effect size measure for LL */
		rrisk,						/* relative risk */
		logRatio,					/* Log Ratio */
		oddsRatio;					/* Odds Ratio */
float	zeroCorrectionDiff = 0.000000000000000001,	/* correction for zero in DIFF calc */
		zeroCorrectionLogRatio = 0.5;	/* correction for zero in Log Ratio calc */

/* input lines should be of the form:
TOTAL     100000  100000 
this         231     388 
that         164     355 
these         22     148 
those         51      52 
*/

	/* grab effect size type from command line */
	strcpy(effectSize, argv[1]);

    /* open tagsummary file */
    ifp = fopen(argv[2], "r");
    while (!feof(ifp))
        {
        /* read tag */
        fscanf(ifp, "%s ", s1);
        if (feof(ifp)) break;
        strcpy(taglist[num_desc].tag, s1);
        /* read description */
        fgets(s1, 100, ifp);
        /* get rid of eoln */
        s1len = strlen(s1);
        if (s1[s1len-1] == '\n') s1[s1len-1] = '\0';
        strcpy(taglist[num_desc].desc, s1);
        num_desc++;
        }
    fclose(ifp);
/*    for (i = 0; i < num_desc; i++)
        {
        printf("%-10s %s\n", taglist[i].tag, taglist[i].desc);
        }
*/
	/* read total line at the start */
    fscanf(stdin, "%s %f %f", text, &totf1, &totf2);
    /* write back out totals and header line */
    fprintf(stdout, "%-10s %6.0f %6.0f \n", text, totf1, totf2);
    fprintf(stdout, "%-10s %6s %8s %6s %8s   %8s", "Item", "O1", "%1", "O2", "%2", "LL");
    if (!strcmp(effectSize, "-1"))
	    fprintf(stdout, " %8s", "%DIFF");
	else if (!strcmp(effectSize, "-2"))
	    fprintf(stdout, " %8s", "Bayes");
	else if (!strcmp(effectSize, "-3"))
	    fprintf(stdout, " %11s", "ELL");
	else if (!strcmp(effectSize, "-4"))
	    fprintf(stdout, " %8s", "RRisk");
	else if (!strcmp(effectSize, "-5"))
	    fprintf(stdout, " %8s", "LogRatio");
	else if (!strcmp(effectSize, "-6"))
	    fprintf(stdout, " %8s", "OddsRatio");
	else
	    fprintf(stdout, " %8s %8s %11s %8s %8s %8s", "%DIFF", "Bayes", "ELL", "RRisk", "LogRatio", "OddsRatio");
    fprintf(stdout, "\n");


	/* now read input and write output line by line */
    while (!feof(stdin))
        {
        fscanf(stdin, "%s %f %f", text, &x, &y);
		if (feof(stdin)) continue;
        ln = 0;
/*        c = getc(stdin);
        while (c != '\n' && c != EOF)
            {
            text[ln] = c; ln++;
            c = getc(stdin);
            }
        text[ln] = '\0';
*/
        lineno++;

		/* initialise stuff */
        sum_obs = sum_tot = 0;
		chi_sq_val = 0; chisq1 = chisq2 = chisq3 = 0;
        likelihood = 0; chi_sq_fail = FALSE;
/*        totf1 = totf2 = totf3 = 100000;*/

		/* calculate line totals */
        sum_obs = x + y;
        sum_tot = totf1 + totf2;

		/* calculate expected frequencies */
        expect1 = totf1 * sum_obs / sum_tot;
        expect2 = totf2 * sum_obs / sum_tot;
        /* calculate relative frequencies (percentage) */
		relfreq1 = x * 100 / totf1 ;
		relfreq2 = y * 100 / totf2;
        /* calculate normalised frequencies */
		normfreq1 = x / totf1 ;
		normfreq2 = y / totf2;
		
		/* decide if overused or underused */
		if (relfreq1 < relfreq2)
			d = '-';
		else
			d = '+';
			
		/* legacy chi-squared code */
        if (expect1 < EXPECT_CUT_OFF || expect2 < EXPECT_CUT_OFF)
            {
            chi_sq_fail = TRUE;
            }
        else
            {
            chisq1 = (x - expect1) * (x - expect1) / expect1;
            chi_sq_val += chisq1;
            chisq2 = (y - expect2) * (y - expect2) / expect2;
            chi_sq_val += chisq2;
            }
            
        /* log likelihood calculation */
        /* calculate Oi ln(Oi/Ei) */
        if (x != 0)
            like1 = x * log(x / expect1);
        else
            like1 = 0;
        if (y != 0)
            like2 = y * log(y / expect2);
        else
            like2 = 0;
        /* sum these */
        likelihood += like1 + like2;
        iloglikelihood = 2 * likelihood;
		/* check for very small negative numbers (rounding errors) */
		if (iloglikelihood < 0)
			iloglikelihood = 0;
		if (chi_sq_fail) chi_sq_val = -1;

		/* calculate %DIFF */
		if (normfreq2 == 0)
			diff = (normfreq1 - normfreq2) * 100 / zeroCorrectionDiff ;
		else
			diff = (normfreq1 - normfreq2) * 100 / normfreq2 ;

		/* calculate Bayes Factor BIC */
		bayes = iloglikelihood - (df * log(sum_tot));
		
		/* Calculate effect size measure for LL */
		minexpect = expect1;
		if (expect2 < expect1)
			minexpect = expect2;
		ell = iloglikelihood / (sum_tot * log(minexpect)) ;
		
		/* Calculate Relative Risk */
		rrisk = normfreq1 / normfreq2 ;
		
		/* Calculate Log Ratio */
		/* if either value is zero then use normalised version of zero correction value */
		if (normfreq1 == 0)
			top = zeroCorrectionLogRatio / totf1;
		else
			top = normfreq1;
		if (normfreq2 == 0)
			bottom = zeroCorrectionLogRatio / totf2;
		else
			bottom = normfreq2;
		logRatio = log2(top / bottom);
/*		if (normfreq2 == 0)
			logRatio = log2(normfreq1 / zeroCorrectionLogRatio) ;
		else
			logRatio = log2(normfreq1 / normfreq2) ;
*/

		/* Calculate Odds Ratio */
		oddsRatio = (x/(totf1-x)) / (y/(totf2-y)) ;

		/* add descriptions if input is a semantic tag frequency list */
        /* find matching description - first try exact match */
        strcpy(s1, text);
        s1len = strlen(s1);
        strcpy(s2, "");
        for (i = 0; i < num_desc; i++)
            {
            if (!strcmp(s1, taglist[i].tag))
                {
                strcpy(s2, taglist[i].desc);
                break;
                }
            }

        /* find matching description by removing all but one of plus and minus chars */
		if (s2[0] == '\0')
			{
	        strcpy(s1, text);
	        s1len = strlen(s1);
	        /* remove +/- from the end of the tag */
	        i = s1len - 1;
        	while (s1[i] == '+' || s1[i] == '-')
        	    {
        	    i--;
        	    if (i == 0) break;
        	    }
        	i++;
		if (i < s1len) i++;
        	s1[i] = '\0';
        	strcpy(s2, "");
        	for (i = 0; i < num_desc; i++)
        	    {
        	    if (!strcmp(s1, taglist[i].tag))
        	        {
        	        strcpy(s2, taglist[i].desc);
        	        break;
        	        }
        	    }
		}

        /* find matching description by removing plus and minus chars */
		if (s2[0] == '\0')
			{
	        strcpy(s1, text);
	        s1len = strlen(s1);
	        /* remove +/- from the end of the tag */
	        i = s1len - 1;
        	while (s1[i] == '+' || s1[i] == '-')
        	    {
        	    i--;
        	    if (i == 0) break;
        	    }
        	i++;
        	s1[i] = '\0';
        	strcpy(s2, "");
        	for (i = 0; i < num_desc; i++)
        	    {
        	    if (!strcmp(s1, taglist[i].tag))
        	        {
        	        strcpy(s2, taglist[i].desc);
        	        break;
        	        }
        	    }
			}


		/* finally, write the output line */
        if (c != EOF)
        	{
        	if (!strcmp(effectSize, "-1"))
	            fprintf(stdout, "%-10s %6.0f %8.2f %6.0f %8.2f %c %8.2f %8.2f     %s\n", text, x, relfreq1, y, relfreq2, d, iloglikelihood, diff, s2);
	        else if (!strcmp(effectSize, "-2"))
	            fprintf(stdout, "%-10s %6.0f %8.2f %6.0f %8.2f %c %8.2f %8.2f     %s\n", text, x, relfreq1, y, relfreq2, d, iloglikelihood, bayes, s2);
	        else if (!strcmp(effectSize, "-3"))
	            fprintf(stdout, "%-10s %6.0f %8.2f %6.0f %8.2f %c %8.2f %11.5f     %s\n", text, x, relfreq1, y, relfreq2, d, iloglikelihood, ell, s2);
	        else if (!strcmp(effectSize, "-4"))
	            fprintf(stdout, "%-10s %6.0f %8.2f %6.0f %8.2f %c %8.2f %8.2f     %s\n", text, x, relfreq1, y, relfreq2, d, iloglikelihood, rrisk, s2);
	        else if (!strcmp(effectSize, "-5"))
	            fprintf(stdout, "%-10s %6.0f %8.2f %6.0f %8.2f %c %8.2f %8.2f     %s\n", text, x, relfreq1, y, relfreq2, d, iloglikelihood, logRatio, s2);
	        else if (!strcmp(effectSize, "-6"))
	            fprintf(stdout, "%-10s %6.0f %8.2f %6.0f %8.2f %c %8.2f %8.2f     %s\n", text, x, relfreq1, y, relfreq2, d, iloglikelihood, oddsRatio, s2);
	        else
	            fprintf(stdout, "%-10s %6.0f %8.2f %6.0f %8.2f %c %8.2f %8.2f %8.2f %11.5f %8.2f %8.2f %8.2f     %s\n", text, x, relfreq1, y, relfreq2, d, iloglikelihood, diff, bayes, ell, rrisk, logRatio, oddsRatio, s2);
            }
        }
        
    /* exit */
	return 0;

}	/* end of main */


