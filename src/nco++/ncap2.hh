/* $Header: /data/zender/nco_20150216/nco/src/nco++/ncap2.hh,v 1.9 2006-03-15 14:01:16 hmb Exp $ */

/* Purpose: netCDF arithmetic processor definitions and function prototypes for ncap.c, ncap_utl.c, ncap_lex.l, and ncap_yacc.y */

/* Copyright (C) 1995--2005 Charlie Zender
   This software may be modified and/or re-distributed under the terms of the GNU General Public License (GPL) Version 2
   See http://www.gnu.ai.mit.edu/copyleft/gpl.html for full license text */

/* Usage:
   #include "ncap.h" *//* netCDF arithmetic processor-specific definitions (symbol table, ...) */

#ifndef NCAP_H /* Header file has not yet been defined in current source file */
#define NCAP_H

#ifdef HAVE_CONFIG_H
#include <config.h> /* Autotools tokens */
#endif /* !HAVE_CONFIG_H */

/* Standard header files */
#include <math.h> /* sin cos cos sin 3.14159 */
#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* atof, atoi, malloc, getopt */
#include <string.h> /* strcmp. . . */
#include <time.h> /* machine time */
#include <unistd.h> /* POSIX stuff */
#include <string>
/* 3rd party vendors */
#include <netcdf.h> /* netCDF definitions and C library */
#include "nco_netcdf.h" /* NCO wrappers for libnetcdf.a */

/* Personal headers */
#include "libnco.h" /* netCDF Operator (NCO) library */
//#include "libnco++.hh" /* netCDF Operator (NCO) C++ library */

// defines custom "template" lists
#include "Ncap2.hh"
#include <vector>
#include "NcapVector.hh"
#include "NcapVarVector.hh"
#include "NcapVar.hh"

#include <antlr/AST.hpp> /* nneeded for ast_ind struct */

/* Define symbol table */

typedef struct{ /* sym_sct */
  char *nm; /* [sng] Symbol name */
  double (*fnc_dbl)(double); /* [fnc] Double-valued function */
  float (*fnc_flt)(float); /* [fnc] Float-valued function */
} sym_sct;

/* Name list structure ncap.c
(for subscript lists) */
typedef struct{ /* nm_lst_sct */
  nm_id_sct *lst; /* [sct] List element */
  int nbr; /* [nbr] Number of structures in list */
} nm_lst_sct;


/* Structure to hold AST pointers to indices in hyperslabs -only temporary */
typedef struct{
  ANTLR_USE_NAMESPACE(antlr)RefAST ind[3];
} ast_lmt_sct;   



typedef struct{ /* prs_sct */
public:
  char *fl_in; /* [sng] Input data file */
  int in_id; /* [id] Input data file ID */
  char *fl_out; /* [sng] Output data file */
  int out_id; /* [id] Output data file ID */

  NcapVector<dmn_sct*> *ptr_dmn_in_vtr;  //Vector of dimensions in input file nb doesn't change
  NcapVector<dmn_sct*> *ptr_dmn_out_vtr; //Vector of dimensions in output file file
  NcapVector<sym_sct*> *ptr_sym_vtr;     //Vector of functions nb doesn't change
  NcapVarVector *ptr_var_vtr;            // list of attributes & variables
  bool ntl_scn; /* [flg] Initial scan of script */
} prs_sct;


/* Begin funtions in ncap_utl.c */

var_sct *                  /* O [sct] initialized variable */
ncap_var_init(
const char * const var_nm, /* I [sng] variable name constant */
prs_sct *prs_arg);          /* I/O  vectors of atts,vars,dims, filenames */


int                /* O  [bool] bool - ture if sucessful */
ncap_var_write     /*   [fnc] Write var to output file prs_arg->fl_out */ 
(var_sct *var,     /* I  [sct] variable to be written - freed at end */  
prs_sct *prs_arg); /* I/O vectors of atts & vars & file names  */


var_sct *                /* O [sct] variable containing attribute */
ncap_att_init(           /*   [fnc] Grab an attribute from input file */
const char *const va_nm, /* I [sng] att name of form var_nm&att_nm */ 
prs_sct *prs_arg);        /* I/O vectors of atts & vars & file names  */

sym_sct *                    /* O [sct] return sym_sct */
ncap_sym_init                /*  [fnc] populate & return a symbol table structure */
(const char * const sym_nm,  /* I [sng] symbol name */
 double (*fnc_dbl)(double),  /* I [fnc_dbl] Pointer to double function */
 float (*fnc_flt)(float));    /* I [fnc_flt] Pointer to float  function */



var_sct *   /* O [sct] Remainder of modulo operation of input variables (var_1%var_2) */
ncap_var_var_mod /* [fnc] Remainder (modulo) operation of two variables */
(var_sct *var_1, /* I [sc,t] Variable structure containing field */
 var_sct *var_2); /* I [sct] Variable structure containing divisor */


var_sct *         /* O [sct] Empowerment of input variables (var_1^var_2) */
ncap_var_var_pwr  /* [fnc] Empowerment of two variables */ 
(var_sct *var_1,  /* I [sct] Variable structure containing base */
 var_sct *var_2); /* I [sct] Variable structure containing exponent */


var_sct *           /* O [sct] Resultant variable (actually is var_in) */
ncap_var_fnc(       /* Apply function to var */   
var_sct *var_in,    /* I/O [sng] input variable */ 
sym_sct *app);       /* I [fnc_ptr] to apply to variable */


var_sct *         /* O [sct] Resultant variable (actually is var) */
ncap_var_abs(     /* Purpose: Find absolute value of each element of var */
var_sct *var);    /* I/O [sct] input variable */


nm_id_sct *            /* O [sct] new copy of xtr_lst */
nco_var_lst_copy(      /*   [fnc] Purpose: Copy xtr_lst and return new list */
nm_id_sct *xtr_lst,    /* I  [sct] input list */ 
int lst_nbr);           /* I  [nbr] number of elements in list */


nm_id_sct *             /* O [sct] New list */
nco_var_lst_sub(        /* [fnc] subract elements of lst_b from )lst */
nm_id_sct *xtr_lst,     /* I [sct] input list */   
int *nbr_xtr,           /* I/O [ptr] size of xtr_lst and new list */
nm_id_sct *xtr_lst_b,   /* I [sct] list to be subtracted */   
int nbr_lst_b);          /* I [nbr] size eof xtr_lst_b */ 


nm_id_sct *            /* O [sct] -- new list */
nco_var_lst_add(       /* [fnc]  add elemenst of lst_a to lst */
nm_id_sct *xtr_lst,    /* I [sct] input list */ 
int *nbr_xtr,          /* I/O [ptr] -- size of xtr_lst & new output list */ 
nm_id_sct *xtr_lst_a,  /* I [sct] list of elemenst to be added to new list */
int nbr_lst_a);         /* I [nbr] size of xtr_lst_a */


nm_id_sct *               /* O [sct] List of dimensions associated with input variable list */
nco_dmn_lst               /* [fnc] Create list of all dimensions in file  */
(const int nc_id,         /* I [id] netCDF input-file ID */
 int * const nbr_dmn);    /* O [nbr] Number of dimensions in  list */


nm_id_sct *                /* O [sct] output list */ 
nco_att_lst_mk      
(const int in_id,         /* I [id] of input file */
 const int out_id,        /* I [id] id of output file */
 NcapVarVector &var_vtr,  /* I [vec] vector of vars & att */
 int *nbr_lst);            /* O [ptr] size of output list */



nco_bool /* O [flg] Variables now conform */
ncap_var_stretch /* [fnc] Stretch variables */
(var_sct **var_1, /* I/O [ptr] First variable */
 var_sct **var_2); /* I/O [ptr] Second variable */

nco_bool          /* O [flg] true if var has been stretched */
ncap_att_stretch  /* stretch a single valued attribute from 1 to sz */
(var_sct* var,    /* I/O [sct] variable */       
 long nw_sz);      /* I [nbr] new var size */

var_sct *         /* O [sct] Sum of input variables (var1+var2) */
ncap_var_var_op   /* [fnc] Add two variables */
(var_sct *var1,  /* I [sct] Input variable structure containing first operand */
 var_sct *var2,  /* I [sct] Input variable structure containing second operand */
 int op);        /* Operation +-% */

  
bool            /* O [flg] true if all var elemenst are true */
ncap_var_lgcl   /* [fnc] calculate a aggregate bool value from a variable */
(var_sct* var);  /* I [sct] input variable */


var_sct*         /* O [sct] casting variable has its own private dims */ 
ncap_mk_cst(   /* [fnc] create casting var from a list of dims */
char **sbs_lst,  /* I [sng] Array of dimension subscripts */
int lst_nbr,     /* I [nbr] size of above list */  
prs_sct *prs_arg);



var_sct*
ncap_do_cst(
var_sct* var,
var_sct* var_cst,
bool bntlscn);



/* End funtions in ncap_utl.c */


#endif /* NCAP_H */




