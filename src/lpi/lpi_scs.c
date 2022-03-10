/**
 *  __   __  ___   _______  _______  _______  ______    _______
 * |  | |  ||   | |       ||       ||   _   ||    _ |  |       |
 * |  |_|  ||   | |  _____||_     _||  |_|  ||   | ||  |_     _|
 * |       ||   | | |_____   |   |  |       ||   |_||_   |   |
 * |       ||   | |_____  |  |   |  |       ||    __  |  |   |
 *  |     | |   |  _____| |  |   |  |   _   ||   |  | |  |   |
 *   |___|  |___| |_______|  |___|  |__| |__||___|  |_|  |___|
 */

 /**
  * @file   lpi_scs.c
  * @ingroup LPIS
  * @brief  SCIP LP interface for SCS.
  * @author vistart<i@vistart.me>
  * @copyright Copyright (c) 2021 - 2022 vistart
  * @link https://vistart.me/
  * @link https://github.com/vistart/scip
  * @license https://vistart.me/license/
  */

  /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "lpi/lpi.h"
#include "scip/pub_message.h"
#include "scs.h"

#define LPINAME            "SCS"                             /**< name of the LPI interface */
#define LPIINFINITY        1e+20                             /**< infinity value */
#define ABS(x) ((x)>0?(x):-(x))                              /**< get absolute value of x */
#define LPIINFINITESIMAL   1e-10                             /**< infinitily small value */
#define ISLPIINFINITESIMAL(x) (ABS(x)<LPIINFINITESIMAL)      /**< determine whether the x is infinitesimal */

/* globally turn off lint warnings: */
/*lint --e{715}*/

struct SCIP_Column
{
    SCIP_Real obj;
    SCIP_Real lb;
    SCIP_Real ub;
    char* name;
};

struct SCIP_Columns
{
    struct SCIP_Column** columns_ptr;
    int                  ncols;
};

struct SCIP_Row
{
    SCIP_Real  lhs;
    SCIP_Real  rhs;
    char* name;
    SCIP_Real* objs;
};

struct SCIP_Rows
{
    struct SCIP_Row** rows_ptr;
    int                nrows;
};

/** LP interface */
struct SCIP_LPi
{
    //SPxSCIP*              scs;                /**< our SCS implementation */
    ScsData* scsdata;
    ScsCone* scscone;
    ScsSettings* scsstgs;
    ScsSolution* scssol;
    ScsInfo* scsinfo;
    ScsWork* scswork;            /**< our SCS work structure *//* 暂时不启用 */
    int* cstat;              /**< array for storing column basis status *//* 暂时不启用 */
    int* rstat;              /**< array for storing row basis status *//* 暂时不启用 */
    int                   cstatsize;          /**< size of cstat array *//* 暂时不启用 */
    int                   rstatsize;          /**< size of rstat array *//* 暂时不启用 */
    // SCIP_PRICING          pricing;            /**< current pricing strategy */
    // SCIP_Bool             solved;             /**< was the current LP solved? */
    //SLUFactor*            factorization;      /**< factorization possibly needed for basis inverse */
    // SCIP_Real             rowrepswitch;       /**< use row representation if number of rows divided by number of columns exceeds this value */
    // SCIP_Real             conditionlimit;     /**< maximum condition number of LP basis counted as stable (-1.0: no limit) */
    // SCIP_Bool             checkcondition;     /**< should condition number of LP basis be checked for stability? */
    SCIP_MESSAGEHDLR* messagehdlr;        /**< messagehdlr handler to printing messages, or NULL *//* 暂时不启用 */
    //int                   nrows;              /**< number of rows */
    //int                   ncols;              /**< number of columns */
    SCIP_OBJSEN           objsen;             /**< objective sense */
    const char* name;               /**< problem name */
    struct SCIP_Columns* columns;
    struct SCIP_Rows* rows;
};

/** LPi state stores basis information */
struct SCIP_LPiState
{
    int                   ncols;              /**< number of LP columns */
    int                   nrows;              /**< number of LP rows */
    //COLPACKET*            packcstat;          /**< column basis status in compressed form */
    //ROWPACKET*            packrstat;          /**< row basis status in compressed form */
};

/** LPi norms to store dual steepest edge */
struct SCIP_LPiNorms
{
    int                   nrows;              /**< number of stored norms corresponding to rows */
    int                   ncols;              /**< number of stored norms corresponding to cols */
    SCIP_Real* norms;              /**< norms to be (re)stored */
};

/*
 * Local Methods
 */
 /** error handling method */
static
void errorMessageAbort(
    void
)
{  /*lint --e{2707}*/
    SCIPerrorMessage("SCS is not ready to use (LPS=scs).\n");
    SCIPerrorMessage("Ensure <lp/solvefreq = -1>; note that continuous variables might require an LP-solver.\n");
    SCIPABORT();
}

/** error handling method */
static
void errorMessage(
    void
)
{
    SCIPerrorMessage("SCS is not ready to use (LPS=scs).\n");
    SCIPerrorMessage("Ensure <lp/solvefreq = -1>; note that continuous variables might require an LP-solver.\n");
}

/*
 * LP Interface Methods
 */

 /*
  * Miscellaneous Methods
  */
  /**
   * 获取指定列的下界，即某个变量的下界。
   * @param lpi 指向线性求解器接口结构体的指针。
   * @param col 列号。从 0 开始。
   * @return 指定列的下界值。
   * 注意！下界有可能为负无穷大。取得值后需与 -LPIINFINITY 比较。
   */
const SCIP_Real get_column_lower_bound_real(
    SCIP_LPI* lpi,
    int col
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    return lpi->columns->columns_ptr[col]->lb;
}

/**
 * 设置指定列的下界，即某个变量的下界。
 * 注意！如果设置值小于负无穷大（-LPIINFINITY），则只记录负无穷大。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @param val 指定列新下界值。
 * @return 设置成功。
 */
SCIP_RETCODE set_column_lower_bound_real(
    SCIP_LPI* lpi,
    int col,
    SCIP_Real val
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    if (val < -LPIINFINITY)
    {
        val = -LPIINFINITY;
    }
    lpi->columns->columns_ptr[col]->lb = val;
    return SCIP_OKAY;
}

/**
 * 获取指定列的上界，即某个变量的上界。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @return 指定列的上界值。
 * 注意！上界有可能为正无穷大。取得值后需与 LPIINFINITY 比较。
 */
SCIP_Real get_column_upper_bound_real(
    SCIP_LPI* lpi,
    int col
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    return lpi->columns->columns_ptr[col]->ub;
}

/**
 * 设置指定列的上界，即某个变量的上界。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @param val 指定列新上界值。
 * @return 设置成功。
 */
SCIP_RETCODE set_column_upper_bound_real(
    SCIP_LPI* lpi,
    int col,
    SCIP_Real val
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    if (val > LPIINFINITY)
    {
        val = LPIINFINITY;
    }
    lpi->columns->columns_ptr[col]->ub = val;
    return SCIP_OKAY;
}

/**
 * 获得指定列的目标值系数，即某个变量在目标函数中的系数。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @return 指定列的目标值系数。
 */
SCIP_Real get_column_obj_real(
    SCIP_LPI* lpi,
    int col
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    return lpi->columns->columns_ptr[col]->obj;
}

/**
 * 设置指定列的目标值系数，即某个变量在目标函数中的系数。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @param val 指定列的目标值系数。
 * @return 设置成功。
 */
SCIP_RETCODE set_column_obj_real(
    SCIP_LPI* lpi,
    int col,
    SCIP_Real val
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    lpi->columns->columns_ptr[col]->obj = val;
    return SCIP_OKAY;
}

/**
 * 获得指定列的名称，即变量名。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @return 指定列的名称。
 */
const char* get_column_name(
    SCIP_LPI* lpi,
    int col
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    return lpi->columns->columns_ptr[col]->name;
}

/**
 * 设置指定列的名称，即变量名。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @param val 新的列名。
 * @return 设置成功。
 */
SCIP_RETCODE set_column_name(
    SCIP_LPI* lpi,
    int col,
    char* val
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    lpi->columns->columns_ptr[col]->name = val;
    return SCIP_OKAY;
}

/**
 * 获取列总数。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @return 列总数。
 */
int get_ncols(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    return lpi->columns->ncols;
}

/**
 * 调试打印列信息。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @return 执行成功。
 */
SCIP_RETCODE debug_print_column(
    SCIP_LPI* lpi,
    int col
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->columns_ptr != NULL);
    SCIPdebugMessage("Col[%d]: %20s, obj: %8.2f, (%8.2f, %8.2f)\n", col, get_column_name(lpi, col),
        get_column_obj_real(lpi, col), get_column_lower_bound_real(lpi, col), get_column_upper_bound_real(lpi, col));
    return SCIP_OKAY;
}

/**
 * 调试打印所有列信息。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @return 执行成功。
 */
SCIP_RETCODE debug_print_all_columns(
    SCIP_LPI* lpi
)
{
    SCIPdebugMessage("calling debugPrintAllColumns.\n");
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    for (int i = 0; i < lpi->columns->ncols; i++)
    {
        debug_print_column(lpi, i);
    }
    return SCIP_OKAY;
}

/**
 * 释放指定列。
 * 注意：释放操作后指针指向 NULL，并不会缩小columns_ptr的尺寸，也不会改变其它任何columns_ptr指针的值，仅为方便清空列（clear_columns）而设。
 * 因此，如果要删除某个列，则应当在删除后主动挪动后续指针值，并修改列总数（ncols）。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @return 释放成功。
 */
SCIP_RETCODE free_column(
    SCIP_LPI* lpi,
    int col
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(col < lpi->columns->ncols);
    if (lpi->columns->columns_ptr[col] != NULL)
    {
        free(lpi->columns->columns_ptr[col]);
    }
    return SCIP_OKAY;
}

/**
 * 重新设置列集合尺寸。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param newsize 新尺寸。
 * @return 设置成功。
 */
SCIP_RETCODE resize_columns(
    SCIP_LPI* lpi,
    int newsize
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    if (lpi->columns->columns_ptr == NULL)
    {
        lpi->columns->columns_ptr = (struct SCIP_Column**)malloc(sizeof(struct SCIP_Column*) * newsize);
    }
    else if (newsize > 0)
    {
        lpi->columns->columns_ptr = realloc(lpi->columns->columns_ptr, sizeof(struct SCIP_Column*) * newsize);
    }
    lpi->columns->ncols = newsize;
    return SCIP_OKAY;
}

/**
 * 初始化列。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param col 列号。从 0 开始。
 * @return 初始化成功。
 */
SCIP_RETCODE init_column(
    SCIP_LPI* lpi,
    int col
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->columns_ptr != NULL);
    assert(lpi->columns->ncols > col);
    /*
    if (lpi->columns->columns_ptr[col] != NULL)
    {
        free(lpi->columns->columns_ptr[col]);
    }*/
    lpi->columns->columns_ptr[col] = (struct SCIP_Column*)malloc(sizeof(struct SCIP_Column));
    return SCIP_OKAY;
}

/**
 * 初始化所有列。
 * 注意！初始化成功后不能直接设置或获取任何列信息，而需要先重新调整列尺寸。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @return 初始化成功。
 */
SCIP_RETCODE init_columns(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    /*
    if (lpi->columns != NULL)
    {
        free(lpi->columns);
    }*/
    lpi->columns = (struct SCIP_Columns*)malloc(sizeof(struct SCIP_Columns));
    lpi->columns->columns_ptr = NULL;
    lpi->columns->ncols = 0;
    return SCIP_OKAY;
}

/**
 * 清空所有列。
 * 注意！您需要自行保证所有列信息结构体地址真实有效。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @return 清空成功。
 */
SCIP_RETCODE clear_columns(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    if (lpi->columns == NULL)
    {
        return SCIP_OKAY;
    }
    for (int i = 0; i < lpi->columns->ncols; i++)
    {
        free_column(lpi, i);
    }
    init_columns(lpi);
    return SCIP_OKAY;
}

SCIP_Real get_row_lhs_real(
    SCIP_LPI* lpi,
    int row
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->nrows > 0);
    assert(lpi->rows->nrows > row);
    assert(lpi->rows->rows_ptr != NULL);
    assert(lpi->rows->rows_ptr[row] != NULL);
    return lpi->rows->rows_ptr[row]->lhs;
}

SCIP_RETCODE set_row_lhs_real(
    SCIP_LPI* lpi,
    int row,
    SCIP_Real val
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->nrows > 0);
    assert(lpi->rows->nrows > row);
    assert(lpi->rows->rows_ptr != NULL);
    assert(lpi->rows->rows_ptr[row] != NULL);
    if (val < -LPIINFINITY)
    {
        val = -LPIINFINITY;
    }
    lpi->rows->rows_ptr[row]->lhs = val;
    return SCIP_OKAY;
}

SCIP_Real get_row_rhs_real(
    SCIP_LPI* lpi,
    int row
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->nrows > 0);
    assert(lpi->rows->nrows > row);
    assert(lpi->rows->rows_ptr != NULL);
    assert(lpi->rows->rows_ptr[row] != NULL);
    return lpi->rows->rows_ptr[row]->rhs;
}

SCIP_RETCODE set_row_rhs_real(
    SCIP_LPI* lpi,
    int row,
    SCIP_Real val
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->nrows > 0);
    assert(lpi->rows->nrows > row);
    assert(lpi->rows->rows_ptr != NULL);
    assert(lpi->rows->rows_ptr[row] != NULL);
    if (val > LPIINFINITY)
    {
        val = LPIINFINITY;
    }
    lpi->rows->rows_ptr[row]->rhs = val;
    return SCIP_OKAY;
}

const char* get_row_name(
    SCIP_LPI* lpi,
    int row
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->nrows > 0);
    assert(lpi->rows->nrows > row);
    assert(lpi->rows->rows_ptr != NULL);
    assert(lpi->rows->rows_ptr[row] != NULL);
    return lpi->rows->rows_ptr[row]->name;
}

SCIP_RETCODE set_row_name(
    SCIP_LPI* lpi,
    int row,
    char* name
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->nrows > 0);
    assert(lpi->rows->nrows > row);
    assert(lpi->rows->rows_ptr != NULL);
    assert(lpi->rows->rows_ptr[row] != NULL);
    lpi->rows->rows_ptr[row]->name = name;
    return SCIP_OKAY;
}

SCIP_Real get_row_obj_real(
    SCIP_LPI* lpi,
    int row,
    int col
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->nrows > row);
    assert(lpi->rows->rows_ptr != NULL);
    return lpi->rows->rows_ptr[row]->objs[col];
}

SCIP_RETCODE set_row_obj_real(
    SCIP_LPI* lpi,
    int row,
    int col,
    SCIP_Real val
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->nrows > row);
    assert(lpi->rows->rows_ptr != NULL);
    lpi->rows->rows_ptr[row]->objs[col] = val;
    return SCIP_OKAY;
}

/**
 * 获取行总数。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @return 行总数。
 */
int get_nrows(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    return lpi->rows->nrows;
}

/**
 * 释放指定行。
 * 注意：释放操作后指针指向 NULL，并不会缩小rows_ptr的尺寸，也不会改变其它任何rows_ptr指针的值，仅为方便清空列（clear_rows）而设。
 * 因此，如果要删除某个行，则应当在删除后主动挪动后续指针值，并修改行总数（nrows）。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param row 行号。从 0 开始。
 * @return 释放成功。
 */
SCIP_RETCODE free_row(
    SCIP_LPI* lpi,
    int row
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(row < lpi->rows->nrows);
    if (lpi->rows->rows_ptr[row] != NULL)
    {
        free(lpi->rows->rows_ptr[row]);
    }
    return SCIP_OKAY;
}

SCIP_RETCODE resize_rows(
    SCIP_LPI* lpi,
    int newsize
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    if (lpi->rows->rows_ptr == NULL)
    {
        lpi->rows->rows_ptr = (struct SCIP_Row**)malloc(sizeof(struct SCIP_Row*) * newsize);
    }
    else if (newsize > 0)
    {
        lpi->rows->rows_ptr = realloc(lpi->rows->rows_ptr, sizeof(struct SCIP_Row*) * newsize);
    }
    lpi->rows->nrows = newsize;
    return SCIP_OKAY;
}

SCIP_RETCODE free_rowobj(
    SCIP_LPI* lpi,
    int row
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->rows_ptr != NULL);
    assert(lpi->rows->nrows > row);
    if (lpi->rows->rows_ptr[row]->objs != NULL)
    {
        free(lpi->rows->rows_ptr[row]->objs);
    }
    lpi->rows->rows_ptr[row]->objs = (SCIP_Real*)malloc(0);
    return SCIP_OKAY;
}

SCIP_RETCODE resize_row_objs(
    SCIP_LPI* lpi,
    int row,
    int newsize
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->rows_ptr != NULL);
    assert(lpi->rows->nrows > row);
    lpi->rows->rows_ptr[row]->objs = realloc(lpi->rows->rows_ptr[row]->objs, sizeof(SCIP_Real) * newsize);
    return SCIP_OKAY;
}

SCIP_RETCODE init_row(
    SCIP_LPI* lpi,
    int row
)
{
    assert(lpi != NULL);
    assert(lpi->rows != NULL);
    assert(lpi->rows->rows_ptr != NULL);
    assert(lpi->rows->nrows > row);
    lpi->rows->rows_ptr[row] = (struct SCIP_Row*)malloc(sizeof(struct SCIP_Row));
    lpi->rows->rows_ptr[row]->objs = (SCIP_Real*)malloc(0);
    resize_row_objs(lpi, row, get_ncols(lpi));
    memset(lpi->rows->rows_ptr[row]->objs, 0, sizeof(SCIP_Real) * get_ncols(lpi));
    return SCIP_OKAY;
}

SCIP_RETCODE init_rows(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    lpi->rows = (struct SCIP_Rows*)malloc(sizeof(struct SCIP_Rows));
    lpi->rows->rows_ptr = NULL;
    lpi->rows->nrows = 0;
    return SCIP_OKAY;
}

SCIP_RETCODE clear_rows(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    if (lpi->rows == NULL)
    {
        return SCIP_OKAY;
    }
    for (int i = 0; i < lpi->rows->nrows; i++)
    {
        free_row(lpi, i);
    }
    init_rows(lpi);
    return SCIP_OKAY;
}

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
static char scsname[100];
const char* SCIPlpiGetSolverName(
    void
)
{
    snprintf(scsname, 100, "%s %s", LPINAME, scs_version()); /*lint !e778*/
    return scsname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
    void
)
{
    return "Linear Programming Solver using Splitting Conic Solver Developed By Zhao Vistart.";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
    SCIP_LPI* lpi                 /**< 指向线性求解器接口结构体的指针 */
)
{  /*lint --e{715}*/
    return (void*)lpi->scswork;
}

/** pass integrality information to LP solver */
SCIP_RETCODE SCIPlpiSetIntegralityInformation(
    SCIP_LPI* lpi,                /**< 指向线性求解器接口结构体的指针 */
    int                   ncols,              /**< length of integrality array */
    int* intInfo             /**< integrality array (0: continuous, 1: integer).
                                              May be NULL iff ncols is 0.  */
)
{ /*lint --e{715}*/
    SCIPerrorMessage("SCIPlpiSetIntegralityInformation() has not been implemented yet.\n");
    return SCIP_LPERROR;
}

/** informs about availability of a primal simplex solving method */
SCIP_Bool SCIPlpiHasPrimalSolve(
    void
)
{
    return TRUE;
}

/** informs about availability of a dual simplex solving method */
SCIP_Bool SCIPlpiHasDualSolve(
    void
)
{
    return TRUE;
}

/** informs about availability of a barrier solving method */
SCIP_Bool SCIPlpiHasBarrierSolve(
    void
)
{
    return FALSE;
}

/**@} */

ScsMatrix* ConstructPMatrix(int n)
{
    double* P_x = (double*)calloc(1, 0);
    int* P_i = (int*)calloc(1, 0);
    int* P_p = (int*)calloc(1, sizeof(int));
    ScsMatrix* m = (ScsMatrix*)calloc(1, sizeof(ScsMatrix));
    m->x = P_x;
    m->i = P_i;
    m->p = P_p;
    m->m = n;
    m->n = n;
    return m;
}


/*
 * LPI Creation and Destruction Methods
 */

 /**@name LPI Creation and Destruction Methods */
 /**@{ */

 /** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
    SCIP_LPI** lpi,                /**< 指向线性求解器接口结构体的指针 */
    SCIP_MESSAGEHDLR* messagehdlr,        /**< message handler to use for printing messages, or NULL */
    const char* name,               /**< problem name */
    SCIP_OBJSEN           objsen              /**< objective sense */
)
{  /*lint --e{715}*/
    //assert(sizeof(SCIP_REAL) == sizeof(scs_float)); /** 检查 SCIP 实数是否与 scs 实数类型相一致，一致才能计算。 */
    SCIPdebugMessage("SCIPlpiCreate()\n");

    assert(lpi != NULL);
    assert(name != NULL);
    SCIP_ALLOC(BMSallocMemory(lpi));
    SCIPdebugMessage("Name: %s\n", name);
    (*lpi)->name = name;
    SCIPdebugMessage("ObjSen: %d\n", objsen);
    (*lpi)->objsen = objsen;
    SCIPdebugMessage("Note that the SCIP is creating an SCS work...\n");
    (*lpi)->scscone = (ScsCone*)calloc(1, sizeof(ScsCone));
    (*lpi)->scsdata = (ScsData*)calloc(1, sizeof(ScsData));
    (*lpi)->scsstgs = (ScsSettings*)calloc(1, sizeof(ScsSettings));
    (*lpi)->scssol = (ScsSolution*)calloc(1, sizeof(ScsSolution));
    (*lpi)->scsinfo = (ScsInfo*)calloc(1, sizeof(ScsInfo));

    /* Utility to set default settings */
    /** TODO: 修改参数 */
    scs_set_default_settings((*lpi)->scsstgs);

    /* Modify tolerances */
    /** TODO: 修改参数 */
    (*lpi)->scsstgs->eps_abs = 1e-9;
    (*lpi)->scsstgs->eps_rel = 1e-9;

    SCIPdebugMessage("size of scs_int = %lu, size of scs_float = %lu\n", sizeof(scs_int), sizeof(scs_float));
    // 初始化列在前，初始化行在后。
    init_columns(*lpi);
    init_rows(*lpi);
    return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
    SCIP_LPI** lpi                 /**< 指向线性求解器接口结构体的指针 */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    SCIPdebugMessage("SCIPlpiFree()\n");

    /* Free allocated memory */
    free((*lpi)->scscone);
    free((*lpi)->scsdata);
    free((*lpi)->scsstgs);
    free((*lpi)->scsinfo);
    /* SCS allocates sol->x,y,s if NULL on entry, need to be freed */
    free((*lpi)->scssol->x);
    free((*lpi)->scssol->y);
    free((*lpi)->scssol->s);
    free((*lpi)->scssol);

    BMSfreeMemory(lpi);

    return SCIP_OKAY;
}

/**@} */

/*
 * Modification Methods
 */

 /**@name Modification Methods */
 /**@{ */

 /** copies LP data with column matrix into LP solver */
SCIP_RETCODE SCIPlpiLoadColLP(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_OBJSEN           objsen,             /**< objective sense */
    int                   ncols,              /**< number of columns */
    const SCIP_Real* obj,                /**< objective function values of columns */
    const SCIP_Real* lb,                 /**< lower bounds of columns */
    const SCIP_Real* ub,                 /**< upper bounds of columns */
    char** colnames,           /**< column names, or NULL */
    int                   nrows,              /**< number of rows */
    const SCIP_Real* lhs,                /**< left hand sides of rows */
    const SCIP_Real* rhs,                /**< right hand sides of rows */
    char** rownames,           /**< row names, or NULL */
    int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
    const int* beg,                /**< start index of each column in ind- and val-array */
    const int* ind,                /**< row indices of constraint matrix entries */
    const SCIP_Real* val                 /**< values of constraint matrix entries */
)
{  /*lint --e{715}*/

#ifndef NDEBUG
    {
        SCIPdebugMessage("SCIPlpiLoadColLP:\n");
        for (int j = 0; j < nnonz; j++) {
            assert(!ISLPIINFINITESIMAL(val[j]));
            printf("val[%d]: %f\n", j, val[j]);
        }
    }
#endif

    assert(lpi != NULL);
    assert(lhs != NULL);
    assert(rhs != NULL);
    assert(obj != NULL);
    assert(lb != NULL);
    assert(ub != NULL);
    assert(beg != NULL);
    assert(ind != NULL);
    assert(val != NULL);
    assert(get_ncols(lpi) >= 0);
    assert(get_nrows(lpi) >= 0);

    return SCIP_OKAY;
}
/**
SCIP_RETCODE addConstraintByCol(
        SCIP_LPI*             lpi,
        int                   index,
        int                   ncols,
        char*                 colname,
        const SCIP_Real       obj,
        const SCIP_Real       lb,
        const SCIP_Real       ub
)
{
    assert(lpi != NULL);
    assert(lpi->constraints != NULL);
    assert(lpi->bounds != NULL);
    int n = 0;
    // 添加小于等于条件
    if (!SCIPlpiIsInfinity(lpi, -lb))
    {
        // 添加 b 向量
        lpi->bounds = (SCIP_Real*)realloc(lpi->bounds, sizeof(SCIP_Real) * (++lpi->nconstraints));
        lpi->bounds[lpi->nconstraints - 1] = -lb;
        // 添加 A 矩阵
        lpi->constraints = (SCIP_Real**)realloc(lpi->constraints, sizeof(SCIP_Real*) * lpi->nconstraints); // 新增一行，用于保存新的变量边界。
        lpi->constraints[lpi->nconstraints - 1] = (SCIP_Real*)malloc(sizeof(SCIP_Real) * lpi->n);
        memset(lpi->constraints[lpi->nconstraints - 1], 0, lpi->n);
        lpi->constraints[lpi->nconstraints - 1][index] = -obj;
    }
    // 添加大于等于条件
    if (!SCIPlpiIsInfinity(lpi, ub))
    {
        // 添加 b 向量
        lpi->bounds = (SCIP_Real*)realloc(lpi->bounds, sizeof(SCIP_Real) * (++lpi->nconstraints));
        lpi->bounds[lpi->nconstraints - 1] = ub;
        // 添加 A 矩阵
        lpi->constraints = (SCIP_Real**)realloc(lpi->constraints, sizeof(SCIP_Real*) * lpi->nconstraints); // 新增一行，用于保存新的变量边界。
        lpi->constraints[lpi->nconstraints - 1] = (SCIP_Real*)malloc(sizeof(SCIP_Real) * lpi->n);
        memset(lpi->constraints[lpi->nconstraints - 1], 0, lpi->n);
        lpi->constraints[lpi->nconstraints - 1][index] = obj;
    }

    if (!SCIPlpiIsInfinity(lpi, -ub))
    { // 添加小于等于条件
        int old_size = lpi->nconstraints;
        resizeConstraints(lpi, lpi->nconstraints + ncols);
        for (int i = 0; i < ncols; i++)
        {
            lpi->constraints[old_size + i] = (i == index) ? ub : 0;
        }
    }
    if (!SCIPlpiIsInfinity(lpi, lb))
    {
        int old_size = lpi->nconstraints;
        resizeConstraints(lpi, lpi->nconstraints + ncols);
        for (int i = 0; i < ncols; i++)
        {
            lpi->constraints[old_size + i] = (i == index) ? -lb : 0;
        }
    }
    return SCIP_OKAY;
}*/

/** adds columns to the LP *//* variables */
SCIP_RETCODE SCIPlpiAddCols(
    SCIP_LPI* lpi,      /**< LP interface structure */
    int                   ncols,    /**< number of columns to be added */
    const SCIP_Real* obj,      /**< objective function values of new columns */
    const SCIP_Real* lb,       /**< lower bounds of new columns */
    const SCIP_Real* ub,       /**< upper bounds of new columns */
    char** colnames, /**< column names, or NULL */
    int                   nnonz,    /**< number of nonzero elements to be added to the constraint matrix */
    const int* beg,      /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
    const int* ind,      /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
    const SCIP_Real* val       /**< values of constraint matrix entries, or NULL if nnonz == 0 */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(get_ncols(lpi) >= 0);
    assert(obj != NULL);
    assert(lb != NULL);
    assert(ub != NULL);
    assert(nnonz == 0 || beg != NULL);
    assert(nnonz == 0 || ind != NULL);
    assert(nnonz == 0 || val != NULL);
    assert(nnonz >= 0);
    assert(ncols >= 0);

#ifndef NDEBUG
    if (nnonz > 0) {
        const int nrows = get_nrows(lpi);
        for (int j = 0; j < nnonz; ++j)
        {
            assert(0 <= ind[j] && ind[j] < nrows);
            assert(!ISLPIINFINITESIMAL(val[j]));
            SCIPdebugMessage("beg[i], ind[i], val[i]: (%d, %d, %f)\n", beg[j], ind[j], val[j]);
        }
    }
#endif
    /**
     * TODO:
     * 1. 获取现在的c向量和n值；
     * 2. 如果某个向量索引已存在，则替换；不存在，则新增。
     * 3. 第1条和第2条的单元测试。
     */
     /** int start, last; */
     /**
     for (int i = 0; i < ncols; ++i) {
         if (nnonz > 0) {
             start = beg[i];
             last = (i == ncols - 1 ? nnonz : beg[i + 1]);
         }
         if (ncols > lpi->n) {
             // 当前索引超出向量尺寸，应当扩容。
             lpi->n = ncols;
             lpi->scsdata->c = realloc(lpi->scsdata->c, sizeof(double) * lpi->n);
         }
         lpi->scsdata->c[i] = obj[i];
     }
     lpi->scsdata->P = ConstructPMatrix(lpi->n);
     for (int i = 0; i < ncols; ++i)
     {
         SCIPdebugMessage("c[%d]: %s, obj: %f, [%f, %f]\n", i, colnames[i], lpi->scsdata->c[i], lb[i], ub[i]);
         addConstraintByCol(lpi, i, ncols, colnames[i], obj[i], lb[i], ub[i]);
     }
     /**
      * TODO:
      * 1. 接收变量边界。
      * 2. 第1条的单元测试。
      *
 #ifdef SCIP_DEBUG
     SCIPdebugMessage("Columns added:\n");
     for (int i = 0; i < lpi->nconstraints; i++)
     {
         SCIPdebugMessage("Cons[%d]:\n", i);
         for (int j = 0; j < lpi->n; j++)
         {
             SCIPdebugMessage(" %6.2f", lpi->constraints[i][j]);
         }
         SCIPdebugMessage(" %6.2f", lpi->bounds[i]);
         SCIPdebugMessage("\n");
     }
 #endif*/
    int oldncols = get_ncols(lpi);
    resize_columns(lpi, oldncols + ncols);
    for (int i = 0; i < ncols; i++)
    {
        init_column(lpi, oldncols + i);
        set_column_obj_real(lpi, oldncols + i, obj[i]);
        set_column_lower_bound_real(lpi, oldncols + i, lb[i]);
        set_column_upper_bound_real(lpi, oldncols + i, ub[i]);
        set_column_name(lpi, oldncols + i, colnames[i]);
    }
#ifdef SCIP_DEBUG
    debug_print_all_columns(lpi);
#endif
    return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   firstcol,           /**< first column to be deleted */
    int                   lastcol             /**< last column to be deleted */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(get_ncols(lpi) >= 0);

    //lpi->ncols -= lastcol - firstcol + 1;
    assert(get_ncols(lpi) >= 0);

    return SCIP_OKAY;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position
 *
 * @TODO
 *
 */
SCIP_RETCODE SCIPlpiDelColset(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int* dstat               /**< deletion status of columns
                                          *   input:  1 if column should be deleted, 0 if not
                                          *   output: new position of column, -1 if column was deleted */
)
{  /*lint --e{715}*/
    int cnt = 0;

    assert(lpi != NULL);
    assert(dstat != NULL);
    assert(get_ncols(lpi) >= 0);

    for (int j = 0; j < get_ncols(lpi); ++j)
    {
        if (dstat[j])
        {
            ++cnt;
            dstat[j] = -1;
        }
        else
            dstat[j] = cnt;
    }
    //lpi->ncols -= cnt;
    assert(get_ncols(lpi) >= 0);

    return SCIP_OKAY;
}


/**
 * 添加约束条件，按SCIP提供的行。
 *
 */
 /**
SCIP_RETCODE addConstraintsByRow(
    SCIP_LPI* lpi,
    const int             row_index,
    const SCIP_Real lhs,
    const SCIP_Real rhs,
    const int* beg,
    const int* ind,
    const double* val
)
{
    SCIPdebugMessage("calling addConstraintsByRow()\n");
    assert(lpi != NULL);
    assert(lpi->constraints != NULL);
    assert(lpi->bounds != NULL);
    if (!SCIPlpiIsInfinity(lpi, -lhs))
    {
        lpi->constraints = (SCIP_Real**)realloc(lpi->constraints, sizeof(SCIP_Real*) * (++lpi->nconstraints)); // 新增一行，用于保存新的变量边界。
        lpi->constraints[lpi->nconstraints - 1] = (SCIP_Real*)malloc(sizeof(SCIP_Real) * lpi->n);
        lpi->bounds = (SCIP_Real*)realloc(lpi->bounds, sizeof(SCIP_Real) * lpi->nconstraints);
        lpi->bounds[lpi->nconstraints - 1] = -lhs;
        for (int j = beg[row_index]; j < beg[row_index] + lpi->n; ++j)
        {
            lpi->constraints[lpi->nconstraints - 1][j - beg[row_index]] = -val[j];
        }
    }
    if (!SCIPlpiIsInfinity(lpi, rhs))
    {
        lpi->constraints = (SCIP_Real**)realloc(lpi->constraints, sizeof(SCIP_Real*) * (++lpi->nconstraints)); // 新增一行，用于保存新的变量边界。
        lpi->constraints[lpi->nconstraints - 1] = (SCIP_Real*)malloc(sizeof(SCIP_Real) * lpi->n);
        lpi->bounds = (SCIP_Real*)realloc(lpi->bounds, sizeof(SCIP_Real) * lpi->nconstraints);
        lpi->bounds[lpi->nconstraints - 1] = rhs;
        for (int j = beg[row_index]; j < beg[row_index] + lpi->n; ++j)
        {
            lpi->constraints[lpi->nconstraints - 1][j - beg[row_index]] = val[j];
        }
    }
#ifdef SCIP_DEBUG
    for (int i = 0; i < lpi->n; ++i)
    {
        //SCIPdebugMessage("lpi->constraints[%d]: %f\n", row_index * lpi->n + i, lpi->constraints[row_index * lpi->n + i]);
    }
#endif
    return SCIP_OKAY;
}*/

/** adds rows to the LP *//* constraints */
SCIP_RETCODE SCIPlpiAddRows(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   nrows,              /**< number of rows to be added */
    const SCIP_Real* lhs,                /**< left hand sides of new rows */
    const SCIP_Real* rhs,                /**< right hand sides of new rows */
    char** rownames,           /**< row names, or NULL */
    int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
    const int* beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
    const int* ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
    const SCIP_Real* val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
)
{  /*lint --e{715}*/
    SCIPdebugMessage("calling SCIPlpiAddRows()\n");
    assert(lpi != NULL);
    assert(get_nrows(lpi) >= 0);
    assert(lhs != NULL);
    assert(rhs != NULL);
    assert(nnonz == 0 || beg != NULL);
    assert(nnonz == 0 || ind != NULL);
    assert(nnonz == 0 || val != NULL);

#ifndef NDEBUG
    /* perform check that no new columns are added - this is forbidden */
    {
        for (int j = 0; j < nnonz; ++j)
        {
            assert(!ISLPIINFINITESIMAL(val[j]));
            assert(0 <= ind[j] && ind[j] < get_ncols(lpi));
        }
    }
#endif
    /*
    lpi->m = nrows;
    if (!lpi->constraints)
    {
        lpi->constraints = (double*)calloc(1, sizeof(double) * lpi->n * lpi->m);
    }
    for (int i = 0; i < nrows; ++i)
    {
        SCIPdebugMessage("r[%d]: %s,", i, rownames[i]);
        SCIPdebugMessage("%f <= ", lhs[i]);
        for (int j = beg[i]; j < beg[i] + lpi->n; ++j)
        {
            SCIPdebugMessage(" %f ", val[j]);
            //lpi->constraints[j] = val[j];
        }
        SCIPdebugMessage(" <= %f\n", rhs[i]);
        addConstraintsByRow(lpi, i, lhs[i], rhs[i], beg, ind, val);
    }
#ifdef SCIP_DEBUG
    SCIPdebugMessage("Rows added:\n");
    for (int i = 0; i < lpi->nconstraints; i++)
    {
        SCIPdebugMessage("Cons[%d]:\n", i);
        for (int j = 0; j < lpi->n; j++)
        {
            SCIPdebugMessage(" %6.2f", lpi->constraints[i][j]);
        }
        SCIPdebugMessage(" %6.2f", lpi->bounds[i]);
        SCIPdebugMessage("\n");
    }
#endif*/
#pragma region 添加行定义，不含涉及变量系数
    int oldnrows = get_nrows(lpi);
    resize_rows(lpi, oldnrows + nrows);
    for (int i = 0; i < nrows; i++)
    {
        init_row(lpi, oldnrows + i);
        set_row_name(lpi, oldnrows + i, rownames[i]);
        set_row_lhs_real(lpi, oldnrows + i, lhs[i]);
        set_row_rhs_real(lpi, oldnrows + i, rhs[i]);
    }
#pragma endregion
#pragma region 为新添加的行更新变量系数
    for (int i = 0; i < nrows; i++)
    {
        for (int j = beg[i]; (i + 1 < nrows) ? (j < beg[i + 1]) : (j < nnonz); j++)
        {
            set_row_obj_real(lpi, i, ind[j], val[j]);
        }
    }
#pragma endregion
    return SCIP_OKAY;
}

/** deletes all rows in the given range from LP
 *
 * @TODO
 */
SCIP_RETCODE SCIPlpiDelRows(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   firstrow,           /**< first row to be deleted */
    int                   lastrow             /**< last row to be deleted */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(get_nrows(lpi) >= 0);

    //lpi->nrows -= lastrow - firstrow + 1;
    assert(get_nrows(lpi) >= 0);

    return SCIP_OKAY;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position
 *
 * @TODO
 */
SCIP_RETCODE SCIPlpiDelRowset(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int* dstat               /**< deletion status of rows
                                          *   input:  1 if row should be deleted, 0 if not
                                          *   output: new position of row, -1 if row was deleted */
)
{  /*lint --e{715}*/
    int cnt = 0;

    assert(lpi != NULL);
    assert(dstat != NULL);
    assert(get_nrows(lpi) >= 0);

    for (int i = 0; i < get_nrows(lpi); ++i)
    {
        if (dstat[i])
        {
            ++cnt;
            dstat[i] = -1;
        }
        else
            dstat[i] = cnt;
    }
    //lpi->nrows -= cnt;
    assert(get_nrows(lpi) >= 0);

    return SCIP_OKAY;
}

/** clears the whole LP
 *
 * @TODO
 */
SCIP_RETCODE SCIPlpiClear(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(get_nrows(lpi) >= 0);
    assert(get_ncols(lpi) >= 0);
    // 先清理行，后清理列。
    clear_rows(lpi);
    clear_columns(lpi);
    return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   ncols,              /**< number of columns to change bounds for */
    const int* ind,                /**< column indices or NULL if ncols is zero */
    const SCIP_Real* lb,                 /**< values for the new lower bounds or NULL if ncols is zero */
    const SCIP_Real* ub                  /**< values for the new upper bounds or NULL if ncols is zero */
)
{  /*lint --e{715}*/

    assert(ncols == 0 || (ind != NULL && lb != NULL && ub != NULL));

    if (ncols <= 0) {
        return SCIP_OKAY;
    }
    for (int j = 0; j < ncols; ++j)
    {
        if (SCIPlpiIsInfinity(lpi, lb[j]))
        {
            SCIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", ind[j]);
            return SCIP_LPERROR;
        }
        if (SCIPlpiIsInfinity(lpi, -ub[j]))
        {
            SCIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", ind[j]);
            return SCIP_LPERROR;
        }
    }

    return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiChgSides(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   nrows,              /**< number of rows to change sides for */
    const int* ind,                /**< row indices */
    const SCIP_Real* lhs,                /**< new values for left hand sides */
    const SCIP_Real* rhs                 /**< new values for right hand sides */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(ind != NULL);
    assert(lhs != NULL);
    assert(rhs != NULL);
    return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiChgCoef(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   row,                /**< row number of coefficient to change */
    int                   col,                /**< column number of coefficient to change */
    SCIP_Real             newval              /**< new value of coefficient */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_OBJSEN           objsen              /**< new objective sense */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return SCIP_OKAY;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   ncols,              /**< number of columns to change objective value for */
    const int* ind,                /**< column indices to change objective value for */
    const SCIP_Real* obj                 /**< new objective values for columns */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(ind != NULL);
    assert(obj != NULL);
    return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   row,                /**< row number to scale */
    SCIP_Real             scaleval            /**< scaling multiplier */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return SCIP_OKAY;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiScaleCol(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   col,                /**< column number to scale */
    SCIP_Real             scaleval            /**< scaling multiplier */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return SCIP_OKAY;
}

/**@} */




/*
 * Data Accessing Methods
 */

 /**@name Data Accessing Methods */
 /**@{ */

 /** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int* nrows               /**< pointer to store the number of rows */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(nrows != NULL);
    assert(get_nrows(lpi) >= 0);

    *nrows = get_nrows(lpi);

    return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int* ncols               /**< pointer to store the number of cols */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(ncols != NULL);
    assert(get_ncols(lpi) >= 0);

    *ncols = get_ncols(lpi);

    return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int* nnonz               /**< pointer to store the number of nonzeros */
)
{  /*lint --e{715}*/
    assert(nnonz != NULL);
    assert(lpi != NULL);
    assert(lpi->scsdata != NULL);
    assert(lpi->scsdata->A != NULL);
    assert(lpi->scsdata->A->p != NULL);
    double* p = lpi->scsdata->A->p;
    *nnonz = 0;
    while (p != NULL) {
        *nnonz++;
        p++;
    }
    return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetCols(
    SCIP_LPI* lpi,       /**< LP interface structure */
    int                   firstcol,  /**< first column to get from LP */
    int                   lastcol,   /**< last column to get from LP */
    SCIP_Real* lb,        /**< buffer to store the lower bound vector, or NULL */
    SCIP_Real* ub,        /**< buffer to store the upper bound vector, or NULL */
    int* nnonz,     /**< pointer to store the number of nonzero elements returned, or NULL */
    int* beg,    /**< buffer to store start index of each column in ind- and val-array, or NULL */
    int* ind,    /**< buffer to store row indices of constraint matrix entries, or NULL */
    SCIP_Real* val     /**< buffer to store values of constraint matrix entries, or NULL */
)
{  /*lint --e{715}*/
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetRows(
    SCIP_LPI* lpi,      /**< LP interface structure */
    int                   firstrow, /**< first row to get from LP */
    int                   lastrow,  /**< last row to get from LP */
    SCIP_Real* lhs,      /**< buffer to store left hand side vector, or NULL */
    SCIP_Real* rhs,      /**< buffer to store right hand side vector, or NULL */
    int* nnonz,    /**< pointer to store the number of nonzero elements returned, or NULL */
    int* beg,      /**< buffer to store start index of each row in ind- and val-array, or NULL */
    int* ind,      /**< buffer to store column indices of constraint matrix entries, or NULL */
    SCIP_Real* val       /**< buffer to store values of constraint matrix entries, or NULL */
)
{  /*lint --e{715}*/
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets column names */
SCIP_RETCODE SCIPlpiGetColNames(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   firstcol,           /**< first column to get name from LP */
    int                   lastcol,            /**< last column to get name from LP */
    char** colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1)
                                              or NULL if namestoragesize is zero */
    char* namestorage,        /**< storage for col names or NULL if namestoragesize is zero */
    int                   namestoragesize,    /**< size of namestorage
                                              (if 0, storageleft returns the storage needed) */
    int* storageleft      /**< amount of storage left (if < 0 the namestorage was not big enough)
                                           or NULL if namestoragesize is zero */
)
{ /*lint --e{715}*/
    assert(lpi != NULL);
    assert(colnames != NULL || namestoragesize == 0);
    assert(namestorage != NULL || namestoragesize == 0);
    assert(namestoragesize >= 0);
    assert(storageleft != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets row names */
SCIP_RETCODE SCIPlpiGetRowNames(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   firstrow,           /**< first row to get name from LP */
    int                   lastrow,            /**< last row to get name from LP */
    char** rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1)
                                              or NULL if namestoragesize is zero */
    char* namestorage,        /**< storage for row names or NULL if namestoragesize is zero */
    int                   namestoragesize,    /**< size of namestorage
                                              (if 0, -storageleft returns the storage needed) */
    int* storageleft      /**< amount of storage left (if < 0 the namestorage was not big enough)
                                           or NULL if namestoragesize is zero */
)
{ /*lint --e{715}*/
    assert(lpi != NULL);
    assert(rownames != NULL || namestoragesize == 0);
    assert(namestorage != NULL || namestoragesize == 0);
    assert(namestoragesize >= 0);
    assert(storageleft != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_OBJSEN* objsen              /**< pointer to store objective sense */
)
{  /*lint --e{715}*/
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/**
 * 从线性规划问题对象中获得一批目标函数系数。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param firstcol 起始列号。
 * @param lastcol 结束列号。
 * @param vals 指定列对应的系数。
 * @return 获取成功。
 */
SCIP_RETCODE SCIPlpiGetObj(
    SCIP_LPI*  lpi,
    int        firstcol,
    int        lastcol,
    SCIP_Real* vals
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(0 <= firstcol && firstcol <= lastcol && lastcol < get_ncols(lpi));
    assert(vals != NULL);
    for (int i = firstcol; i <= lastcol; i++)
    {
        vals[i - firstcol] = get_column_obj_real(lpi, i);
    }
    return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   firstcol,           /**< first column to get bounds for */
    int                   lastcol,            /**< last column to get bounds for */
    SCIP_Real* lbs,                /**< array to store lower bound values, or NULL */
    SCIP_Real* ubs                 /**< array to store upper bound values, or NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(firstcol <= lastcol);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   firstrow,           /**< first row to get sides for */
    int                   lastrow,            /**< last row to get sides for */
    SCIP_Real* lhss,               /**< array to store left hand side values, or NULL */
    SCIP_Real* rhss                /**< array to store right hand side values, or NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(firstrow <= lastrow);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   row,                /**< row number of coefficient */
    int                   col,                /**< column number of coefficient */
    SCIP_Real* val                 /**< pointer to store the value of the coefficient */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(val != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/**@} */




/*
 * Solving Methods
 */

 /**@name Solving Methods */
 /**@{ */

 /** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_Bool             crossover           /**< perform crossover */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   col,                /**< column to apply strong branching on */
    SCIP_Real             psol,               /**< fractional current primal solution value of column */
    int                   itlim,              /**< iteration limit for strong branchings */
    SCIP_Real* down,               /**< stores dual bound after branching column down */
    SCIP_Real* up,                 /**< stores dual bound after branching column up */
    SCIP_Bool* downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                          *   otherwise, it can only be used as an estimate value */
    SCIP_Bool* upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                          *   otherwise, it can only be used as an estimate value */
    int* iter                /**< stores total number of strong branching iterations, or -1;
                                              may be NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(down != NULL);
    assert(up != NULL);
    assert(downvalid != NULL);
    assert(upvalid != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPlpiStrongbranchesFrac(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int* cols,               /**< columns to apply strong branching on */
    int                   ncols,              /**< number of columns */
    SCIP_Real* psols,              /**< fractional current primal solution values of columns */
    int                   itlim,              /**< iteration limit for strong branchings */
    SCIP_Real* down,               /**< stores dual bounds after branching columns down */
    SCIP_Real* up,                 /**< stores dual bounds after branching columns up */
    SCIP_Bool* downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                          *   otherwise, they can only be used as an estimate values */
    SCIP_Bool* upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                          *   otherwise, they can only be used as an estimate values */
    int* iter                /**< stores total number of strong branching iterations, or -1;
                                              may be NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(cols != NULL);
    assert(psols != NULL);
    assert(down != NULL);
    assert(up != NULL);
    assert(downvalid != NULL);
    assert(upvalid != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPlpiStrongbranchInt(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   col,                /**< column to apply strong branching on */
    SCIP_Real             psol,               /**< current integral primal solution value of column */
    int                   itlim,              /**< iteration limit for strong branchings */
    SCIP_Real* down,               /**< stores dual bound after branching column down */
    SCIP_Real* up,                 /**< stores dual bound after branching column up */
    SCIP_Bool* downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                          *   otherwise, it can only be used as an estimate value */
    SCIP_Bool* upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                          *   otherwise, it can only be used as an estimate value */
    int* iter                /**< stores total number of strong branching iterations, or -1;
                                              may be NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(down != NULL);
    assert(up != NULL);
    assert(downvalid != NULL);
    assert(upvalid != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPlpiStrongbranchesInt(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int* cols,               /**< columns to apply strong branching on */
    int                   ncols,              /**< number of columns */
    SCIP_Real* psols,              /**< current integral primal solution values of columns */
    int                   itlim,              /**< iteration limit for strong branchings */
    SCIP_Real* down,               /**< stores dual bounds after branching columns down */
    SCIP_Real* up,                 /**< stores dual bounds after branching columns up */
    SCIP_Bool* downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                          *   otherwise, they can only be used as an estimate values */
    SCIP_Bool* upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                          *   otherwise, they can only be used as an estimate values */
    int* iter                /**< stores total number of strong branching iterations, or -1;
                                              may be NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(cols != NULL);
    assert(psols != NULL);
    assert(down != NULL);
    assert(up != NULL);
    assert(downvalid != NULL);
    assert(upvalid != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}
/**@} */




/*
 * Solution Information Methods
 */

 /**@name Solution Information Methods */
 /**@{ */

 /** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** gets information about primal and dual feasibility of the current LP solution
 *
 *  The feasibility information is with respect to the last solving call and it is only relevant if SCIPlpiWasSolved()
 *  returns true. If the LP is changed, this information might be invalidated.
 *
 *  Note that @a primalfeasible and @a dualfeasible should only return true if the solver has proved the respective LP
 *  to be feasible. Thus, the return values should be equal to the values of SCIPlpiIsPrimalFeasible() and
 *  SCIPlpiIsDualFeasible(), respectively. Note that if feasibility cannot be proved, they should return false (even if
 *  the problem might actually be feasible).
 */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_Bool* primalfeasible,     /**< pointer to store primal feasibility status */
    SCIP_Bool* dualfeasible        /**< pointer to store dual feasibility status */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(primalfeasible != NULL);
    assert(dualfeasible != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff current LP solution is stable
 *
 *  This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
 *  infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
 *  version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
 *  SCIPlpiIsStable() should return false.
 */
SCIP_Bool SCIPlpiIsStable(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_Bool* success             /**< pointer to store, whether the instability could be ignored */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(success != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_Real* objval              /**< stores the objective value */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(objval != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets primal and dual solution vectors for feasible LPs
 *
 *  Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
 *  SCIPlpiIsOptimal() returns true.
 */
SCIP_RETCODE SCIPlpiGetSol(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_Real* objval,             /**< stores the objective value, may be NULL if not needed */
    SCIP_Real* primsol,            /**< primal solution vector, may be NULL if not needed */
    SCIP_Real* dualsol,            /**< dual solution vector, may be NULL if not needed */
    SCIP_Real* activity,           /**< row activity vector, may be NULL if not needed */
    SCIP_Real* redcost             /**< reduced cost vector, may be NULL if not needed */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_Real* ray                 /**< primal ray */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(ray != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_Real* dualfarkas          /**< dual Farkas row multipliers */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(dualfarkas != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
    SCIP_LPI* lpi,             /**< LP interface structure */
    int* iterations       /**< pointer to store the number of iterations of the last solve call */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(iterations != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for @p quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
    SCIP_Real* quality             /**< pointer to store quality number */
)
{ /*lint --e{715}*/
    assert(lpi != NULL);
    assert(quality != NULL);

    *quality = SCIP_INVALID;

    return SCIP_OKAY;
}

/**@} */




/*
 * LP Basis Methods
 */

 /**@name LP Basis Methods */
 /**@{ */

 /** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int* cstat,              /**< array to store column basis status, or NULL */
    int* rstat               /**< array to store row basis status, or NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
    SCIP_LPI* lpi,                /**< LP interface structure */
    const int* cstat,              /**< array with column basis status */
    const int* rstat               /**< array with row basis status */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(cstat != NULL);
    assert(rstat != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiGetBasisInd(
    SCIP_LPI* lpi,          /**< LP interface structure */
    int* bind          /**< pointer to store basis indices ready to keep number of rows entries */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(bind != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** get row of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvRow(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   r,                  /**< row number */
    SCIP_Real* coef,               /**< pointer to store the coefficients of the row */
    int* inds,               /**< array to store the non-zero indices, or NULL */
    int* ninds               /**< pointer to store the number of non-zero indices, or NULL
                                          *   (-1: if we do not store sparsity information) */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(coef != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvCol(
    SCIP_LPI* lpi,          /**< LP interface structure */
    int                   c,            /**< column number of B^-1; this is NOT the number of the column in the LP;
                                          *  you have to call SCIPlpiGetBasisInd() to get the array which links the
                                          *  B^-1 column numbers to the row and column numbers of the LP!
                                          *  c must be between 0 and nrows-1, since the basis has the size
                                          *  nrows * nrows */
    SCIP_Real* coef,         /**< pointer to store the coefficients of the column */
    int* inds,         /**< array to store the non-zero indices, or NULL */
    int* ninds         /**< pointer to store the number of non-zero indices, or NULL
                                          *   (-1: if we do not store sparsity information) */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(coef != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** get row of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP
 *  solver uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be
 *  negated; see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvARow(
    SCIP_LPI* lpi,            /**< LP interface structure */
    int                   r,              /**< row number */
    const SCIP_Real* binvrow,        /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
    SCIP_Real* coef,           /**< vector to return coefficients of the row */
    int* inds,           /**< array to store the non-zero indices, or NULL */
    int* ninds           /**< pointer to store the number of non-zero indices, or NULL
                                          *   (-1: if we do not store sparsity information) */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(coef != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** get column of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP
 *  solver uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be
 *  negated; see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvACol(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   c,                  /**< column number */
    SCIP_Real* coef,               /**< vector to return coefficients of the column */
    int* inds,               /**< array to store the non-zero indices, or NULL */
    int* ninds               /**< pointer to store the number of non-zero indices, or NULL
                                          *   (-1: if we do not store sparsity information) */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(coef != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/**@} */




/*
 * LP State Methods
 */

 /**@name LP State Methods */
 /**@{ */

 /** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
    SCIP_LPI* lpi,                /**< LP interface structure */
    BMS_BLKMEM* blkmem,             /**< block memory */
    SCIP_LPISTATE** lpistate            /**< pointer to LPi state information (like basis information) */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(blkmem != NULL);
    assert(lpistate != NULL);
    assert(blkmem != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
    SCIP_LPI* lpi,                /**< LP interface structure */
    BMS_BLKMEM* blkmem,             /**< block memory */
    const SCIP_LPISTATE* lpistate            /**< LPi state information (like basis information), or NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(blkmem != NULL);
    assert(lpistate != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
    SCIP_LPI* lpi,                /**< LP interface structure */
    BMS_BLKMEM* blkmem,             /**< block memory */
    SCIP_LPISTATE** lpistate            /**< pointer to LPi state information (like basis information) */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpistate != NULL);
    assert(blkmem != NULL);
    return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_LPISTATE* lpistate            /**< LP state information (like basis information), or NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessageAbort();
    return FALSE;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
    SCIP_LPI* lpi,                /**< LP interface structure */
    const char* fname               /**< file name */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(fname != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
    SCIP_LPI* lpi,                /**< LP interface structure */
    const char* fname               /**< file name */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(fname != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/**@} */




/*
 * LP Pricing Norms Methods
 */

 /**@name LP Pricing Norms Methods */
 /**@{ */

 /** stores LPi pricing norms information
  *  @todo should we store norm information?
  */
SCIP_RETCODE SCIPlpiGetNorms(
    SCIP_LPI* lpi,                /**< LP interface structure */
    BMS_BLKMEM* blkmem,             /**< block memory */
    SCIP_LPINORMS** lpinorms            /**< pointer to LPi pricing norms information */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(blkmem != NULL);
    assert(lpinorms != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** loads LPi pricing norms into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetNorms()
 */
SCIP_RETCODE SCIPlpiSetNorms(
    SCIP_LPI* lpi,                /**< LP interface structure */
    BMS_BLKMEM* blkmem,             /**< block memory */
    const SCIP_LPINORMS* lpinorms            /**< LPi pricing norms information, or NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
    SCIP_LPI* lpi,                /**< LP interface structure */
    BMS_BLKMEM* blkmem,             /**< block memory */
    SCIP_LPINORMS** lpinorms            /**< pointer to LPi pricing norms information, or NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/**@} */




/*
 * Parameter Methods
 */

 /**@name Parameter Methods */
 /**@{ */

 /** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_LPPARAM          type,               /**< parameter number */
    int* ival                /**< buffer to store the parameter value */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(ival != NULL);
    return SCIP_PARAMETERUNKNOWN;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_LPPARAM          type,               /**< parameter number */
    int                   ival                /**< parameter value */
)
{  /*lint --e{715}*/
    SCIPdebugMessage("calling SCIPlpiSetIntpar()\n");
    assert(lpi != NULL);

    switch (type)
    {
    case SCIP_LPPAR_FROMSCRATCH:
        assert(ival == TRUE || ival == FALSE);
        break;
    case SCIP_LPPAR_FASTMIP:
        assert(ival == TRUE || ival == FALSE);
        break;
    case SCIP_LPPAR_REFACTOR:
        assert(ival == TRUE || ival == FALSE);
        break;
    case SCIP_LPPAR_LPINFO:
        assert(ival == TRUE || ival == FALSE);
        break;
    case SCIP_LPPAR_LPITLIM:
        assert(ival >= 0);
        /* -1 <= ival, -1 meaning no time limit, 0 stopping immediately */
        if (ival >= INT_MAX) {
            ival = -1;
        }
        break;
    case SCIP_LPPAR_PRESOLVING:
        assert(ival == TRUE || ival == FALSE);
        break;
    case SCIP_LPPAR_PRICING:
        break;
    case SCIP_LPPAR_SCALING:
        assert(ival == TRUE || ival == FALSE);
        break;
    case SCIP_LPPAR_TIMING:
        assert(ival >= 0 && ival < 3);
        break;
    case SCIP_LPPAR_RANDOMSEED:
        break;
    case SCIP_LPPAR_POLISHING:
        assert(ival >= 0 && ival < 3);
        break;
    default:
        return SCIP_PARAMETERUNKNOWN;
    }  /*lint !e788*/

    return SCIP_OKAY;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_LPPARAM          type,               /**< parameter number */
    SCIP_Real* dval                /**< buffer to store the parameter value */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(dval != NULL);
    return SCIP_PARAMETERUNKNOWN;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_LPPARAM          type,               /**< parameter number */
    SCIP_Real             dval                /**< parameter value */
)
{  /*lint --e{715}*/
    SCIPdebugMessage("calling SCIPlpiSetRealpar()\n");
    assert(lpi != NULL);
    switch (type)
    {
    case SCIP_LPPAR_FEASTOL:
        break;
    case SCIP_LPPAR_DUALFEASTOL:
        break;
    case SCIP_LPPAR_OBJLIM:
        break;
    case SCIP_LPPAR_LPTILIM:
        break;
    case SCIP_LPPAR_ROWREPSWITCH:
        break;
    case SCIP_LPPAR_CONDITIONLIMIT:
        break;
    case SCIP_LPPAR_BARRIERCONVTOL:
        break;
    default:
        return SCIP_PARAMETERUNKNOWN;
    }  /*lint !e788*/
    return SCIP_OKAY;
}

/** interrupts the currently ongoing lp solve or disables the interrupt */
SCIP_RETCODE SCIPlpiInterrupt(
    SCIP_LPI* lpi,            /**< LP interface structure */
    SCIP_Bool             interrupt       /**< TRUE if interrupt should be set, FALSE if it should be disabled */
)
{
    /*lint --e{715}*/
    assert(lpi != NULL);

    return SCIP_OKAY;
}

/**@} */

/*
 * Numerical Methods
 */

 /**@name Numerical Methods */
 /**@{ */

 /** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return LPIINFINITY;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_Real             val                 /**< value to be checked for infinity */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    if (val >= LPIINFINITY) {
        return TRUE;
    }
    return FALSE;
}

/**@} */

/*
 * File Interface Methods
 */

 /**@name File Interface Methods */
 /**@{ */

 /** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
    SCIP_LPI* lpi,                /**< LP interface structure */
    const char* fname               /**< file name */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(fname != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
    SCIP_LPI* lpi,                /**< LP interface structure */
    const char* fname               /**< file name */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    errorMessage();
    return SCIP_PLUGINNOTFOUND;
}

/**@} */
