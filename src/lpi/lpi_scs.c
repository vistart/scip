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
#include "scip/bitencode.h"
#include "scip/pub_message.h"
#include "scs.h"

#define LPINAME            "SCS"                             /**< name of the LPI interface */
#define LPIINFINITY        1e+20                             /**< infinity value */
#define ABS(x) ((x)>0?(x):-(x))                              /**< get absolute value of x */
#define LPIINFINITESIMAL   1e-10                             /**< infinitily small value */
#define ISLPIINFINITESIMAL(x) (ABS(x)<LPIINFINITESIMAL)      /**< determine whether the x is infinitesimal */

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

/* globally turn off lint warnings: */
/*lint --e{715}*/

struct SCIP_Column
{
    SCIP_Real obj;
    SCIP_Real lb;
    SCIP_Real ub;
    char*     name;
    int       intInfo;
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
    ScsData*              scsdata;
    ScsCone*              scscone;
    ScsSettings*          scsstgs;
    ScsSolution*          scssol;
    ScsInfo*              scsinfo;
    ScsWork*              scswork;            /**< our SCS work structure */
    int*                  cstat;              /**< array for storing column basis status */
    int*                  rstat;              /**< array for storing row basis status */
    int                   cstatsize;          /**< size of cstat array */
    int                   rstatsize;          /**< size of rstat array */
    // SCIP_PRICING          pricing;            /**< current pricing strategy */
    // SCIP_Bool             solved;             /**< was the current LP solved? */
    //SLUFactor*            factorization;      /**< factorization possibly needed for basis inverse */
    // SCIP_Real             rowrepswitch;       /**< use row representation if number of rows divided by number of columns exceeds this value */
    // SCIP_Real             conditionlimit;     /**< maximum condition number of LP basis counted as stable (-1.0: no limit) */
    // SCIP_Bool             checkcondition;     /**< should condition number of LP basis be checked for stability? */
    SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL *//* 暂时不启用 */
    //int                   nrows;              /**< number of rows */
    //int                   ncols;              /**< number of columns */
    SCIP_OBJSEN           objsen;             /**< objective sense */
    const char* name;               /**< problem name */
    struct SCIP_Columns*  columns;
    struct SCIP_Rows*     rows;
    int                   nconsbycol;
    SCIP_Real             objlim;
    SCIP_Real             feastol;
    SCIP_Real             dualfeastol;
    SCIP_Real             lptilim;
    SCIP_Real             rowrepswitch;
    SCIP_Real             conditionlimit;
    SCIP_Bool             checkcondition;
    SCIP_Real             markowitz;
    SCIP_Bool             fromscratch;
    SCIP_Bool             lpinfo;
    int                   lpitlim;
    SCIP_Longint          presoving;
    SCIP_PRICING          pricing;
    SCIP_Longint          pricer;
    SCIP_Longint          scaling;
    SCIP_Longint          timing;
    SCIP_Longint          randomseed;
    SCIP_Longint          polishing;
    SCIP_Longint          refactor;
};

/** LPi state stores basis information */
struct SCIP_LPiState
{
    int                   ncols;              /**< number of LP columns */
    int                   nrows;              /**< number of LP rows */
    COLPACKET*            packcstat;          /**< column basis status in compressed form */
    ROWPACKET*            packrstat;          /**< row basis status in compressed form */
};

/** LPi norms to store dual steepest edge */
struct SCIP_LPiNorms
{
    int                   nrows;              /**< number of stored norms corresponding to rows */
    int                   ncols;              /**< number of stored norms corresponding to cols */
    SCIP_Real* norms;                         /**< norms to be (re)stored */
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

int get_column_integrality(
    SCIP_LPI* lpi,
    int col
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    return lpi->columns->columns_ptr[col]->intInfo;
}

SCIP_RETCODE set_column_integrality(
    SCIP_LPI* lpi,
    int col,
    int intInfo
)
{
    assert(lpi != NULL);
    assert(lpi->columns != NULL);
    assert(lpi->columns->ncols > 0);
    assert(lpi->columns->ncols > col);
    assert(lpi->columns->columns_ptr != NULL);
    lpi->columns->columns_ptr[col]->intInfo = intInfo;
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
        lpi->columns->columns_ptr = (struct SCIP_Column**)calloc(sizeof(struct SCIP_Column*), newsize);
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
    lpi->columns->columns_ptr[col] = (struct SCIP_Column*)calloc(sizeof(struct SCIP_Column), 1);
    lpi->columns->columns_ptr[col]->intInfo = 0;
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
    lpi->columns = (struct SCIP_Columns*)calloc(sizeof(struct SCIP_Columns), 1);
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

/**
 * 获取行左边界值。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param row 行号。
 * @return 获取成功。
 */
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

/**
 * 设置行左边界值。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param row 行号。
 * @param val 左边界值。
 * @return 设置成功。
 */
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
        lpi->rows->rows_ptr = (struct SCIP_Row**)calloc(sizeof(struct SCIP_Row*), newsize);
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
    lpi->rows->rows_ptr[row]->objs = (SCIP_Real*)calloc(sizeof(SCIP_Real), 0);
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
    lpi->rows->rows_ptr[row] = (struct SCIP_Row*)calloc(sizeof(struct SCIP_Row), 1);
    lpi->rows->rows_ptr[row]->objs = (SCIP_Real*)calloc(sizeof(SCIP_Real), 0);
    resize_row_objs(lpi, row, get_ncols(lpi));
    memset(lpi->rows->rows_ptr[row]->objs, 0, sizeof(SCIP_Real) * get_ncols(lpi));
    return SCIP_OKAY;
}

SCIP_RETCODE init_rows(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    lpi->rows = (struct SCIP_Rows*)calloc(sizeof(struct SCIP_Rows), 1);
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

SCIP_RETCODE init_state(
    SCIP_LPI* lpi,
    int nrows,
    int ncols
)
{
    assert(lpi != NULL);
    lpi->cstat = (int*)calloc(get_ncols(lpi), sizeof(int));
    lpi->rstat = (int*)calloc(get_nrows(lpi), sizeof(int));
    return SCIP_OKAY;
}

SCIP_RETCODE resize_state_rows(
    SCIP_LPI* lpi,
    int nrows
)
{
    assert(lpi != NULL);
    lpi->rstat = realloc(lpi->rstat, nrows * sizeof(int));
    return SCIP_OKAY;
}

SCIP_RETCODE resize_state_columns(
    SCIP_LPI* lpi,
    int ncols
)
{
    assert(lpi != NULL);
    lpi->cstat = realloc(lpi->cstat, ncols * sizeof(int));
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
    //SCIPdebugMessage("SCS doesn't support the integrality of variables.\n");
    SCIPdebugMessage("calling SCIPlpiSetIntegralityInformation().\n");
    for (int i = 0; i < ncols; i++)
    {
        set_column_integrality(lpi, i, intInfo[i]);
    }
    return SCIP_OKAY;
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

/*
 * LPi state methods
 */
/**
 * 返回需要存储列组合信息的组合数。
 * @param ncols 列数
 * @return 组合数。
 */
static
int colpacketNum(
    int ncols
)
{
    return (ncols + (int)COLS_PER_PACKET - 1) / (int)COLS_PER_PACKET;
}

/**
 * 返回需要存储行组合信息的组合数。
 * @param nrows 行数
 * @return 组合数。
 */
static
int rowpacketNum(
    int nrows
)
{
    return (nrows + (int)ROWS_PER_PACKET - 1) / (int)ROWS_PER_PACKET;
}

/**
 *
 */
static
void lpistatePack(
    SCIP_LPISTATE* lpistate,
    const int* cstat,
    const int* rstat
)
{
    assert(lpistate != NULL);
    assert(lpistate->packcstat != NULL);
    assert(lpistate->packrstat != NULL);

    SCIPencodeDualBit(cstat, lpistate->packcstat, lpistate->ncols);
    SCIPencodeDualBit(rstat, lpistate->packrstat, lpistate->nrows);
}

/**
 *
 */
static
void lpistateUnpack(
    const SCIP_LPISTATE* lpistate,            /**< pointer to LPi state data */
    int* cstat,               /**< buffer for storing basis status of columns in unpacked format */
    int* rstat                /**< buffer for storing basis status of rows in unpacked format */
)
{
    assert(lpistate != NULL);
    assert(lpistate->packcstat != NULL);
    assert(lpistate->packrstat != NULL);

    SCIPdecodeDualBit(lpistate->packcstat, cstat, lpistate->ncols);
    SCIPdecodeDualBit(lpistate->packrstat, rstat, lpistate->nrows);
}

/**
 *
 */
static
SCIP_RETCODE lpistateCreate(
    SCIP_LPISTATE** lpistate,           /**< pointer to LPi state */
    BMS_BLKMEM* blkmem,             /**< block memory */
    int                   ncols,              /**< number of columns to store */
    int                   nrows               /**< number of rows to store */
)
{
    assert(lpistate != NULL);
    assert(blkmem != NULL);
    assert(ncols >= 0);
    assert(nrows >= 0);

    int nColPackets = colpacketNum(ncols);
    int nRowPackets = rowpacketNum(nrows);

    SCIP_ALLOC(BMSallocBlockMemory(blkmem, lpistate));
    SCIP_ALLOC(BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packcstat, nColPackets));
    SCIP_ALLOC(BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packrstat, nRowPackets));

    return SCIP_OKAY;
}

/**
 *
 */
static
void lpistateFree(
    SCIP_LPISTATE** lpistate,           /**< pointer to LPi state information (like basis information) */
    BMS_BLKMEM* blkmem              /**< block memory */
)
{
    assert(blkmem != NULL);
    assert(lpistate != NULL);
    assert(*lpistate != NULL);

    int nColPackets = colpacketNum((*lpistate)->ncols);
    int nRowPackets = rowpacketNum((*lpistate)->nrows);

    BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packcstat, nColPackets);
    BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packrstat, nRowPackets);
    BMSfreeBlockMemory(blkmem, lpistate);
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
#ifdef SCIP_DEBUG
    (*lpi)->scsstgs->verbose = 1;
#else
    (*lpi)->scsstgs->verbose = 0;
#endif

    /* Modify tolerances */
    /** TODO: 修改参数 */
    (*lpi)->scsstgs->eps_abs = 1e-9;
    (*lpi)->scsstgs->eps_rel = 1e-9;

    SCIPdebugMessage("size of scs_int = %lu, size of scs_float = %lu\n", sizeof(scs_int), sizeof(scs_float));
    // 初始化列在前，初始化行在后。
    init_columns(*lpi);
    init_rows(*lpi);
    init_state(*lpi, get_nrows(*lpi), get_ncols(*lpi));
    (*lpi)->nconsbycol = 0;
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
    for (int i = lastcol; i >= firstcol; i--)
    {
        free_column(lpi, i);
        // 将后续指针向前挪一个位置。
        for (int j = i; j < get_ncols(lpi) - 1; j++)
        {
            lpi->columns->columns_ptr[j] = lpi->columns->columns_ptr[j + 1];
        }
        resize_columns(lpi, get_ncols(lpi) - 1);
    }

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
    for (int i = lastrow; i >= firstrow; i--)
    {
        free_row(lpi, i);
        for (int j = i; j < get_nrows(lpi) - 1; j++)
        {
            lpi->rows->rows_ptr[j] = lpi->rows->rows_ptr[j + 1];
        }
        resize_rows(lpi, get_nrows(lpi) - 1);
    }
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
    for (int i = 0; i < ncols; i++)
    {
        assert(0 <= ind[i] && ind[i] < get_ncols(lpi));
        if (SCIPlpiIsInfinity(lpi, lb[i]))
        {
            SCIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", ind[i]);
            return SCIP_LPERROR;
        }
        if (SCIPlpiIsInfinity(lpi, -ub[i]))
        {
            SCIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", ind[i]);
            return SCIP_LPERROR;
        }
#ifdef SCIP_DEBUG
        SCIP_Real old_lb = get_column_lower_bound_real(lpi, ind[i]);
        SCIP_Real old_ub = get_column_upper_bound_real(lpi, ind[i]);
#endif
        set_column_lower_bound_real(lpi, ind[i], lb[i]);
        set_column_upper_bound_real(lpi, ind[i], ub[i]);
#ifdef SCIP_DEBUG
        SCIPdebugMessage("the lower bound of column[%d]: %8.2f is replaced with %8.2f\n",
            ind[i], old_lb, get_column_lower_bound_real(lpi, ind[i]));
        SCIPdebugMessage("the upper bound of column[%d]: %8.2f is replaced with %8.2f\n",
            ind[i], old_ub, get_column_upper_bound_real(lpi, ind[i]));
#endif
        assert(get_column_lower_bound_real(lpi, ind[i]) <= get_column_upper_bound_real(lpi, ind[i]));
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
    *nnonz = (int)(lpi->scsdata->A->p[lpi->scsdata->A->n]);
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
 * 从线性规划问题对象中获得一批变量在目标函数中的系数。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param firstcol 起始列号。
 * @param lastcol 结束列号。
 * @param vals 指定列对应的系数。可以通过此参数传出获取结果。
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

/**
 * 从线性规划问题对象中获得一批变量的上下界。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param firstcol 起始列号。
 * @param lastcol 结束列号。
 * @param lbs 下界。可以通过此参数传出获取结果。
 * @param ubs 上界。可以通过此参数传出获取结果。
 * @return 获取成功。
 */
SCIP_RETCODE SCIPlpiGetBounds(
    SCIP_LPI*  lpi,
    int        firstcol,
    int        lastcol,
    SCIP_Real* lbs,
    SCIP_Real* ubs
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(0 <= firstcol && firstcol <= lastcol && lastcol < get_ncols(lpi));
    for (int i = firstcol; i <= lastcol; i++)
    {
        if (lbs != NULL)
        {
            lbs[i - firstcol] = get_column_lower_bound_real(lpi, i);
        }
        if (ubs != NULL)
        {
            ubs[i - firstcol] = get_column_upper_bound_real(lpi, i);
        }
    }
    return SCIP_OKAY;
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

/**
 * 根据 SCIP 提供的列（Column）构建 SCS 的 A 矩阵（约束）和对应的 c 向量（约束上界）。
 * @param lpi 指向线性求解器接口结构体的指针。
 * @param AMatrixOfColumns 仅包含列信息的A矩阵。
 * @param CVector 与A矩阵对应的 c 向量。
 * @paran nvector 行数。
 * @return 执行成功。
 */
SCIP_RETCODE ConstructAMatrixAndCVectorByColumns(
    SCIP_LPI*   lpi,
    scs_float*** AMatrixOfColumns,
    scs_float*** CVector,
    int*         nvector
)
{
    SCIPdebugMessage("calling ConstructAMatrixAndCVectorByColumns...\n");
    assert(lpi != NULL);
    const int ncols = get_ncols(lpi);
    if (lpi->columns == NULL || ncols == 0)
    {
        return SCIP_OKAY;
    }
    *nvector = 0;
    for (int i = 0; i < ncols; i++)
    {
        scs_float lb = get_column_lower_bound_real(lpi, i);
        if (!SCIPlpiIsInfinity(lpi, -lb))
        {
            ++*nvector;
            //*CVector = realloc(*CVector, sizeof(SCIP_Real) * ++*nvector);
            //(*CVector)[*nvector - 1] = -lb;
            //SCIPdebugMessage("(*CVector)[%d]: %8.2f\n", *nvector - 1, (*CVector)[*nvector - 1]);
            //*AMatrixOfColumns = realloc(*AMatrixOfColumns, sizeof(SCIP_Real*) * *nvector);
            //(*AMatrixOfColumns)[*nvector - 1] = (SCIP_Real*)calloc(sizeof(SCIP_Real), ncols);
            //(*AMatrixOfColumns)[*nvector - 1][i] = -get_column_obj_real(lpi, i);
        }
        scs_float ub = get_column_upper_bound_real(lpi, i);
        if (!SCIPlpiIsInfinity(lpi, ub))
        {
            ++*nvector;
            //*CVector = realloc(*CVector, sizeof(SCIP_Real) * ++*nvector);
            //(*CVector)[*nvector - 1] = ub;
            //SCIPdebugMessage("(*CVector)[%d]: %8.2f\n", *nvector - 1, (*CVector)[*nvector - 1]);
            //*AMatrixOfColumns = realloc(*AMatrixOfColumns, sizeof(SCIP_Real*) * *nvector);
            //(*AMatrixOfColumns)[*nvector - 1] = (SCIP_Real*)calloc(sizeof(SCIP_Real), ncols);
            //(*AMatrixOfColumns)[*nvector - 1][i] = get_column_obj_real(lpi, i);
        }
    }
    SCIPdebugMessage("*nvector: %d\n", *nvector);
    int nvector_ptr = 0;
    *AMatrixOfColumns = (scs_float**)calloc(sizeof(scs_float*), *nvector);
    *CVector = (scs_float**)calloc(sizeof(scs_float*), *nvector);
    for (int i = 0; i < ncols; i++)
    {
        scs_float lb = get_column_lower_bound_real(lpi, i);
        if (!SCIPlpiIsInfinity(lpi, -lb))
        {
            (*CVector)[nvector_ptr] = (scs_float*)calloc(1, sizeof(scs_float));
            (*CVector)[nvector_ptr][0] = -lb;
            SCIPdebugMessage("(*CVector)[%d]: %8.2f\n", nvector_ptr, (*CVector)[nvector_ptr][0]);
            (*AMatrixOfColumns)[nvector_ptr] = (scs_float*)calloc(sizeof(scs_float), ncols);
            (*AMatrixOfColumns)[nvector_ptr][i] = -1; //-get_column_obj_real(lpi, i);
            SCIPdebugMessage("(*AMatrixOfColumns)[%d][%d]: %8.2f\n", nvector_ptr, i,
                (*AMatrixOfColumns)[nvector_ptr][i]);
            nvector_ptr++;
        }
        scs_float ub = get_column_upper_bound_real(lpi, i);
        if (!SCIPlpiIsInfinity(lpi, ub))
        {
            (*CVector)[nvector_ptr] = (scs_float*)calloc(1, sizeof(scs_float));
            (*CVector)[nvector_ptr][0] = ub;
            SCIPdebugMessage("(*CVector)[%d]: %8.2f\n", nvector_ptr, (*CVector)[nvector_ptr][0]);
            (*AMatrixOfColumns)[nvector_ptr] = (scs_float*)calloc(sizeof(scs_float), ncols);
            (*AMatrixOfColumns)[nvector_ptr][i] = 1; // get_column_obj_real(lpi, i);
            SCIPdebugMessage("(*AMatrixOfColumns)[%d][%d]: %8.2f\n", nvector_ptr, i,
                (*AMatrixOfColumns)[nvector_ptr][i]);
            nvector_ptr++;
        }
    }
    return SCIP_OKAY;
}

SCIP_RETCODE ConstructAMatrixAndCVectorByRows(
    SCIP_LPI* lpi,
    scs_float*** AMatrixOfRows,
    scs_float*** CVector,
    int*         nvector
)
{
    assert(lpi != NULL);
    const int nrows = get_nrows(lpi);
    if (lpi->rows == NULL || nrows == 0)
    {
        return SCIP_OKAY;
    }
    *AMatrixOfRows = (scs_float**)calloc(sizeof(scs_float*), 0);
    *CVector = (scs_float**)calloc(sizeof(scs_float*), 0);
    *nvector = 0;
    for (int i = 0; i < get_nrows(lpi); i++)
    {
        const scs_float lhs = get_row_lhs_real(lpi, i);
        if (!SCIPlpiIsInfinity(lpi, -lhs))
        {
            scs_float* AVectorLHS = calloc(sizeof(scs_float), nrows);
            for (int j = 0; j < get_ncols(lpi); j++)
            {
                AVectorLHS[j] = -get_row_obj_real(lpi, i, j);
            }
            *AMatrixOfRows = realloc(*AMatrixOfRows, sizeof(scs_float*) * ++*nvector);
            (*AMatrixOfRows)[*nvector - 1] = AVectorLHS;
            *CVector = realloc(*CVector, sizeof(scs_float*) * *nvector);
            (*CVector)[*nvector - 1] = (scs_float*)calloc(1, sizeof(scs_float));
            (*CVector)[*nvector - 1][0] = -lhs;
        }
        const scs_float rhs = get_row_rhs_real(lpi, i);
        if (!SCIPlpiIsInfinity(lpi, rhs))
        {
            scs_float* AVectorRHS = calloc(sizeof(scs_float), nrows);
            for (int j = 0; j < get_ncols(lpi); j++)
            {
                AVectorRHS[j] = get_row_obj_real(lpi, i, j);
            }
            *AMatrixOfRows = realloc(*AMatrixOfRows, sizeof(scs_float*) * ++*nvector);
            (*AMatrixOfRows)[*nvector - 1] = AVectorRHS;
            *CVector = realloc(*CVector, sizeof(scs_float*) * *nvector);
            (*CVector)[*nvector - 1] = (scs_float*)calloc(1, sizeof(scs_float));
            (*CVector)[*nvector - 1][0] = rhs;
        }
    }
    return SCIP_OKAY;
}

SCIP_RETCODE debug_print_matrix_real(
    scs_float** matrix,
    scs_int row,
    scs_int col
)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            SCIPdebugMessage("matrix[%d][%d]: %8.2f", i, j, matrix[i][j]);
        }
        SCIPdebugMessage("\n");
    }
    return SCIP_OKAY;
}

SCIP_RETCODE CombineTwoMatricesByRow(
    scs_float** matrixA,
    scs_float** matrixB,
    scs_float*** matrix,
    const int nrowA,
    const int nrowB,
    const int ncol
)
{
    *matrix = (scs_float**)calloc(nrowA + nrowB, sizeof(scs_float*));
    for (int i = 0; i < nrowA; i++)
    {
        (*matrix)[i] = (scs_float*)calloc(ncol, sizeof(scs_float));
        for (int j = 0; j < ncol; j++)
        {
            (*matrix)[i][j] = matrixA[i][j];
        }
    }
    for (int i = nrowA; i < nrowA + nrowB; i++)
    {
        (*matrix)[i] = (scs_float*)calloc(ncol, sizeof(scs_float));
        for (int j = 0; j < ncol; j++)
        {
            (*matrix)[i][j] = matrixB[i - nrowA][j];
        }
    }
    return SCIP_OKAY;
}

SCIP_RETCODE CompressMatrixByRow(
    scs_float** matrix,
    const int row,
    const int col,
    scs_float** x,
    scs_int** ix,
    scs_int** p
)
{
    int nnonz = 0;
    *x = (scs_float*)calloc(nnonz, sizeof(scs_float));
    *ix = (scs_int*)calloc(nnonz, sizeof(scs_int));
    *p = (scs_int*)calloc(row + 1, sizeof(scs_int));
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if (ISLPIINFINITESIMAL(matrix[i][j]))
            {
                continue;
            }
            *x = realloc(*x, ++nnonz * sizeof(scs_float));
            (*x)[nnonz - 1] = matrix[i][j];
            *ix = realloc(*ix, nnonz * sizeof(scs_int));
            (*ix)[nnonz - 1] = j;
        }
        (*p)[i + 1] = nnonz;
    }
    /**
    int nnonzbycol = 0;
    for (int j = 0; j < col; j++)
    {
        for (int i = 0; i < row; i++)
        {
            if (ISLPIINFINITESIMAL(matrix[i][j]))
            {
                continue;
            }
            nnonzbycol++;
        }
        (*p)[j + 1] = nnonzbycol;
    }*/
    return SCIP_OKAY;
}

SCIP_RETCODE CompressMatrixByColumn(
    scs_float** matrix,
    const int row,
    const int col,
    scs_float** x,
    scs_int** ix,
    scs_int** p
)
{
    int nnonz = 0;
    *x = (scs_float*)calloc(nnonz, sizeof(scs_float));
    *ix = (scs_int*)calloc(nnonz, sizeof(scs_int));
    *p = (scs_int*)calloc(col + 1, sizeof(scs_int));
    for (int i = 0; i < col; i++)
    {
        for (int j = 0; j < row; j++)
        {
            if (ISLPIINFINITESIMAL(matrix[j][i]))
            {
                continue;
            }
            *x = realloc(*x, ++nnonz * sizeof(scs_float));
            (*x)[nnonz - 1] = matrix[j][i];
            *ix = realloc(*ix, nnonz * sizeof(scs_int));
            (*ix)[nnonz - 1] = j;
        }
        (*p)[i + 1] = nnonz;
    }
    return SCIP_OKAY;
}

SCIP_RETCODE InverseMatrix(
    scs_float** origin,
    scs_float*** result,
    int row,
    int col
)
{
    *result = (scs_float**)calloc(col, sizeof(scs_float*));
    for (int j = 0; j < col; j++)
    {
        (*result)[j] = (scs_float*)calloc(row, sizeof(scs_float));
    }
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            (*result)[j][i] = origin[i][j];
        }
    }
    return SCIP_OKAY;
}

SCIP_RETCODE ConstructPMatrix(
    scs_float** Px,
    scs_int** Pi,
    scs_int** Pp,
    int n
)
{
    scs_float** matrix = calloc(n, sizeof(scs_float*));
    for (int i = 0; i < n; i++)
    {
        matrix[i] = (scs_float*)calloc(n, sizeof(scs_float));
    }
    return CompressMatrixByColumn(matrix, n, n, &*Px, &*Pi, &*Pp);
}

SCIP_RETCODE ConstructCVector(
    SCIP_LPI* lpi,
    scs_float** c
)
{
    *c = (scs_float*)calloc(get_ncols(lpi), sizeof(scs_float));
    for (int i = 0; i < get_ncols(lpi); i++)
    {
        (*c)[i] = get_column_obj_real(lpi, i);
    }
    return SCIP_OKAY;
}

SCIP_RETCODE ConstructAMatrix(
    SCIP_LPI* lpi,
    scs_float** Ax,
    scs_int** Ai,
    scs_int** Ap,
    scs_float** b,
    scs_float** c,
    int*  m,
    int*  n
)
{
    scs_float** AMatrixOfColumns;
    scs_float** CVectorOfColumns;
    int nvectorCol = 0;
    ConstructAMatrixAndCVectorByColumns(lpi, &AMatrixOfColumns, &CVectorOfColumns, &nvectorCol);
    assert(AMatrixOfColumns != NULL);
    assert(CVectorOfColumns != NULL);

    debug_print_matrix_real(AMatrixOfColumns, nvectorCol, get_ncols(lpi));
    debug_print_matrix_real(CVectorOfColumns, nvectorCol, 1);

    scs_float** AMatrixOfRows;
    scs_float** CVectorOfRows;
    int nvectorRow = 0;
    ConstructAMatrixAndCVectorByRows(lpi, &AMatrixOfRows, &CVectorOfRows, &nvectorRow);

    debug_print_matrix_real(AMatrixOfRows, nvectorRow, get_ncols(lpi));
    debug_print_matrix_real(CVectorOfRows, nvectorRow, 1);

    *m = nvectorCol + nvectorRow;
    *n = get_ncols(lpi);

    scs_float** AMatrix;
    scs_float** CVector;

    CombineTwoMatricesByRow(AMatrixOfColumns, AMatrixOfRows, &AMatrix, nvectorCol, nvectorRow, *n);
    lpi->nconsbycol = nvectorCol;
    debug_print_matrix_real(AMatrix, *m, *n);
    CombineTwoMatricesByRow(CVectorOfColumns, CVectorOfRows, &CVector, nvectorCol, nvectorRow, 1);
    debug_print_matrix_real(CVector, *m, 1);

    CompressMatrixByColumn(AMatrix, *m, *n, &*Ax, &*Ai, &*Ap);

    scs_float** Ib;
    InverseMatrix(CVector, &Ib, *m, 1);
    *b = Ib[0];
    ConstructCVector(lpi, &*c);
    return SCIP_OKAY;
}

SCIP_RETCODE ConstructScsData(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    scs_float* Ax = NULL;
    scs_int* Ai = NULL;
    scs_int* Ap = NULL;
    int m = 0;
    int n = 0;
    scs_float* b = NULL;
    scs_float* c = NULL;
    ConstructAMatrix(lpi, &Ax, &Ai, &Ap, &b, &c, &m, &n);
    assert(Ax != NULL);
    assert(Ai != NULL);
    assert(Ap != NULL);
    
    lpi->scsdata->m = m;
    lpi->scsdata->n = n;
    scs_float* Px = NULL;
    scs_int* Pi = NULL;
    scs_int* Pp = NULL;
    ConstructPMatrix(&Px, &Pi, &Pp, n);
    assert(Px != NULL);
    assert(Pi != NULL);
    assert(Pp != NULL);
    lpi->scsdata->b = b;
    lpi->scsdata->c = c;
    ScsMatrix A = { Ax, Ai, Ap, m, n };
    lpi->scsdata->A = (ScsMatrix*)calloc(1, sizeof(A));
    memcpy(lpi->scsdata->A, &A, sizeof(A));
    ScsMatrix P = { Px, Pi, Pp, n, n};
    lpi->scsdata->P = (ScsMatrix*)calloc(1, sizeof(P));
    memcpy(lpi->scsdata->P, &P, sizeof(P));
    return SCIP_OKAY;
}

SCIP_RETCODE debug_print_scs_data(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    assert(lpi->scsdata != NULL);
    assert(lpi->scsdata->m > 0);
    assert(lpi->scsdata->n > 0);
    assert(lpi->scsdata->A != NULL);
    assert(lpi->scsdata->A->m == lpi->scsdata->m && lpi->scsdata->A->n == lpi->scsdata->n);
    assert(lpi->scsdata->P != NULL);
    assert(lpi->scsdata->P->m == lpi->scsdata->n && lpi->scsdata->P->n == lpi->scsdata->n);
    assert(lpi->scsdata->b != NULL);
    assert(lpi->scsdata->c != NULL);
    SCIPdebugMessage("SCSData P matrix:\n");
    scs_int nnonz = lpi->scsdata->P->p[lpi->scsdata->P->n];
    if (nnonz) {
        for (int i = 0; i < nnonz; i++)
        {
            SCIPdebugMessage("%8.2f at %lld, ", lpi->scsdata->P->x[i], lpi->scsdata->P->i[i]);
        }
        SCIPdebugMessage("\n");
    } else
    {
        SCIPdebugMessage("P matrix is empty.\n");
    }

    SCIPdebugMessage("SCSData A matrix:\n");
    nnonz = lpi->scsdata->A->p[lpi->scsdata->A->n];
    if (nnonz) {
        for (int i = 0; i < nnonz; i++)
        {
            SCIPdebugMessage("%8.2f at %lld ", lpi->scsdata->A->x[i], lpi->scsdata->A->i[i]);
        }
        SCIPdebugMessage("\n");
    } else
    {
        SCIPdebugMessage("A matrix is empty.\n");
    }
    return SCIP_OKAY;
}

SCIP_RETCODE debug_print_scs_solution(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    SCIPdebugMessage("Primal objective: [%d]%8.2f\n", lpi->objsen, lpi->scsinfo->pobj);
    SCIPdebugMessage("Dual objective: [%d]%8.2f\n", lpi->objsen, lpi->scsinfo->dobj);
    assert(lpi->scssol != NULL);
    if (lpi->scssol->s)
    {
        SCIPdebugMessage("Slack variables:\n");
	    for (int i = 0; i < get_ncols(lpi); i++)
	    {
            SCIPdebugMessage("s[%d]: %8.2f ", i, lpi->scssol->s[i]);
	    }
        SCIPdebugMessage("\n");
    }
    if (lpi->scssol->x)
    {
        SCIPdebugMessage("Primal Solutions:\n");
	    for (int i = 0; i < get_ncols(lpi); i++)
	    {
            SCIPdebugMessage("x[%d]: %8.2f ", i, lpi->scssol->x[i]);
	    }
        SCIPdebugMessage("\n");
    }
    if (lpi->scssol->y)
    {
        SCIPdebugMessage("Dual Solutions:\n");
	    for (int i = 0; i < get_nrows(lpi); i++)
	    {
            SCIPdebugMessage("y[%d]: %8.2f ", i, lpi->scssol->y[i + lpi->nconsbycol]);
	    }
        SCIPdebugMessage("\n");
    }
    return SCIP_OKAY;
}

SCIP_RETCODE scsSolve(
    SCIP_LPI* lpi
)
{
    assert(lpi != NULL);
    ConstructScsData(lpi);
    debug_print_scs_data(lpi);
    lpi->scscone->z = 0;
    lpi->scscone->l = lpi->scsdata->m;
    lpi->scswork = scs_init(lpi->scsdata, lpi->scscone, lpi->scsstgs);
    scs_int exitflag = scs_solve(lpi->scswork, lpi->scssol, lpi->scsinfo, 0);
    debug_print_scs_solution(lpi);
    scs_finish(lpi->scswork);
    return SCIP_OKAY;
}

 /** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return scsSolve(lpi);
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return scsSolve(lpi);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
    SCIP_LPI* lpi,                /**< LP interface structure */
    SCIP_Bool             crossover           /**< perform crossover */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    return SCIPlpiSolveDual(lpi);
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
    return lpi->scsinfo->status_val == SCS_SOLVED || lpi->scsinfo->status_val == SCS_SOLVED_INACCURATE;
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
    assert(lpi->scsinfo != NULL);
    assert(primalfeasible != NULL);
    assert(dualfeasible != NULL);
    *primalfeasible = !SCIPlpiIsInfinity(lpi, lpi->scsinfo->pobj); 
    *dualfeasible = !SCIPlpiIsInfinity(lpi, lpi->scsinfo->dobj);
    return SCIP_OKAY;
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
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val == SCS_UNBOUNDED;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return SCIPlpiIsInfinity(lpi, lpi->scsinfo->pobj);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val == SCS_SOLVED || lpi->scsinfo->status_val == SCS_SOLVED_INACCURATE;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
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
    return FALSE;
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val == SCS_UNBOUNDED;
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return SCIPlpiIsInfinity(lpi, lpi->scsinfo->dobj);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val == SCS_SOLVED || lpi->scsinfo->status_val == SCS_SOLVED_INACCURATE;
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val == SCS_SOLVED;
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
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val == SCS_INFEASIBLE || lpi->scsinfo->status_val == SCS_UNBOUNDED ||
        lpi->scsinfo->status_val == SCS_SOLVED || lpi->scsinfo->status_val == SCS_FAILED;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val == SCS_SIGINT || lpi->scsinfo->status_val == SCS_FAILED
    || lpi->scsinfo->status_val == SCS_UNFINISHED;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val == SCS_SIGINT || lpi->scsinfo->status_val == SCS_FAILED
    || lpi->scsinfo->status_val == SCS_UNFINISHED;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val == SCS_SIGINT || lpi->scsinfo->status_val == SCS_FAILED
    || lpi->scsinfo->status_val == SCS_UNFINISHED;
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
    SCIP_LPI* lpi                 /**< LP interface structure */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    assert(lpi->scsinfo != NULL);
    return lpi->scsinfo->status_val;
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
    assert(lpi->scsinfo != NULL);
    assert(objval != NULL);
    *objval = lpi->objsen == SCIP_OBJSEN_MINIMIZE ? lpi->scsinfo->pobj : -lpi->scsinfo->pobj;
    return SCIP_OKAY;
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
    assert(lpi->scsinfo != NULL);
    assert(lpi->scssol != NULL);
    if (objval) {
        SCIPlpiGetObjval(lpi, &*objval);
    }
    if (primsol) {
        for (int i = 0; i < get_ncols(lpi); i++)
        {
            primsol[i] = lpi->scssol->x[i];
        }
    }
    /**
    dualsol = NULL;
    activity = NULL;
    redcost = NULL;
	*/
    if (dualsol) {
        // @TODO 转换方式存疑。
        /**
        for (int i = lpi->nconsbycol; i < get_nrows(lpi); i++)
        {
            dualsol[i] = lpi->scssol->y[i];
        }*/
        dualsol = (SCIP_Real*)calloc(get_nrows(lpi) - lpi->nconsbycol, sizeof(SCIP_Real));
    }
    if (activity)
    {
	    for (int i = lpi->nconsbycol; i< get_nrows(lpi); i++)
	    {
		    
	    }
    }
    if (redcost)
    {
        for (int i = 0; i < get_ncols(lpi); i++)
        {
            redcost[i] = 0;
        }
    }
    return SCIP_OKAY;
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
    assert(lpi->scsinfo != NULL);
    assert(iterations != NULL);
    *iterations = lpi->scsinfo->iter;
    return SCIP_OKAY;
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
    SCIPdebugMessage("calling SCIPlpiGetBase()...\n");
    assert(lpi != NULL);
    if (rstat != NULL)
    {
        for (int i = 0; i < get_nrows(lpi); i++)
        {
            rstat[i] = lpi->rstat[i];
        }
    }
    if (cstat != NULL)
    {
        for (int i = 0; i < get_ncols(lpi); i++)
        {
            cstat[i] = lpi->cstat[i];
        }
    }
    return SCIP_OKAY;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
    SCIP_LPI* lpi,                /**< LP interface structure */
    const int* cstat,              /**< array with column basis status */
    const int* rstat               /**< array with row basis status */
)
{  /*lint --e{715}*/
    SCIPdebugMessage("calling SCIPlpiSetBase()...\n");
    assert(lpi != NULL);
    const int ncols = get_ncols(lpi);
    const int nrows = get_nrows(lpi);
    assert(cstat != NULL || ncols == 0);
    assert(rstat != NULL || nrows == 0);
    SCIP_CALL(resize_state_columns(lpi, ncols));
    SCIP_CALL(resize_state_rows(lpi, nrows));
    if (rstat != NULL)
    {
        for (int i = 0; i < get_nrows(lpi); i++)
        {
            lpi->rstat[i] = rstat[i];
        }
    }
    if (cstat != NULL)
    {
        for (int i = 0; i < get_ncols(lpi); i++)
        {
            lpi->cstat[i] = cstat[i];
        }
    }
    return SCIP_OKAY;
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
 * dynamic memory arrays
 */

 /** resizes cstat array to have at least num entries */
static
SCIP_RETCODE ensureCstatMem(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   num                 /**< minimal number of entries in array */
)
{
    assert(lpi != NULL);

    if (num > lpi->cstatsize)
    {
        int newsize;

        newsize = MAX(2 * lpi->cstatsize, num);
        SCIP_ALLOC(BMSreallocMemoryArray(&lpi->cstat, newsize));
        lpi->cstatsize = newsize;
    }
    assert(num <= lpi->cstatsize);

    return SCIP_OKAY;
}

/** resizes rstat array to have at least num entries */
static
SCIP_RETCODE ensureRstatMem(
    SCIP_LPI* lpi,                /**< LP interface structure */
    int                   num                 /**< minimal number of entries in array */
)
{
    assert(lpi != NULL);

    if (num > lpi->rstatsize)
    {
        int newsize;

        newsize = MAX(2 * lpi->rstatsize, num);
        SCIP_ALLOC(BMSreallocMemoryArray(&lpi->rstat, newsize));
        lpi->rstatsize = newsize;
    }
    assert(num <= lpi->rstatsize);

    return SCIP_OKAY;
}

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
    SCIPdebugMessage("calling SCIPlpiGetState()...\n");
    assert(lpi != NULL);
    assert(blkmem != NULL);
    assert(lpistate != NULL);
    assert(blkmem != NULL);
    const int ncols = get_ncols(lpi);
    const int nrows = get_nrows(lpi);
    assert(ncols >= 0);
    assert(nrows >= 0);
    SCIP_CALL(lpistateCreate(lpistate, blkmem, ncols, nrows));
    SCIP_CALL(SCIPlpiGetBase(lpi, lpi->cstat, lpi->rstat));
    (*lpistate)->ncols = ncols;
    (*lpistate)->nrows = nrows;
    lpistatePack(*lpistate, lpi->cstat, lpi->rstat);
    return SCIP_OKAY;
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
    SCIPdebugMessage("calling SCIPlpiSetState()...\n");

    assert(lpi != NULL);
    /** assert(blkmem != NULL); */
    assert(lpistate != NULL);

    int lpncols = get_ncols(lpi);
    int lpnrows = get_nrows(lpi);
    assert(lpistate->ncols <= lpncols);
    assert(lpistate->nrows <= lpnrows);
    SCIP_CALL(ensureCstatMem(lpi, lpncols));
    SCIP_CALL(ensureRstatMem(lpi, lpnrows));
    lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);
    for (int i = lpistate->ncols; i < lpncols; i++)
    {
        SCIP_Real lb = get_column_lower_bound_real(lpi, i);
        if (SCIPlpiIsInfinity(lpi, REALABS(lb)))
        {
            lb = get_column_lower_bound_real(lpi, i);
            if (SCIPlpiIsInfinity(lpi, REALABS(lb)))
            {
                lpi->cstat[i] = SCIP_BASESTAT_ZERO;
            } else
            {
                lpi->cstat[i] = SCIP_BASESTAT_UPPER;
            }
        } else
        {
            lpi->cstat[i] = SCIP_BASESTAT_LOWER;
        }
    }
    for (int i = lpistate->nrows; i < lpnrows; i++)
    {
        lpi->rstat[i] = SCIP_BASESTAT_BASIC;
    }
    SCIP_CALL(SCIPlpiSetBase(lpi, lpi->cstat, lpi->rstat));
    return SCIP_OKAY;
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
    if (*lpistate != NULL)
    {
        lpistateFree(lpistate, blkmem);
    }
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
    //errorMessage();
    return SCIP_OKAY;
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
    //errorMessage();
    return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
    SCIP_LPI* lpi,                /**< LP interface structure */
    BMS_BLKMEM* blkmem,             /**< block memory */
    SCIP_LPINORMS** lpinorms            /**< pointer to LPi pricing norms information, or NULL */
)
{  /*lint --e{715}*/
    assert(lpi != NULL);
    //errorMessage();
    return SCIP_OKAY;
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
    switch (type)
    {
    case SCIP_LPPAR_FROMSCRATCH:
        *ival = lpi->fromscratch;
        break;
    case SCIP_LPPAR_REFACTOR:
        *ival = lpi->refactor;
        break;
    case SCIP_LPPAR_LPINFO:
        *ival = lpi->lpinfo;
        break;
    case SCIP_LPPAR_LPITLIM:
        *ival = lpi->lpitlim;
        /* -1 <= ival, -1 meaning no time limit, 0 stopping immediately */
        if (*ival == -1)
            *ival = INT_MAX;
        break;
    case SCIP_LPPAR_PRESOLVING:
        *ival = lpi->presoving;
        break;
    case SCIP_LPPAR_PRICING:
        *ival = (int)lpi->pricing;
        break;
    case SCIP_LPPAR_SCALING:
        *ival = lpi->scaling;
        break;
    case SCIP_LPPAR_TIMING:
        *ival = lpi->timing;
        break;
    case SCIP_LPPAR_RANDOMSEED:
        *ival = lpi->randomseed;
        break;
    case SCIP_LPPAR_POLISHING:
        *ival = lpi->polishing;
        break;
    default:
        return SCIP_PARAMETERUNKNOWN;
    }  /*lint !e788*/

    return SCIP_OKAY;
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
        lpi->fromscratch = (SCIP_Bool)ival;
        break;
    case SCIP_LPPAR_LPINFO:
        assert(ival == TRUE || ival == FALSE);
        lpi->lpinfo = (SCIP_Bool)ival;
        break;
    case SCIP_LPPAR_LPITLIM:
        assert(ival >= 0);
        if (ival >= INT_MAX)
        {
            ival = -1;
        }
        lpi->lpitlim = ival;
        break;
    case SCIP_LPPAR_PRESOLVING:
        assert(ival == TRUE || ival == FALSE);
        lpi->presoving = ival;
        break;
    case SCIP_LPPAR_PRICING:
        lpi->pricing = (SCIP_PRICING)ival;
        // switch (lpi->pricing) : lpi->pricer = ...
        break;
    case SCIP_LPPAR_SCALING:
        assert(ival >= 0 && ival <= 2);
        lpi->scaling = ival;
        break;
    case SCIP_LPPAR_TIMING:
        assert(ival >= 0 && ival < 3);
        lpi->timing = ival;
        break;
    case SCIP_LPPAR_RANDOMSEED:
        lpi->randomseed = (unsigned long)ival;
        break;
    case SCIP_LPPAR_POLISHING:
        assert(ival >= 0 && ival < 3);
        lpi->polishing = ival;
        break;
    case SCIP_LPPAR_REFACTOR:
        assert(ival >= 0);
        lpi->refactor = ival;
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
    switch (type)
    {
    case SCIP_LPPAR_FEASTOL:
        *dval = lpi->feastol;
        break;
    case SCIP_LPPAR_DUALFEASTOL:
        *dval = lpi->dualfeastol;
        break;
    case SCIP_LPPAR_OBJLIM:
        *dval = lpi->objlim;
        break;
    case SCIP_LPPAR_LPTILIM:
        *dval = lpi->lptilim;
        break;
    case SCIP_LPPAR_ROWREPSWITCH:
        *dval = lpi->rowrepswitch;
        if (*dval >= SCIPlpiInfinity(lpi))
        {
            *dval = -1.0;
        }
        break;
    case SCIP_LPPAR_CONDITIONLIMIT:
        *dval = lpi->conditionlimit;
        break;
    case SCIP_LPPAR_MARKOWITZ:
        /**
        *dval = lpi->markowitz;
        break;*/
    default:
        return SCIP_PARAMETERUNKNOWN;
    }  /*lint !e788*/
    return SCIP_OKAY;
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
        assert(dval > 0.0);
        lpi->feastol = dval;
        break;
    case SCIP_LPPAR_DUALFEASTOL:
        assert(dval > 0.0);
        lpi->dualfeastol = dval;
        break;
    case SCIP_LPPAR_OBJLIM:
        lpi->objlim = dval;
        break;
    case SCIP_LPPAR_LPTILIM:
        assert(dval > 0.0);
        lpi->lptilim = dval;
        break;
    case SCIP_LPPAR_ROWREPSWITCH:
        assert(dval > 0.0 || ISLPIINFINITESIMAL(dval) || ISLPIINFINITESIMAL(dval + 1.0));
        if (ISLPIINFINITESIMAL(dval + 1))
        {
            lpi->rowrepswitch = SCIPlpiInfinity(lpi);
        } else
        {
            lpi->rowrepswitch = dval;
        }
        break;
    case SCIP_LPPAR_CONDITIONLIMIT:
        lpi->conditionlimit = dval;
        lpi->checkcondition = dval > 0.0 || ISLPIINFINITESIMAL(dval);
        break;
    case SCIP_LPPAR_BARRIERCONVTOL:
        /**
        if (dval < 1e-4)
        {
            dval = 1e-4;
        } else if (dval > 0.9999)
        {
            dval = 0.9999;
        }
        lpi->markowitz = dval;
        break;*/
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
