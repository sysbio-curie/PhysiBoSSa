/* A Bison parser, made by GNU Bison 3.4.1.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2019 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with RC or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 2 "RunConfigGrammar.y"

/* 
   MaBoSS (Markov Boolean Stochastic Simulator)
   Copyright (C) 2011 Institut Curie, 26 rue d'Ulm, Paris, France
   
   MaBoSS is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   
   MaBoSS is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA 
*/

/*
   Module:
     RunConfigGrammar.l

   Authors:
     Eric Viara <viara@sysra.com>
     Gautier Stoll <gautier.stoll@curie.fr>
 
   Date:
     January-March 2011
     updated to add extended initial conditions: October 2014
*/

#include "BooleanNetwork.h"
#include "RunConfig.h"

extern int RClex();

static void RCerror(const char *s);
static Network* network;
static RunConfig* config;

extern std::string RC_error_head();

#line 115 "RunConfigGrammar.tab.c"

# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif


/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int RCdebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum RCtokentype
  {
    VARIABLE = 258,
    SYMBOL = 259,
    DOUBLE = 260,
    INTEGER = 261,
    LOGAND = 262,
    LOGOR = 263,
    LOGXOR = 264,
    LOGNOT = 265,
    EQUAL = 266,
    NOT_EQUAL = 267,
    NODE = 268,
    GTEQ = 269,
    LTEQ = 270
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 47 "RunConfigGrammar.y"

  Expression* expr;
  char* str;
  double d;
  long long l;
  std::vector<const Node*>* node_list;
  std::vector<Expression*>* expr_list;
  IStateGroup::ProbaIState* istate_expr;
  std::vector<IStateGroup::ProbaIState*>* istate_expr_list;
  ArgumentList* arg_list;


#line 184 "RunConfigGrammar.tab.c"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE RClval;

int RCparse (void);





#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 RCtype_uint8;
#else
typedef unsigned char RCtype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 RCtype_int8;
#else
typedef signed char RCtype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 RCtype_uint16;
#else
typedef unsigned short RCtype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 RCtype_int16;
#else
typedef short RCtype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about RClval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if ! defined RCoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined RCoverflow || YYERROR_VERBOSE */


#if (! defined RCoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union RCalloc
{
  RCtype_int16 RCss_alloc;
  YYSTYPE RCvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union RCalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (RCtype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T RCnewbytes;                                            \
        YYCOPY (&RCptr->Stack_alloc, Stack, RCsize);                    \
        Stack = &RCptr->Stack_alloc;                                    \
        RCnewbytes = RCstacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        RCptr += RCnewbytes / sizeof (*RCptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T RCi;                         \
          for (RCi = 0; RCi < (Count); RCi++)   \
            (Dst)[RCi] = (Src)[RCi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  16
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   102

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  34
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  24
/* YYNRULES -- Number of rules.  */
#define YYNRULES  57
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  109

#define YYUNDEFTOK  2
#define YYMAXUTOK   270

/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by RClex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                                \
  ((unsigned) (YYX) <= YYMAXUTOK ? RCtranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by RClex.  */
static const RCtype_uint8 RCtranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    26,     2,     2,     2,     2,     2,     2,
      22,    23,    27,    24,    21,    25,    18,    28,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    33,    16,
      29,    17,    30,    32,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    19,     2,    20,    31,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const RCtype_uint16 RCrline[] =
{
       0,    88,    88,    91,    96,    97,    98,    99,   102,   112,
     143,   158,   164,   170,   177,   182,   189,   200,   205,   212,
     222,   226,   230,   235,   241,   246,   254,   258,   263,   271,
     275,   279,   283,   287,   293,   297,   301,   307,   311,   315,
     321,   325,   329,   333,   337,   343,   347,   351,   357,   361,
     367,   371,   377,   381,   385,   391,   395,   401
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const RCtname[] =
{
  "$end", "error", "$undefined", "VARIABLE", "SYMBOL", "DOUBLE",
  "INTEGER", "LOGAND", "LOGOR", "LOGXOR", "LOGNOT", "EQUAL", "NOT_EQUAL",
  "NODE", "GTEQ", "LTEQ", "';'", "'='", "'.'", "'['", "']'", "','", "'('",
  "')'", "'+'", "'-'", "'!'", "'*'", "'/'", "'<'", "'>'", "'^'", "'?'",
  "':'", "$accept", "translation_unit", "decl", "runconfig_decl",
  "node_attr_decl", "symbol_istate_list", "symbol_list",
  "istate_expression_list", "istate_expression", "expression_list",
  "var_decl", "primary_expression", "argument_expression_list",
  "postfix_expression", "unary_expression", "multiplicative_expression",
  "additive_expression", "relational_expression", "equality_expression",
  "logical_and_expression", "logical_or_expression",
  "logical_xor_expression", "conditional_expression", "expression", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const RCtype_uint16 RCtoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,    59,    61,    46,    91,
      93,    44,    40,    41,    43,    45,    33,    42,    47,    60,
      62,    94,    63,    58
};
# endif

#define YYPACT_NINF -70

#define RCpact_value_is_default(Yystate) \
  (!!((Yystate) == (-70)))

#define YYTABLE_NINF -1

#define RCtable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const RCtype_int8 RCpact[] =
{
      45,    -6,    -3,   -70,    46,    13,   -70,   -70,   -70,     8,
     -70,    41,    41,    68,   -70,    16,   -70,   -70,    83,   -70,
       6,   -70,   -70,    41,    41,    41,    41,    41,   -70,   -70,
     -70,    42,    50,   -11,    65,    55,    80,    -7,   -70,    73,
      74,    75,   -70,    87,    76,    17,   -70,    71,   -70,   -70,
     -70,    41,    41,    41,    41,    41,    41,    41,    41,    41,
      41,    41,    41,    41,    41,    41,   -70,   -70,    41,   -70,
      49,   -70,   -16,   -70,   -70,   -70,   -70,    42,    42,    50,
      50,    50,    50,   -11,   -11,    65,    55,    80,    80,    62,
      81,    52,   -70,    77,    41,   -70,    41,   -70,   -70,    49,
      49,   -70,   -70,   -70,    58,   -70,   -70,    49,   -70
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const RCtype_uint8 RCdefact[] =
{
       0,     0,     0,     7,     0,     0,     2,     5,     6,     0,
       4,     0,     0,     0,    12,     0,     1,     3,     0,    22,
       0,    21,    20,     0,     0,     0,     0,     0,    26,    29,
      34,    37,    40,    45,    48,    50,    52,    55,    57,     0,
       0,     0,    11,     0,     0,     0,    33,     0,    30,    31,
      32,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    19,     8,     0,    13,
       0,    28,     0,    24,    23,    35,    36,    38,    39,    44,
      43,    41,    42,    46,    47,    49,    51,    53,    54,     0,
       0,     0,    14,     0,     0,    27,     0,     9,    10,     0,
       0,    25,    56,    15,     0,    17,    16,     0,    18
};

  /* YYPGOTO[NTERM-NUM].  */
static const RCtype_int8 RCpgoto[] =
{
     -70,   -70,    93,   -70,   -70,   -70,   -70,   -70,     0,   -70,
     -70,   -69,   -70,   -70,   -17,    27,     2,    24,    39,    40,
      22,   -70,     5,   -12
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const RCtype_int8 RCdefgoto[] =
{
      -1,     5,     6,     7,     8,     9,    15,    91,    92,   104,
      10,    28,    72,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const RCtype_uint8 RCtable[] =
{
      40,    93,    63,    55,    56,    94,    46,    95,    48,    49,
      50,    11,    47,    16,    12,    13,     1,     2,    57,    58,
      19,    20,    21,    22,    64,    65,    18,    23,    45,     3,
      93,   105,     4,    73,    75,    76,    42,    43,   108,    24,
      71,    25,    26,    27,    19,    20,    21,    22,     1,     2,
      14,    23,    19,    89,    21,    22,    90,    79,    80,    81,
      82,     3,    61,    24,     4,    25,    26,    27,    98,    51,
      52,    24,    41,    99,    53,    54,    59,    60,   106,   107,
      77,    78,   101,    83,    84,    87,    88,    44,    62,    66,
      67,    69,    68,    70,    74,    96,   100,    97,    17,   103,
      85,   102,    86
};

static const RCtype_uint8 RCcheck[] =
{
      12,    70,     9,    14,    15,    21,    23,    23,    25,    26,
      27,    17,    24,     0,    17,    18,     3,     4,    29,    30,
       3,     4,     5,     6,    31,    32,    18,    10,    22,    16,
      99,   100,    19,    45,    51,    52,    20,    21,   107,    22,
      23,    24,    25,    26,     3,     4,     5,     6,     3,     4,
       4,    10,     3,    65,     5,     6,    68,    55,    56,    57,
      58,    16,     7,    22,    19,    24,    25,    26,    16,    27,
      28,    22,     4,    21,    24,    25,    11,    12,    20,    21,
      53,    54,    94,    59,    60,    63,    64,     4,     8,    16,
      16,     4,    17,    17,    23,    33,    19,    16,     5,    99,
      61,    96,    62
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const RCtype_uint8 RCstos[] =
{
       0,     3,     4,    16,    19,    35,    36,    37,    38,    39,
      44,    17,    17,    18,     4,    40,     0,    36,    18,     3,
       4,     5,     6,    10,    22,    24,    25,    26,    45,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    56,    57,
      57,     4,    20,    21,     4,    22,    48,    57,    48,    48,
      48,    27,    28,    24,    25,    14,    15,    29,    30,    11,
      12,     7,     8,     9,    31,    32,    16,    16,    17,     4,
      17,    23,    46,    57,    23,    48,    48,    49,    49,    50,
      50,    50,    50,    51,    51,    52,    53,    54,    54,    57,
      57,    41,    42,    45,    21,    23,    33,    16,    16,    21,
      19,    57,    56,    42,    43,    45,    20,    21,    45
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const RCtype_uint8 RCr1[] =
{
       0,    34,    35,    35,    36,    36,    36,    36,    37,    38,
      38,    39,    40,    40,    41,    41,    42,    43,    43,    44,
      45,    45,    45,    45,    46,    46,    47,    47,    47,    48,
      48,    48,    48,    48,    49,    49,    49,    50,    50,    50,
      51,    51,    51,    51,    51,    52,    52,    52,    53,    53,
      54,    54,    55,    55,    55,    56,    56,    57
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const RCtype_uint8 RCr2[] =
{
       0,     2,     1,     2,     1,     1,     1,     1,     4,     6,
       6,     3,     1,     3,     1,     3,     4,     1,     3,     4,
       1,     1,     1,     3,     1,     3,     1,     4,     3,     1,
       2,     2,     2,     2,     1,     3,     3,     1,     3,     3,
       1,     3,     3,     3,     3,     1,     3,     3,     1,     3,
       1,     3,     1,     3,     3,     1,     5,     1
};


#define RCerrok         (RCerrstatus = 0)
#define RCclearin       (RCchar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto RCacceptlab
#define YYABORT         goto RCabortlab
#define YYERROR         goto RCerrorlab


#define YYRECOVERING()  (!!RCerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (RCchar == YYEMPTY)                                        \
      {                                                           \
        RCchar = (Token);                                         \
        RClval = (Value);                                         \
        YYPOPSTACK (RClen);                                       \
        RCstate = *RCssp;                                         \
        goto RCbackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        RCerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (RCdebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (RCdebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      RC_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
RC_symbol_value_print (FILE *RCo, int RCtype, YYSTYPE const * const RCvaluep)
{
  FILE *RCoutput = RCo;
  YYUSE (RCoutput);
  if (!RCvaluep)
    return;
# ifdef YYPRINT
  if (RCtype < YYNTOKENS)
    YYPRINT (RCo, RCtoknum[RCtype], *RCvaluep);
# endif
  YYUSE (RCtype);
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
RC_symbol_print (FILE *RCo, int RCtype, YYSTYPE const * const RCvaluep)
{
  YYFPRINTF (RCo, "%s %s (",
             RCtype < YYNTOKENS ? "token" : "nterm", RCtname[RCtype]);

  RC_symbol_value_print (RCo, RCtype, RCvaluep);
  YYFPRINTF (RCo, ")");
}

/*------------------------------------------------------------------.
| RC_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
RC_stack_print (RCtype_int16 *RCbottom, RCtype_int16 *RCtop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; RCbottom <= RCtop; RCbottom++)
    {
      int RCbot = *RCbottom;
      YYFPRINTF (stderr, " %d", RCbot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (RCdebug)                                                  \
    RC_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
RC_reduce_print (RCtype_int16 *RCssp, YYSTYPE *RCvsp, int RCrule)
{
  unsigned long RClno = RCrline[RCrule];
  int RCnrhs = RCr2[RCrule];
  int RCi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             RCrule - 1, RClno);
  /* The symbols being reduced.  */
  for (RCi = 0; RCi < RCnrhs; RCi++)
    {
      YYFPRINTF (stderr, "   $%d = ", RCi + 1);
      RC_symbol_print (stderr,
                       RCstos[RCssp[RCi + 1 - RCnrhs]],
                       &RCvsp[(RCi + 1) - (RCnrhs)]
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (RCdebug)                          \
    RC_reduce_print (RCssp, RCvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int RCdebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef RCstrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define RCstrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
RCstrlen (const char *RCstr)
{
  YYSIZE_T RClen;
  for (RClen = 0; RCstr[RClen]; RClen++)
    continue;
  return RClen;
}
#  endif
# endif

# ifndef RCstpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define RCstpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
RCstpcpy (char *RCdest, const char *RCsrc)
{
  char *RCd = RCdest;
  const char *RCs = RCsrc;

  while ((*RCd++ = *RCs++) != '\0')
    continue;

  return RCd - 1;
}
#  endif
# endif

# ifndef RCtnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for RCerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from RCtname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
RCtnamerr (char *RCres, const char *RCstr)
{
  if (*RCstr == '"')
    {
      YYSIZE_T RCn = 0;
      char const *RCp = RCstr;

      for (;;)
        switch (*++RCp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++RCp != '\\')
              goto do_not_strip_quotes;
            else
              goto append;

          append:
          default:
            if (RCres)
              RCres[RCn] = *RCp;
            RCn++;
            break;

          case '"':
            if (RCres)
              RCres[RCn] = '\0';
            return RCn;
          }
    do_not_strip_quotes: ;
    }

  if (! RCres)
    return RCstrlen (RCstr);

  return (YYSIZE_T) (RCstpcpy (RCres, RCstr) - RCres);
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
RCsyntax_error (YYSIZE_T *RCmsg_alloc, char **RCmsg,
                RCtype_int16 *RCssp, int RCtoken)
{
  YYSIZE_T RCsize0 = RCtnamerr (YY_NULLPTR, RCtname[RCtoken]);
  YYSIZE_T RCsize = RCsize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *RCformat = YY_NULLPTR;
  /* Arguments of RCformat. */
  char const *RCarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int RCcount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in RCchar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated RCchar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (RCtoken != YYEMPTY)
    {
      int RCn = RCpact[*RCssp];
      RCarg[RCcount++] = RCtname[RCtoken];
      if (!RCpact_value_is_default (RCn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int RCxbegin = RCn < 0 ? -RCn : 0;
          /* Stay within bounds of both RCcheck and RCtname.  */
          int RCchecklim = YYLAST - RCn + 1;
          int RCxend = RCchecklim < YYNTOKENS ? RCchecklim : YYNTOKENS;
          int RCx;

          for (RCx = RCxbegin; RCx < RCxend; ++RCx)
            if (RCcheck[RCx + RCn] == RCx && RCx != YYTERROR
                && !RCtable_value_is_error (RCtable[RCx + RCn]))
              {
                if (RCcount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    RCcount = 1;
                    RCsize = RCsize0;
                    break;
                  }
                RCarg[RCcount++] = RCtname[RCx];
                {
                  YYSIZE_T RCsize1 = RCsize + RCtnamerr (YY_NULLPTR, RCtname[RCx]);
                  if (RCsize <= RCsize1 && RCsize1 <= YYSTACK_ALLOC_MAXIMUM)
                    RCsize = RCsize1;
                  else
                    return 2;
                }
              }
        }
    }

  switch (RCcount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        RCformat = S;                       \
      break
    default: /* Avoid compiler warnings. */
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T RCsize1 = RCsize + RCstrlen (RCformat);
    if (RCsize <= RCsize1 && RCsize1 <= YYSTACK_ALLOC_MAXIMUM)
      RCsize = RCsize1;
    else
      return 2;
  }

  if (*RCmsg_alloc < RCsize)
    {
      *RCmsg_alloc = 2 * RCsize;
      if (! (RCsize <= *RCmsg_alloc
             && *RCmsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *RCmsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *RCp = *RCmsg;
    int RCi = 0;
    while ((*RCp = *RCformat) != '\0')
      if (*RCp == '%' && RCformat[1] == 's' && RCi < RCcount)
        {
          RCp += RCtnamerr (RCp, RCarg[RCi++]);
          RCformat += 2;
        }
      else
        {
          RCp++;
          RCformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
RCdestruct (const char *RCmsg, int RCtype, YYSTYPE *RCvaluep)
{
  YYUSE (RCvaluep);
  if (!RCmsg)
    RCmsg = "Deleting";
  YY_SYMBOL_PRINT (RCmsg, RCtype, RCvaluep, RClocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (RCtype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int RCchar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE RClval;
/* Number of syntax errors so far.  */
int RCnerrs;


/*----------.
| RCparse.  |
`----------*/

int
RCparse (void)
{
    int RCstate;
    /* Number of tokens to shift before error messages enabled.  */
    int RCerrstatus;

    /* The stacks and their tools:
       'RCss': related to states.
       'RCvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow RCoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    RCtype_int16 RCssa[YYINITDEPTH];
    RCtype_int16 *RCss;
    RCtype_int16 *RCssp;

    /* The semantic value stack.  */
    YYSTYPE RCvsa[YYINITDEPTH];
    YYSTYPE *RCvs;
    YYSTYPE *RCvsp;

    YYSIZE_T RCstacksize;

  int RCn;
  int RCresult;
  /* Lookahead token as an internal (translated) token number.  */
  int RCtoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE RCval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char RCmsgbuf[128];
  char *RCmsg = RCmsgbuf;
  YYSIZE_T RCmsg_alloc = sizeof RCmsgbuf;
#endif

#define YYPOPSTACK(N)   (RCvsp -= (N), RCssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int RClen = 0;

  RCssp = RCss = RCssa;
  RCvsp = RCvs = RCvsa;
  RCstacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  RCstate = 0;
  RCerrstatus = 0;
  RCnerrs = 0;
  RCchar = YYEMPTY; /* Cause a token to be read.  */
  goto RCsetstate;


/*------------------------------------------------------------.
| RCnewstate -- push a new state, which is found in RCstate.  |
`------------------------------------------------------------*/
RCnewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  RCssp++;


/*--------------------------------------------------------------------.
| RCnewstate -- set current state (the top of the stack) to RCstate.  |
`--------------------------------------------------------------------*/
RCsetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", RCstate));
  YY_ASSERT (0 <= RCstate && RCstate < YYNSTATES);
  *RCssp = (RCtype_int16) RCstate;

  if (RCss + RCstacksize - 1 <= RCssp)
#if !defined RCoverflow && !defined YYSTACK_RELOCATE
    goto RCexhaustedlab;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T RCsize = (YYSIZE_T) (RCssp - RCss + 1);

# if defined RCoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *RCvs1 = RCvs;
        RCtype_int16 *RCss1 = RCss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if RCoverflow is a macro.  */
        RCoverflow (YY_("memory exhausted"),
                    &RCss1, RCsize * sizeof (*RCssp),
                    &RCvs1, RCsize * sizeof (*RCvsp),
                    &RCstacksize);
        RCss = RCss1;
        RCvs = RCvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= RCstacksize)
        goto RCexhaustedlab;
      RCstacksize *= 2;
      if (YYMAXDEPTH < RCstacksize)
        RCstacksize = YYMAXDEPTH;

      {
        RCtype_int16 *RCss1 = RCss;
        union RCalloc *RCptr =
          (union RCalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (RCstacksize));
        if (! RCptr)
          goto RCexhaustedlab;
        YYSTACK_RELOCATE (RCss_alloc, RCss);
        YYSTACK_RELOCATE (RCvs_alloc, RCvs);
# undef YYSTACK_RELOCATE
        if (RCss1 != RCssa)
          YYSTACK_FREE (RCss1);
      }
# endif

      RCssp = RCss + RCsize - 1;
      RCvsp = RCvs + RCsize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long) RCstacksize));

      if (RCss + RCstacksize - 1 <= RCssp)
        YYABORT;
    }
#endif /* !defined RCoverflow && !defined YYSTACK_RELOCATE */

  if (RCstate == YYFINAL)
    YYACCEPT;

  goto RCbackup;


/*-----------.
| RCbackup.  |
`-----------*/
RCbackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  RCn = RCpact[RCstate];
  if (RCpact_value_is_default (RCn))
    goto RCdefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (RCchar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      RCchar = RClex ();
    }

  if (RCchar <= YYEOF)
    {
      RCchar = RCtoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      RCtoken = YYTRANSLATE (RCchar);
      YY_SYMBOL_PRINT ("Next token is", RCtoken, &RClval, &RClloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  RCn += RCtoken;
  if (RCn < 0 || YYLAST < RCn || RCcheck[RCn] != RCtoken)
    goto RCdefault;
  RCn = RCtable[RCn];
  if (RCn <= 0)
    {
      if (RCtable_value_is_error (RCn))
        goto RCerrlab;
      RCn = -RCn;
      goto RCreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (RCerrstatus)
    RCerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", RCtoken, &RClval, &RClloc);

  /* Discard the shifted token.  */
  RCchar = YYEMPTY;

  RCstate = RCn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++RCvsp = RClval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END
  goto RCnewstate;


/*-----------------------------------------------------------.
| RCdefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
RCdefault:
  RCn = RCdefact[RCstate];
  if (RCn == 0)
    goto RCerrlab;
  goto RCreduce;


/*-----------------------------.
| RCreduce -- do a reduction.  |
`-----------------------------*/
RCreduce:
  /* RCn is the number of a rule to reduce with.  */
  RClen = RCr2[RCn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  RCval = RCvsp[1-RClen];


  YY_REDUCE_PRINT (RCn);
  switch (RCn)
    {
  case 2:
#line 89 "RunConfigGrammar.y"
    {
}
#line 1350 "RunConfigGrammar.tab.c"
    break;

  case 3:
#line 92 "RunConfigGrammar.y"
    {
}
#line 1357 "RunConfigGrammar.tab.c"
    break;

  case 8:
#line 103 "RunConfigGrammar.y"
    {
  NetworkState network_state;
  double value = (RCvsp[-1].expr)->eval(NULL, network_state);
  config->setParameter((RCvsp[-3].str), value);
  free((RCvsp[-3].str));
  delete (RCvsp[-1].expr);
}
#line 1369 "RunConfigGrammar.tab.c"
    break;

  case 9:
#line 113 "RunConfigGrammar.y"
    {
  Node* node = network->getNode((RCvsp[-5].str));
  NetworkState network_state;
  double value = (RCvsp[-1].expr)->eval(NULL, network_state);
  if (!strcasecmp((RCvsp[-3].str), "istate")) {
    if (backward_istate) {
      node->setIState((bool)value);
    } else {
      if (value < 0) {
	new IStateGroup(network, node);
      } else {
	new IStateGroup(network, node, (RCvsp[-1].expr));
      }
    }
  } else if (!strcasecmp((RCvsp[-3].str), "is_internal")) {
    node->isInternal((bool)value);
  } else if (!strcasecmp((RCvsp[-3].str), "refstate")) {
    if (value < 0) {
      node->setReference(false);
    } else {
      node->setReferenceState((NodeState)value);
    }
  } else {
    throw BNException(std::string(RC_error_head() + "invalid node attribute: ") + (RCvsp[-3].str) + ", valid attributes are: istate or is_internal");
  }

  free((RCvsp[-5].str));
  free((RCvsp[-3].str));
  delete (RCvsp[-1].expr);
}
#line 1404 "RunConfigGrammar.tab.c"
    break;

  case 10:
#line 144 "RunConfigGrammar.y"
    {
  if (strcasecmp((RCvsp[-3].str), "istate")) {
    throw BNException(std::string(RC_error_head() + "invalid node group attribute: ") + (RCvsp[-3].str) + ", valid attribute is istate");
  }
  std::string error_msg;
  new IStateGroup(network, (RCvsp[-5].node_list), (RCvsp[-1].istate_expr_list), error_msg);
  if (error_msg.length() > 0) {
    throw BNException(std::string(RC_error_head() + error_msg));
  }
  
  free((RCvsp[-3].str));
}
#line 1421 "RunConfigGrammar.tab.c"
    break;

  case 11:
#line 159 "RunConfigGrammar.y"
    {
  (RCval.node_list) = (RCvsp[-1].node_list);
}
#line 1429 "RunConfigGrammar.tab.c"
    break;

  case 12:
#line 165 "RunConfigGrammar.y"
    {
  (RCval.node_list) = new std::vector<const Node*>();
  (RCval.node_list)->push_back(network->getNode((RCvsp[0].str)));
  free((RCvsp[0].str));
}
#line 1439 "RunConfigGrammar.tab.c"
    break;

  case 13:
#line 171 "RunConfigGrammar.y"
    {
  (RCval.node_list) = (RCvsp[-2].node_list);
  (RCval.node_list)->push_back(network->getNode((RCvsp[0].str)));
}
#line 1448 "RunConfigGrammar.tab.c"
    break;

  case 14:
#line 178 "RunConfigGrammar.y"
    {
  (RCval.istate_expr_list) = new std::vector<IStateGroup::ProbaIState*>();
  (RCval.istate_expr_list)->push_back((RCvsp[0].istate_expr));
}
#line 1457 "RunConfigGrammar.tab.c"
    break;

  case 15:
#line 183 "RunConfigGrammar.y"
    {
  (RCval.istate_expr_list) = (RCvsp[-2].istate_expr_list);
  (RCval.istate_expr_list)->push_back((RCvsp[0].istate_expr));
}
#line 1466 "RunConfigGrammar.tab.c"
    break;

  case 16:
#line 190 "RunConfigGrammar.y"
    {
  (RCval.istate_expr) = new IStateGroup::ProbaIState((RCvsp[-3].expr), (RCvsp[-1].expr_list));
  delete (RCvsp[-3].expr);
  for (std::vector<Expression*>::iterator it = (RCvsp[-1].expr_list)->begin(); it != (RCvsp[-1].expr_list)->end(); ++it) {
    delete *it;
  }
  delete (RCvsp[-1].expr_list);
}
#line 1479 "RunConfigGrammar.tab.c"
    break;

  case 17:
#line 201 "RunConfigGrammar.y"
    {
  (RCval.expr_list) = new std::vector<Expression*>();
  (RCval.expr_list)->push_back((RCvsp[0].expr));
}
#line 1488 "RunConfigGrammar.tab.c"
    break;

  case 18:
#line 206 "RunConfigGrammar.y"
    {
  (RCval.expr_list) = (RCvsp[-2].expr_list);
  (RCval.expr_list)->push_back((RCvsp[0].expr));
}
#line 1497 "RunConfigGrammar.tab.c"
    break;

  case 19:
#line 213 "RunConfigGrammar.y"
    {
  const Symbol* symbol = network->getSymbolTable()->getOrMakeSymbol((RCvsp[-3].str));
  free((RCvsp[-3].str));
  NetworkState dummy_state;
  network->getSymbolTable()->setSymbolValue(symbol, (RCvsp[-1].expr)->eval(NULL, dummy_state));
  delete (RCvsp[-1].expr);
}
#line 1509 "RunConfigGrammar.tab.c"
    break;

  case 20:
#line 223 "RunConfigGrammar.y"
    {
  (RCval.expr) = new ConstantExpression((RCvsp[0].l));
}
#line 1517 "RunConfigGrammar.tab.c"
    break;

  case 21:
#line 227 "RunConfigGrammar.y"
    {
  (RCval.expr) = new ConstantExpression((RCvsp[0].d));
}
#line 1525 "RunConfigGrammar.tab.c"
    break;

  case 22:
#line 231 "RunConfigGrammar.y"
    {
  (RCval.expr) = new SymbolExpression(network->getSymbolTable(), network->getSymbolTable()->getOrMakeSymbol((RCvsp[0].str)));
  free((RCvsp[0].str));
}
#line 1534 "RunConfigGrammar.tab.c"
    break;

  case 23:
#line 236 "RunConfigGrammar.y"
    {
  (RCval.expr) = new ParenthesisExpression((RCvsp[-1].expr));
}
#line 1542 "RunConfigGrammar.tab.c"
    break;

  case 24:
#line 242 "RunConfigGrammar.y"
    {
  (RCval.arg_list) = new ArgumentList();
  (RCval.arg_list)->push_back((RCvsp[0].expr));
}
#line 1551 "RunConfigGrammar.tab.c"
    break;

  case 25:
#line 247 "RunConfigGrammar.y"
    {
  (RCval.arg_list) = (RCvsp[-2].arg_list);
  (RCval.arg_list)->push_back((RCvsp[0].expr));
}
#line 1560 "RunConfigGrammar.tab.c"
    break;

  case 26:
#line 255 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1568 "RunConfigGrammar.tab.c"
    break;

  case 27:
#line 259 "RunConfigGrammar.y"
    {
  (RCval.expr) = new FuncCallExpression((RCvsp[-3].str), (RCvsp[-1].arg_list));
  free((RCvsp[-3].str));
}
#line 1577 "RunConfigGrammar.tab.c"
    break;

  case 28:
#line 264 "RunConfigGrammar.y"
    {
  (RCval.expr) = new FuncCallExpression((RCvsp[-2].str), NULL);
  free((RCvsp[-2].str));
}
#line 1586 "RunConfigGrammar.tab.c"
    break;

  case 29:
#line 272 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1594 "RunConfigGrammar.tab.c"
    break;

  case 30:
#line 276 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1602 "RunConfigGrammar.tab.c"
    break;

  case 31:
#line 280 "RunConfigGrammar.y"
    {
  (RCval.expr) = new SubExpression(new ConstantExpression(0.0), (RCvsp[0].expr));
}
#line 1610 "RunConfigGrammar.tab.c"
    break;

  case 32:
#line 284 "RunConfigGrammar.y"
    {
  (RCval.expr) = new NotLogicalExpression((RCvsp[0].expr));
}
#line 1618 "RunConfigGrammar.tab.c"
    break;

  case 33:
#line 288 "RunConfigGrammar.y"
    {
  (RCval.expr) = new NotLogicalExpression((RCvsp[0].expr));
}
#line 1626 "RunConfigGrammar.tab.c"
    break;

  case 34:
#line 294 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1634 "RunConfigGrammar.tab.c"
    break;

  case 35:
#line 298 "RunConfigGrammar.y"
    {
  (RCval.expr) = new MulExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1642 "RunConfigGrammar.tab.c"
    break;

  case 36:
#line 302 "RunConfigGrammar.y"
    {
  (RCval.expr) = new DivExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1650 "RunConfigGrammar.tab.c"
    break;

  case 37:
#line 308 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1658 "RunConfigGrammar.tab.c"
    break;

  case 38:
#line 312 "RunConfigGrammar.y"
    {
  (RCval.expr) = new AddExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1666 "RunConfigGrammar.tab.c"
    break;

  case 39:
#line 316 "RunConfigGrammar.y"
    {
  (RCval.expr) = new SubExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1674 "RunConfigGrammar.tab.c"
    break;

  case 40:
#line 322 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1682 "RunConfigGrammar.tab.c"
    break;

  case 41:
#line 326 "RunConfigGrammar.y"
    {
  (RCval.expr) = new LetterExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1690 "RunConfigGrammar.tab.c"
    break;

  case 42:
#line 330 "RunConfigGrammar.y"
    {
  (RCval.expr) = new GreaterExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1698 "RunConfigGrammar.tab.c"
    break;

  case 43:
#line 334 "RunConfigGrammar.y"
    {
  (RCval.expr) = new LetterOrEqualExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1706 "RunConfigGrammar.tab.c"
    break;

  case 44:
#line 338 "RunConfigGrammar.y"
    {
  (RCval.expr) = new GreaterOrEqualExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1714 "RunConfigGrammar.tab.c"
    break;

  case 45:
#line 344 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1722 "RunConfigGrammar.tab.c"
    break;

  case 46:
#line 348 "RunConfigGrammar.y"
    {
  (RCval.expr) = new EqualExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1730 "RunConfigGrammar.tab.c"
    break;

  case 47:
#line 352 "RunConfigGrammar.y"
    {
  (RCval.expr) = new NotEqualExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1738 "RunConfigGrammar.tab.c"
    break;

  case 48:
#line 358 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1746 "RunConfigGrammar.tab.c"
    break;

  case 49:
#line 362 "RunConfigGrammar.y"
    {
  (RCval.expr) = new AndLogicalExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1754 "RunConfigGrammar.tab.c"
    break;

  case 50:
#line 368 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1762 "RunConfigGrammar.tab.c"
    break;

  case 51:
#line 372 "RunConfigGrammar.y"
    {
  (RCval.expr) = new OrLogicalExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1770 "RunConfigGrammar.tab.c"
    break;

  case 52:
#line 378 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1778 "RunConfigGrammar.tab.c"
    break;

  case 53:
#line 382 "RunConfigGrammar.y"
    {
  (RCval.expr) = new XorLogicalExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1786 "RunConfigGrammar.tab.c"
    break;

  case 54:
#line 386 "RunConfigGrammar.y"
    {
  (RCval.expr) = new XorLogicalExpression((RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1794 "RunConfigGrammar.tab.c"
    break;

  case 55:
#line 392 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1802 "RunConfigGrammar.tab.c"
    break;

  case 56:
#line 396 "RunConfigGrammar.y"
    {
  (RCval.expr) = new CondExpression((RCvsp[-4].expr), (RCvsp[-2].expr), (RCvsp[0].expr));
}
#line 1810 "RunConfigGrammar.tab.c"
    break;

  case 57:
#line 402 "RunConfigGrammar.y"
    {
  (RCval.expr) = (RCvsp[0].expr);
}
#line 1818 "RunConfigGrammar.tab.c"
    break;


#line 1822 "RunConfigGrammar.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter RCchar, and that requires
     that RCtoken be updated with the new translation.  We take the
     approach of translating immediately before every use of RCtoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering RCchar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", RCr1[RCn], &RCval, &RCloc);

  YYPOPSTACK (RClen);
  RClen = 0;
  YY_STACK_PRINT (RCss, RCssp);

  *++RCvsp = RCval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int RClhs = RCr1[RCn] - YYNTOKENS;
    const int RCi = RCpgoto[RClhs] + *RCssp;
    RCstate = (0 <= RCi && RCi <= YYLAST && RCcheck[RCi] == *RCssp
               ? RCtable[RCi]
               : RCdefgoto[RClhs]);
  }

  goto RCnewstate;


/*--------------------------------------.
| RCerrlab -- here on detecting error.  |
`--------------------------------------*/
RCerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  RCtoken = RCchar == YYEMPTY ? YYEMPTY : YYTRANSLATE (RCchar);

  /* If not already recovering from an error, report this error.  */
  if (!RCerrstatus)
    {
      ++RCnerrs;
#if ! YYERROR_VERBOSE
      RCerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR RCsyntax_error (&RCmsg_alloc, &RCmsg, \
                                        RCssp, RCtoken)
      {
        char const *RCmsgp = YY_("syntax error");
        int RCsyntax_error_status;
        RCsyntax_error_status = YYSYNTAX_ERROR;
        if (RCsyntax_error_status == 0)
          RCmsgp = RCmsg;
        else if (RCsyntax_error_status == 1)
          {
            if (RCmsg != RCmsgbuf)
              YYSTACK_FREE (RCmsg);
            RCmsg = (char *) YYSTACK_ALLOC (RCmsg_alloc);
            if (!RCmsg)
              {
                RCmsg = RCmsgbuf;
                RCmsg_alloc = sizeof RCmsgbuf;
                RCsyntax_error_status = 2;
              }
            else
              {
                RCsyntax_error_status = YYSYNTAX_ERROR;
                RCmsgp = RCmsg;
              }
          }
        RCerror (RCmsgp);
        if (RCsyntax_error_status == 2)
          goto RCexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (RCerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (RCchar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (RCchar == YYEOF)
            YYABORT;
        }
      else
        {
          RCdestruct ("Error: discarding",
                      RCtoken, &RClval);
          RCchar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto RCerrlab1;


/*---------------------------------------------------.
| RCerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
RCerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label RCerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (RClen);
  RClen = 0;
  YY_STACK_PRINT (RCss, RCssp);
  RCstate = *RCssp;
  goto RCerrlab1;


/*-------------------------------------------------------------.
| RCerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
RCerrlab1:
  RCerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      RCn = RCpact[RCstate];
      if (!RCpact_value_is_default (RCn))
        {
          RCn += YYTERROR;
          if (0 <= RCn && RCn <= YYLAST && RCcheck[RCn] == YYTERROR)
            {
              RCn = RCtable[RCn];
              if (0 < RCn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (RCssp == RCss)
        YYABORT;


      RCdestruct ("Error: popping",
                  RCstos[RCstate], RCvsp);
      YYPOPSTACK (1);
      RCstate = *RCssp;
      YY_STACK_PRINT (RCss, RCssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++RCvsp = RClval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", RCstos[RCn], RCvsp, RClsp);

  RCstate = RCn;
  goto RCnewstate;


/*-------------------------------------.
| RCacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
RCacceptlab:
  RCresult = 0;
  goto RCreturn;


/*-----------------------------------.
| RCabortlab -- YYABORT comes here.  |
`-----------------------------------*/
RCabortlab:
  RCresult = 1;
  goto RCreturn;


#if !defined RCoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| RCexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
RCexhaustedlab:
  RCerror (YY_("memory exhausted"));
  RCresult = 2;
  /* Fall through.  */
#endif


/*-----------------------------------------------------.
| RCreturn -- parsing is finished, return the result.  |
`-----------------------------------------------------*/
RCreturn:
  if (RCchar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      RCtoken = YYTRANSLATE (RCchar);
      RCdestruct ("Cleanup: discarding lookahead",
                  RCtoken, &RClval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (RClen);
  YY_STACK_PRINT (RCss, RCssp);
  while (RCssp != RCss)
    {
      RCdestruct ("Cleanup: popping",
                  RCstos[*RCssp], RCvsp);
      YYPOPSTACK (1);
    }
#ifndef RCoverflow
  if (RCss != RCssa)
    YYSTACK_FREE (RCss);
#endif
#if YYERROR_VERBOSE
  if (RCmsg != RCmsgbuf)
    YYSTACK_FREE (RCmsg);
#endif
  return RCresult;
}
#line 407 "RunConfigGrammar.y"


void runconfig_setNetwork(Network* _network)
{
  network = _network;
}

void runconfig_setConfig(RunConfig* _config)
{
  config = _config;
}

#include "lex.RC.cc"

void RC_scan_expression(const char* str)
{
    RC_switch_to_buffer(RC_scan_string(str));
}
