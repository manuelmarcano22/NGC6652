/*
 * This file is part of the ESO Recipe Execution Tool
 * Copyright (C) 2001-2017 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>

#include "ltdl.h"

#include <cpl.h>
#include <cpl_version.h>        /* should be fixed */

#include "er_macros.h"
#include "er_help.h"
#include "er_stringarray.h"
#include "er_params.h"
#include "er_paramutils.h"
#include "er_plugin.h"
#include "er_fileutils.h"
#include "er_stringutils.h"
#ifdef ENABLE_PYTHON_RECIPES
#include "er_python.h"
#endif


/**
 * @defgroup er_main Main function for EsoRex
 *
 * Module containing the top-level function for the EsoRex application.
 *
 */


/**@{*/

/**
 * @brief   Fills a list with EsoRex-specific parameters.
 *
 * @param param_list  List of parameters
 */

static void
er_init_parameters(cpl_parameterlist *param_list)
{
    cpl_parameter *p;


    /*
     * Processing different EsoRex options
     */

    /*
     *  help
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".help",
                                 CPL_TYPE_BOOL,
                                 "Display this help and exit. "
                                 "If a recipe name is also given, then "
                                 "help will be given for it as well.",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "help");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_CFG);
    cpl_parameterlist_append (param_list, p);

    /*
     *  version
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".version",
                                 CPL_TYPE_BOOL,
                                 "Display version information and exit.",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "version");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_CFG);
    cpl_parameterlist_append (param_list, p);

    /*
     *      check-sof-exist
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".check-sof-exist",
                                 CPL_TYPE_BOOL,
                                 "When TRUE, all the input files must exist  "
                                 "and be readable before calling the recipe. ",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "check-sof-exist");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_CHECK_SOF_EXIST");
    cpl_parameterlist_append (param_list, p);

    /*
     *  config
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".config",
                                 CPL_TYPE_STRING,
                                 "Configuration file to be used for EsoRex.",
                                 PACKAGE_RESOURCE, NULL);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "config");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_CONFIG");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_CFG);
    cpl_parameterlist_append (param_list, p);

    /*
     *  create-config
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".create-config",
                                 /*
                             CPL_TYPE_BOOL,
                                  */
                                 CPL_TYPE_STRING,
                                 "Creates a configuration file for Esorex. "
                                 "If set to TRUE a config file 'esorex.rc' is created "
                                 "in the '.esorex' directory in $HOME of the user."
                                 "If a filename is specified, a config file will be "
                                 "created accordingly. "
                                 "If a recipe is specified in the command line, then the "
                                 "configuration file will "
                                 "be created for the recipe instead (called "
                                 "'recipename.rc')"
                                 " Note that an existing file will be "
                                 "overwritten, but a backup file will be "
                                 "copied to 'filename.rc.bak' in "
                                 "the same directory. If the filename ends with "
                                 "extension .json then a machine-readable "
                                 "JSON format will be used",
                                 PACKAGE_RESOURCE, "FALSE");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "create-config");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_CREATE_CONFIG");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_CFG);
    cpl_parameterlist_append (param_list, p);

    /*
     *  link-dir
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".link-dir",
                                 CPL_TYPE_STRING,
                                 "The directory in which a symbolic link "
                                 "to each of the product files should be "
                                 "written. The "
                                 "enable/disable switch to control "
                                 "whether the link is actually made is "
                                 "the '--suppress-link' option.",
                                 PACKAGE_RESOURCE, "/tmp");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "link-dir");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_LINK_DIR");
    cpl_parameterlist_append (param_list, p);

    /*
     *      log-dir
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".log-dir",
                                 CPL_TYPE_STRING,
                                 "Directory where to place the logfile.",
                                 PACKAGE_RESOURCE, ".");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "log-dir");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_LOG_DIR");
    cpl_parameterlist_append (param_list, p);

    /*
     *      log-file
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".log-file",
                                 CPL_TYPE_STRING,
                                 "Filename of logfile.",
                                 PACKAGE_RESOURCE, PACKAGE ".log");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "log-file");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_LOG_FILE");
    cpl_parameterlist_append (param_list, p);

    /*
     *      log-level
     */

    p = cpl_parameter_new_enum (PACKAGE_RESOURCE ".log-level",
                                CPL_TYPE_STRING,
                                "Controls the severity level of messages "
                                "that will be printed to the logfile.",
                                PACKAGE_RESOURCE,
                                "info", 5, "debug", "info", "warning",
                                "error", "off");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "log-level");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_LOG_LEVEL");
    cpl_parameterlist_append (param_list, p);


    /*
     * no-datamd5 and no-checksum
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".no-datamd5",
                                 CPL_TYPE_BOOL,
                                 "Disables the computation of the MD5 data hash "
                                 "for FITS product files.", PACKAGE_RESOURCE,
                                 FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "no-datamd5");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_NO_DATAMD5");
    cpl_parameterlist_append (param_list, p);


    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".no-checksum",
                                 CPL_TYPE_BOOL,
                                 "Disables the computation of the standard "
                                 "FITS product checksums.", PACKAGE_RESOURCE,
                                 FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "no-checksum");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_NO_CHECKSUM");
    cpl_parameterlist_append (param_list, p);


    /*
     *  man-page
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".man-page",
                                 CPL_TYPE_BOOL,
                                 "Display a manual page for the specified "
                                 "recipe, and then exit. Note that this option "
                                 "only applies to recipes, and that it does "
                                 "nothing for " PACKAGE
                                 " by itself. See also "
                                 "the '--help' option.", PACKAGE_RESOURCE,
                                 FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "man-page");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_MAN_PAGE");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_CFG);
    cpl_parameterlist_append (param_list, p);

    /*
     *  memcheck
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".mem-check",
                                 CPL_TYPE_BOOL,
                                 "Report on memory status at completion "
                                 "of recipe execution.",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "mem-check");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_MEM_CHECK");
    cpl_parameterlist_append (param_list, p);

    /*
     *      msg-level
     */

    p = cpl_parameter_new_enum (PACKAGE_RESOURCE ".msg-level",
                                CPL_TYPE_STRING,
                                "Controls the severity level of messages "
                                "that will be printed to the terminal.",
                                PACKAGE_RESOURCE,
                                "info", 5, "debug", "info", "warning",
                                "error", "off");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "msg-level");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_MSG_LEVEL");
    cpl_parameterlist_append (param_list, p);

    /*
     *  output-dir
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".output-dir",
                                 CPL_TYPE_STRING,
                                 "The directory where the product "
                                 "files should be finally moved to "
                                 "(all products are first created in the "
                                 "current dir).",
                                 PACKAGE_RESOURCE, ".");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "output-dir");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_OUTPUT_DIR");
    cpl_parameterlist_append (param_list, p);

    /*
     *  output-prefix
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".output-prefix",
                                 CPL_TYPE_STRING,
                                 "Prefix applied to any output file. "
                                 "For example, specifying 'pre' would translate "
                                 "'filename.fits' to 'pre_0000.fits'. See also "
                                 "the '--suppress-prefix' option.",
                                 PACKAGE_RESOURCE, "out");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "output-prefix");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_OUTPUT_PREFIX");
    cpl_parameterlist_append (param_list, p);

    /*
     *      output-readonly  (replaces in some ways the old "output-overwrite")
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".output-readonly",
                                 CPL_TYPE_BOOL,
                                 "When TRUE, any existing product files "
                                 "in the specified output directory will be "
                                 "set to read-only, for user, group and other. "
                                 "If FALSE, then EsoRex will use the default "
                                 "permissions for that account/directory. "
                                 "destroy any pre-existing files. "
                                 "This option exists for the Paranal operations "
                                 "environment. This option can additionally be "
                                 "used to prevent EsoRex from overwriting "
                                 "pre-existing files.",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "output-readonly");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_OUTPUT_READONLY");
    cpl_parameterlist_append (param_list, p);

    /*
     *  paf-config
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".paf-config",
                                 CPL_TYPE_STRING,
                                 "Configuration file for creation of PAF files.",
                                 PACKAGE_RESOURCE, "");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "paf-config");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_PAF_CONFIG");
    cpl_parameterlist_append (param_list, p);

    /*
     *  params
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".params",
                                 CPL_TYPE_BOOL,
                                 "List the input parameters and their current "
                                 "settings (whether from the command line or a "
                                 "configuration file) for the " PACKAGE
                                 " application. Parameters are "
                                 "labelled using the parameter's alias. "
                                 "If a recipe is also specified, then the "
                                 "list of its parameters will also be generated "
                                 "in the same way.",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "params");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_PARAMS");
    cpl_parameterlist_append (param_list, p);

    /*
     *  products-sof
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".products-sof",
                                 CPL_TYPE_STRING,
                                 "Output file which contains the FITS files "
                                 "created by the recipe. If the filename ends with "
                                 "extension .json then a machine-readable JSON format"
                                 "will be used",
                                 PACKAGE_RESOURCE, "");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "products-sof");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_PRODUCTS_SOF");
    cpl_parameterlist_append (param_list, p);

    /*
     *  recipes
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".recipes",
                                 CPL_TYPE_BOOL,
                                 "Display a list of all available recipes "
                                 "(that are available in the directory tree "
                                 "specified with '--recipe-dir').",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "recipes");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_RECIPES");
    cpl_parameterlist_append (param_list, p);

    /*
     *  recipe-config
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".recipe-config",
                                 CPL_TYPE_STRING,
                                 "Configuration file for any selected recipe.",
                                 PACKAGE_RESOURCE, NULL);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "recipe-config");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_RECIPE_CONFIG");
    cpl_parameterlist_append (param_list, p);

    /*
     *  recipe-dir
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".recipe-dir",
                                 CPL_TYPE_STRING,
                                 "Directory containing recipe libraries. Note "
                                 "that " PACKAGE
                                 " will recursively search not "
                                 "only the specified directory, but all "
                                 "sub-directories below it as well. "
                                 "Multiple directory heads may be "
                                 "specified, by separating the "
                                 "starting paths with colons (:).",
                                 PACKAGE_RESOURCE, ".");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "recipe-dir");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_PLUGIN_DIR");
    cpl_parameterlist_append (param_list, p);

    /*
     *  show hidden
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".show-hidden",
                                 CPL_TYPE_BOOL,
                                 "When TRUE recipe parameters which are "
                                 "declared as hidden by the recipe are shown "
                                 "in the output of help and man page options.",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "show-hidden");
    cpl_parameterlist_append (param_list, p);

    /*
     *      suppress_link
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".suppress-link",
                                 CPL_TYPE_BOOL,
                                 "When TRUE, no symbolic link is created "
                                 "to the output product. However, "
                                 "if FALSE, then a symbolic link is created "
                                 "in the directory specified with the "
                                 "option '--link-dir' for each product "
                                 "that is created by the recipe.",
                                 PACKAGE_RESOURCE, TRUE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "suppress-link");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_SUPPRESS_LINK");
    cpl_parameterlist_append (param_list, p);

    /*
     *      suppress_prefix
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".suppress-prefix",
                                 CPL_TYPE_BOOL,
                                 "When TRUE, the original name of the "
                                 "output product, as produced by the "
                                 "recipe, is maintained. "
                                 "If FALSE, then the name of the output "
                                 "file is changed to the \"prefix_number\" "
                                 "format. The prefix can be altered using the "
                                 "'--output-prefix' option.",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "suppress-prefix");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_SUPPRESS_PREFIX");
    cpl_parameterlist_append (param_list, p);

    /*
     *  time
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".time",
                                 CPL_TYPE_BOOL,
                                 "Measure and show the recipe's execution time.",
                                 PACKAGE_RESOURCE, FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "time");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_TIME");
    cpl_parameterlist_append (param_list, p);

    /*
     *      unload-plugin
     */

    p = cpl_parameter_new_value (PACKAGE_RESOURCE ".unload-plugin",
                                 CPL_TYPE_BOOL,
                                 "When TRUE, the plugin is unloaded "
                                 "after execution. "
                                 "If FALSE, the plugin is not unloaded "
                                 "after processing, so that a software "
                                 "like, e.g. valgrind, can be used "
                                 "for debugging the executed recipe. ",
                                 PACKAGE_RESOURCE, TRUE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "unload-plugin");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_ENV,
                             PACKAGE_ENV "_UNLOAD_PLUGIN");
    cpl_parameterlist_append (param_list, p);


    /*
     * The parameter list is now filled. The remaining step is to set
     * the tag to indicate that the sources of the current values for
     * all parameters in this list are the EsoRex defaults.
     */

    p = cpl_parameterlist_get_first (param_list);
    while (p != NULL)
    {
        char *cptr, def_str[] = PACKAGE " default";

        cptr = def_str;             /* save source of this param */
        er_manage_sources(1,cpl_parameter_get_name(p),&cptr);

        /* Get the next parameter in the list */
        p = cpl_parameterlist_get_next (param_list);
    }                            /* End of loop through all parameters in the list */

}                               /* End of er_init_parameters() */


/**
 * @brief
 *   Write application header to the standard output.
 */

static void
er_print_header(void)
{
    fprintf(stdout, "\n     ***** %s, version %s  *****\n\n",
            PACKAGE_NAME, PACKAGE_VERSION);
    fflush(stdout);

    return;
}


inline static int
_esorex_config_value_get(FILE *fp, const char *key, char *value, size_t sz)
{

    int has_key = 0;

    char buffer[MAXSTRLENCONF];

    while (fgets(buffer, MAXSTRLENCONF, fp) != NULL) {

        char *s = buffer;

        while (isspace(*s)) {
            ++s;
        }

        if ((s[0] != '\0') && (s[0] != '#') &&
                (strncmp(s, key, strlen(key)) == 0)) {

            s = strchr(s, '=');

            if (s != NULL) {

                char v[MAXSTRLENCONF];

                int nitems = sscanf(s + 1, "%s", v);

                if (nitems == 1) {
                    has_key = 1;
                    strncpy(value, v, sz);
                    break;
                }
                else {
                    has_key = -1;
                    break;
                }

            }
            else {
                has_key = -1;
                break;
            }

        }

    }

    return has_key;
}


inline static void
_esorex_config_files_set(int argc, char *argv[],
                         char **syscfg, char **cfg)
{

    const char *env_config = getenv("ESOREX_CONFIG");

    int iarg;

    const char *_home   = getenv("HOME");
    const char *_cfg = NULL;

    char *_path   = NULL;

    struct stat sb;


    /*
     *  File path buffer size. Add an extra character for the directory
     *  separator, and another one for the string terminator.
     */

    size_t sz  = strlen(GLOBAL_RC_NAME) + 2;
    size_t _sz = 0;


    if (_home == NULL) {

        struct passwd *pw = getpwuid(getuid());
        _home = pw->pw_dir;

    }

    _sz = strlen(_home) + strlen(GLOBAL_RC_DIR) + 1;
    sz += (_sz > strlen(SYSCONFDIR)) ? _sz : strlen(SYSCONFDIR);

    _path = static_cast<char*>(malloc(sz));


    /*
     * Try to find a configuration file in the user's home directory.
     * If it does not exist, try the system wide as a fallback.
     */

    snprintf(_path, sz, "%s/%s/%s", _home, GLOBAL_RC_DIR, GLOBAL_RC_NAME);

    if ((stat(_path, &sb) == 0) && S_ISREG(sb.st_mode)) {
        _cfg = _path;
    }
    else {
        snprintf(_path, sz, "%s/%s", SYSCONFDIR, GLOBAL_RC_NAME);

        if ((stat(_path, &sb) == 0) && S_ISREG(sb.st_mode)) {
            _cfg = _path;
        }
    }

    if (_cfg != NULL) {
        *syscfg = strdup(_cfg);
    }

    free(_path);


    /*
     * Check if the user has explicitly specified a configuration file,
     * either setting an environment variable or by passing a command
     * line option.
     */

    _cfg = NULL;

    if (env_config != NULL) {
        _cfg = env_config;
    }

    for (iarg = 1; iarg < argc; ++iarg) {

        if ((strlen(argv[iarg]) > 2) &&
                (strncmp(argv[iarg] + 2, "config", 6) == 0)) {

            const char *s = strchr(argv[iarg] + 2, '=');

            if ((s != NULL) && (*(s + 1) != '\0')) {
                _cfg = ++s;
            }

        }

    }

    if (_cfg != NULL) {
        *cfg = strdup(_cfg);
    }

    return;

}


// Does an advance parsing of the command line, environment, and the
// configuration file, to determine whether CPL has to be initialized
// for memory debugging. This needs to be known before any CPL call
// happens!

// FIXME: This is a hack to workaround the fact that the code does
//        not do user input collection/verification before the core services
//        are started.

inline static int
_esorex_memory_mode_get(int argc, char *argv[],
                        const char *syscfg, const char *cfg)
{
    int iarg;
    int set_memory_mode = 0;

    const char *env_memcheck = getenv("ESOREX_MEM_CHECK");


    if ((syscfg != NULL) && (syscfg[0] != '\0')) {

        const char *key = PACKAGE_RESOURCE ".mem-check";

        char value[MAXSTRLENCONF];

        FILE *fp = fopen(syscfg, "r");

        if (fp != NULL) {

            if (_esorex_config_value_get(fp, key, value, MAXSTRLENCONF) == 1) {

                if ((strncmp(value, "TRUE", 4) == 0) ||
                        strncmp(value, "true", 4) == 0) {
                    set_memory_mode = 1;
                }
            }

            fclose(fp);

        }

    }

    if ((cfg != NULL) && (cfg[0] != '\0')) {

        const char *key = PACKAGE_RESOURCE ".mem-check";

        char value[MAXSTRLENCONF];

        FILE *fp = fopen(cfg, "r");

        if (fp != NULL) {

            if (_esorex_config_value_get(fp, key, value, MAXSTRLENCONF) == 1) {

                if ((strncmp(value, "TRUE", 4) == 0) ||
                        strncmp(value, "true", 4) == 0) {
                    set_memory_mode = 1;
                }
            }

            fclose(fp);
        }

    }

    if ((env_memcheck != NULL) &&
            ((strncmp(env_memcheck, "TRUE", 4) == 0) ||
                    (strncmp(env_memcheck, "true", 4) == 0))) {
        set_memory_mode = 1;
    }


    for (iarg = 1; iarg < argc; ++iarg) {

        if ((strlen(argv[iarg]) > 2) &&
                (strncmp(argv[iarg] + 2, "mem-check", 9) == 0)) {

            if (strstr(argv[iarg], "TRUE") || strstr(argv[iarg], "true")) {

                const char *memory_mode = getenv("CPL_MEMORY_MODE");

                if ((memory_mode == NULL) || (strncmp(memory_mode, "0", 1) == 0)) {
                    set_memory_mode = 1;
                }

            }
            else if (strstr(argv[iarg], "FALSE") || strstr(argv[iarg], "false")) {
                set_memory_mode = 0;
            }

        }

    }

    return set_memory_mode;

}


int main(int argc, char *argv[])
{
    const char *fn = "EsoRex";

    char  plugin_name[FILEMAX+2];
    char  *conf_file_global=NULL, *conf_file_local=NULL;

    const char *cdescr;

    cpl_parameter *p = NULL;
    cpl_parameterlist *caller_parameter_list = NULL;
    cpl_msg_severity msg_level = CPL_MSG_ERROR;

    er_stringarray_t *set_of_frames_filenames = NULL;

    int e_code = 0;
    int e_code2 = 0;
    int flag_mem_check = 1;

    unsigned  init_flag = CPL_INIT_DEFAULT;


    /* check that we're not somebody else... */

    if (getuid() != geteuid())
    {
        fprintf(stderr, "effective user ID not the same as real user ID => no EsoRex...\n");
        exit(EXIT_FAILURE);
    }

    /* print version of Esorex */

    er_print_header ();


    _esorex_config_files_set(argc, argv, &conf_file_global, &conf_file_local);

    if (_esorex_memory_mode_get(argc, argv,
                                conf_file_global, conf_file_local) != 0) {
        setenv("CPL_MEMORY_MODE", "1", 0);
    }


    /* Initialization and setup */

    cpl_init (init_flag);                    /* also does cpl_msg_init() stuff */

    cpl_msg_set_domain (PACKAGE);
    cpl_msg_set_time_off ();
    cpl_msg_set_threadid_off ();
    cpl_msg_set_domain_on ();
    cpl_msg_set_component_off ();



    /* Main processing */

    if (lt_dlinit () != 0)
    {
        cpl_msg_error (fn, "Unable to initialise ltdl; aborting program");
        cpl_end ();                  /* stop subsystems of CPL */
        e_code = -1;
        goto clean_up;
    }


    /* Create list of parameters for EsoRex */

    caller_parameter_list = cpl_parameterlist_new ();
    set_of_frames_filenames = er_stringarray_new ();

    er_init_parameters (caller_parameter_list);


    /* Process caller configuration information */

    plugin_name[0] = '\0';
    e_code = params_process_configuration (caller_parameter_list,
                                           conf_file_global, conf_file_local,
                                           argc, argv,
                                           plugin_name, set_of_frames_filenames);


    /* Process Plugin */

    if (e_code == 0 && plugin_name[0] != '\0')
    {
        // FIXME: Error handling is non-local here and hard to trace.
        //        This has to be improved!
        e_code = plugin_process_plugin (caller_parameter_list,
                                        plugin_name, set_of_frames_filenames,
                                        argc, argv);
    }
    else if(e_code !=-99999)  /* entering just Esorex yields the library versions */
    {
        cdescr = cpl_get_description(CPL_DESCRIPTION_DEFAULT);
        (void) printf("\nLibraries used: %s\n\n",cdescr);

        // Shutdown ltdl as this is not done by cleanup
        lt_dlexit();

        goto clean_up;
    }

    msg_level = cpl_msg_get_level();    /* test, if we are in debug mode */
    if (msg_level == CPL_MSG_DEBUG)
    {                       /* if debug level set, */
        e_code2 = cpl_error_get_code ();        /* check the CPL_error subsystem */
        if (e_code2 != CPL_ERROR_NONE)
        {
            cpl_msg_error (fn, "CPL-error '%s' was set in '%s'",
                           cpl_error_get_message (), cpl_error_get_where ());
            if (e_code == 0)
            {
                cpl_msg_error (fn, "recipe %s returned status: %d\n",plugin_name,e_code);
            }
        }
    }

    p = cpl_parameterlist_find (caller_parameter_list,
                                PACKAGE_RESOURCE ".unload-plugin");

    if (cpl_parameter_get_bool(p) != 0) {
        if (lt_dlexit () != 0) {
            cpl_msg_error (fn, "Unable to deinitialize ltdl");
        }
    }
    else {
        cpl_msg_warning(fn, "unloading of plugins is inhibited by program "
                        "option 'unload-plugin'!");
    }

    p = cpl_parameterlist_find (caller_parameter_list, PACKAGE_RESOURCE ".mem-check");
    flag_mem_check = cpl_parameter_get_bool (p);


    /* cleanup allocated structures and memory */

    clean_up:
    if (caller_parameter_list) {
        cpl_parameterlist_delete (caller_parameter_list);
    }

    if (set_of_frames_filenames) {
        er_stringarray_delete (set_of_frames_filenames);
    }

    if (conf_file_local != NULL)
    {
        free(conf_file_local);
    }
    if (conf_file_global != NULL)
    {
        free(conf_file_global);
    }

    (void) er_help_free();
    (void) er_manage_sources(3, "", NULL);  /* free buffer of sources */

#ifdef ENABLE_PYTHON_RECIPES
    /* Free any buffers allocated by er_python_* functions. */
    er_python_cleanup();
#endif

    if (flag_mem_check != 0)    /* see if we should check the memory */
    {
        cpl_memory_dump();
    }

    cpl_end ();         /* stop subsystems of CPL */

    // Reset error code -99999 used for non-executing options
    // like --version

    if (e_code == -99999)
    {
        e_code = 0;
    }

    return (e_code);
}
/**@}*/
/* End of file */
