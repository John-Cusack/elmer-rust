/* Build-time helper: expand the MATC `$` prelude of SOLVER.KEYWORDS.
 *
 * Elmer's CheckKeyword reads SOLVER.KEYWORDS through ReadAndTrim with
 * evaluation on, so `$ for(...)...` lines execute and their printed rows
 * (e.g. `Solver:String: 'Tensor Field 3 Component 7'`) become keyword
 * entries. This probe reproduces that: for every line of the snapshot
 * starting with `$`, evaluate the rest of the line with the vendored
 * libmatc (same init as Load.c:mtc_init_once) and print the result.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *mtc_domath(char *);
void mtc_init(FILE *, FILE *, FILE *);

int main(int argc, char **argv) {
    char line[163840];
    FILE *f;

    if (argc != 2) return 2;
    f = fopen(argv[1], "r");
    if (!f) return 1;

    mtc_init(NULL, stdout, stderr);
    mtc_domath("format( 12,\"rowform\")");

    while (fgets(line, sizeof(line), f)) {
        char *cmd, *out;
        size_t n;
        if (line[0] != '$') continue;
        /* Strip trailing newline (mtc_domath output carries its own). */
        n = strlen(line);
        while (n > 1 && (line[n - 1] == '\n' || line[n - 1] == '\r')) line[--n] = '\0';
        cmd = line + 1;
        out = mtc_domath(cmd);
        if (out) fputs(out, stdout);
    }
    fclose(f);
    return 0;
}
