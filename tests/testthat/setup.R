# Factor level ordering depends on the collation locale. Pin it so the golden
# fixtures (generated under C collation) compare deterministically across machines.
Sys.setlocale("LC_COLLATE", "C")
