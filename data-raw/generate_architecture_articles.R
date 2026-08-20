# Generate vignettes/articles/architecture-*.Rmd from docs/*.md - maintainer-run, NOT part of the
# package build (mirrors data-raw/benchmark_performance.R's own "not part of automated
# testthat/build" convention). docs/*.md is the single source of truth; the articles are a
# pkgdown-facing derived copy (YAML frontmatter + cross-module links resolved to
# architecture-*.html instead of bare *.md filenames), so this replaces what used to be two
# independently hand-edited, silently-diverging copies of the same content.
#
# Run manually whenever docs/*.md changes, before `pkgdown::build_site()`:
#   Rscript data-raw/generate_architecture_articles.R
# then review the diff (`git diff vignettes/articles/`) before committing - this OVERWRITES every
# architecture-*.Rmd file unconditionally.

docs_dir <- "docs"
articles_dir <- file.path("vignettes", "articles")

doc_files <- list.files(docs_dir, pattern = "^[0-9]{2}_.*\\.md$")
if (length(doc_files) == 0) {
  stop("No docs/NN_*.md files found - run this from the package root (soilSIM/).")
}

slug_for <- function(filename) {
  gsub("_", "-", sub("\\.md$", "", sub("^[0-9]+_", "", filename)))
}

title_for <- function(filename) {
  first_line <- readLines(file.path(docs_dir, filename), n = 1)
  if (!startsWith(first_line, "# ")) {
    stop("Expected docs/", filename, " to start with a top-level '# Title' heading, found: ", first_line)
  }
  sub("^# ", "", first_line)
}

# Build the full filename -> (title, slug) lookup table FIRST (two-pass), since every doc's body
# can reference every other doc, not just itself.
doc_meta <- lapply(doc_files, function(fn) list(title = title_for(fn), slug = slug_for(fn)))
names(doc_meta) <- doc_files

cat("Generating", length(doc_files), "architecture articles from", docs_dir, "...\n")

for (fn in doc_files) {
  lines <- readLines(file.path(docs_dir, fn), warn = FALSE)

  # Drop the leading '# Title' line (becomes the YAML title instead) and one following blank
  # line, if present.
  body_lines <- lines[-1]
  if (length(body_lines) > 0 && body_lines[1] == "") {
    body_lines <- body_lines[-1]
  }
  body <- paste(body_lines, collapse = "\n")

  # Resolve every cross-reference to another doc (markdown-link or bare-backtick form) into a
  # pkgdown-style [Title](architecture-slug.html) link, using the TARGET file's own real title as
  # the link text (matching this codebase's existing convention - discards whatever text was
  # originally inside the source [...] link, since it was always just the target's own title or
  # filename anyway).
  for (target_fn in doc_files) {
    target <- doc_meta[[target_fn]]
    replacement <- paste0("[", target$title, "](architecture-", target$slug, ".html)")

    # [any link text](NN_target_file.md) -> replacement
    link_pattern <- paste0("\\[[^]]*\\]\\(", target_fn, "\\)")
    body <- gsub(link_pattern, replacement, body, fixed = FALSE)

    # `NN_target_file.md` (bare backtick-code reference, no markdown link) -> replacement
    backtick_pattern <- paste0("`", target_fn, "`")
    body <- gsub(backtick_pattern, replacement, body, fixed = TRUE)
  }

  meta <- doc_meta[[fn]]
  out <- paste0(
    "---\n",
    "title: \"", meta$title, "\"\n",
    "---\n\n",
    body
  )

  out_path <- file.path(articles_dir, paste0("architecture-", meta$slug, ".Rmd"))
  writeLines(out, out_path, useBytes = TRUE)
  cat("  wrote", out_path, "\n")
}

cat("Done. Review with `git diff vignettes/articles/` before committing.\n")
