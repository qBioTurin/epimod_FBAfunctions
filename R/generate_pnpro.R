#' Generate a PNPRO XML file from an arc table  — debug build
#'
#' @param arc_df     Data-frame with columns
#'                   `place, transition, direction, multiplicity, command`
#' @param pnpro_out  Path of the PNPRO file to write.
#' @export
generate_pnpro <- function(arc_df, pnpro_out) {

  ## ───────────────────────────────────────────────────────────────
  ## 0) quick overview
  ## ───────────────────────────────────────────────────────────────
  message("\n─── generate_pnpro() ──────────────────────────────────────")
  message("   distinct places     : ", length(unique(arc_df$place)))
  message("   distinct transitions: ", length(unique(arc_df$transition)))
  message("   arc rows            : ", nrow(arc_df))

  # ----------------------------------------------------------------
  # 1) unique places & transitions
  # ----------------------------------------------------------------
  places      <- unique(arc_df$place)
  transitions <- unique(arc_df$transition)

  # commands lookup
  command_df <- arc_df[!is.na(arc_df$command) & arc_df$command != "",
                       c("transition", "command")]
  command_df <- command_df[!duplicated(command_df$transition), ]
  commands   <- setNames(command_df$command, command_df$transition)

  # ----------------------------------------------------------------
  # 2) XML header + node blocks
  # ----------------------------------------------------------------
  xml_header <- '<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- This project file has been saved by the New GreatSPN Editor, v.100 --><project name="Generated Project" version="121">
  <gspn name="PetriNet" zoom="50">
    <nodes>'

  ## places (auto-layout)
  places_xml <- ""
  y <- 3.0
  for (place in places) {
    places_xml <- paste0(
      places_xml,
      sprintf('\n      <place label-x="0.0" label-y="0.0" name="%s" x="6.0" y="%s"/>',
              place, y))
    y <- y + 2.0
  }

  ## transitions (auto-layout)
  transitions_xml <- ""
  y <- 7.0;  x <- 4.0
  for (tr in transitions) {

    if (tr %in% names(commands)) {
      cmd <- gsub('"', '&quot;', commands[tr])
      delay_attr <- sprintf('delay="%s"', cmd)
    } else {
      delay_attr <- 'delay="1.0"'
    }

    message("   ↪ transition XML: ", tr, "   (", delay_attr, ")")

    transitions_xml <- paste0(
      transitions_xml,
      sprintf('\n      <transition %s delay-x="0.0" delay-y="0.0" label-x="0.0" label-y="0.0" name="%s" nservers-x="0.5" rotation="0.0" type="EXP" x="%s" y="%s"/>',
              delay_attr, tr, x, y))

    x <- x + 3.0
    if (x > 15) { x <- 4.0; y <- y + 3.0 }
  }

  # ----------------------------------------------------------------
  # 3) edges / arcs
  # ----------------------------------------------------------------
  edges_header <- '\n    </nodes>\n    <edges>'
  edges_xml    <- ""
  skipped      <- 0

  for (i in seq_len(nrow(arc_df))) {
    row <- arc_df[i, ]
    mult_attr <- if (!is.na(row$multiplicity) && row$multiplicity > 1)
                   sprintf(' mult="%d"', row$multiplicity) else ""

    if (is.na(row$direction) ||
        !row$direction %in% c("INPUT", "OUTPUT")) {
      message("   ⚠ skipping arc row ", i,
              " – invalid direction: ", row$direction)
      skipped <- skipped + 1
      next
    }

    edges_xml <- if (row$direction == "INPUT") {
      paste0(edges_xml,
             sprintf('\n      <arc head="%s" kind="INPUT"%s tail="%s"/>',
                     row$transition, mult_attr, row$place))
    } else {
      paste0(edges_xml,
             sprintf('\n      <arc head="%s" kind="OUTPUT"%s tail="%s"/>',
                     row$place,      mult_attr, row$transition))
    }
  }

  message("✔ arcs rendered      : ", nrow(arc_df) - skipped)
  if (skipped) message("⚠ arcs skipped       : ", skipped)

  # ----------------------------------------------------------------
  # 4) footer
  # ----------------------------------------------------------------
  xml_footer <- '\n    </edges>
  </gspn>
  <measures gspn-name="PetriNet" name="Measures" simplified-UI="false">
    <assignments/>
    <greatspn/>
    <formulas>
      <formula comment="Basic statistics of the toolchain execution." language="STAT"/>
      <formula comment="All the basic Petri net measures" language="ALL"/>
    </formulas>
  </measures>
</project>'

  # assemble + write
  xml_content <- paste0(
    xml_header, places_xml, transitions_xml, edges_header, edges_xml,
    xml_footer)

  dir.create(dirname(pnpro_out), recursive = TRUE, showWarnings = FALSE)
  writeLines(xml_content, pnpro_out)
  message("✓ PNPRO written to: ", pnpro_out)

  invisible(pnpro_out)
}

