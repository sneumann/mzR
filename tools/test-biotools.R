source("tools/biotools.R")
roclet <- biotools_roclet()
ns <- c(b = "biotoolsSchema")
schema <- biotools_read_xml("tools/biotools/biotools.xsd")

expect_error <- function(expr, pattern) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(error, "error"),
              grepl(pattern, conditionMessage(error), fixed = TRUE))
}
parse_blocks <- function(lines) {
    roxygen2::parse_text(c(paste0("#' ", lines), "NULL"), env = NULL)
}
process <- function(blocks, path = ".") {
    roxygen2::roclet_process(roclet, blocks, env = NULL, base_path = path)
}
texts <- function(doc, xpath) {
    xml2::xml_text(xml2::xml_find_all(doc, xpath, ns))
}

blocks <- roxygen2::parse_package(".", env = NULL)
doc <- process(blocks)
stopifnot(isTRUE(xml2::xml_validate(doc, schema)))
stopifnot(length(xml2::xml_find_all(doc, "//b:function", ns)) == 4L)
stopifnot(identical(texts(doc, "//b:name")[1], "mzR"))
description <- read.dcf("DESCRIPTION")[1, ]
stopifnot(identical(texts(doc, "//b:version"), description[["Version"]]))
stopifnot(identical(texts(doc, "//b:description"),
                    gsub("[[:space:]]+", " ", description[["Description"]])))
stopifnot(identical(texts(doc, "//b:topic/b:term"), c("Proteomics", "Metabolomics")))
stopifnot(identical(texts(doc, "//b:function[1]/b:input/b:format/b:term"),
                    c("mzML", "mzXML", "NetCDF")))
stopifnot(identical(texts(doc, "//b:function[2]/b:input/b:format/b:term"), "mzIdentML"))
stopifnot(identical(texts(doc, "//b:function[3]/b:output/b:format/b:term"),
                    c("mzML", "mzXML")))
stopifnot(length(xml2::xml_find_all(doc, "//b:publication", ns)) == 4L)
stopifnot(length(xml2::xml_find_all(doc, "//b:credit", ns)) == 7L)
stopifnot("https://orcid.org/0000-0002-6977-7147" %in% texts(doc, "//b:orcidid"))
stopifnot("Primary contact" %in% texts(doc, "//b:typeRole"))
stopifnot(any(grepl("Access spectra, chromatograms", texts(doc, "//b:note"), fixed = TRUE)))

local({
    path <- tempfile("biotools-test-")
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE))
    dir.create(file.path(path, "R"))
    dir.create(file.path(path, "tools", "biotools"), recursive = TRUE)
    file.copy("tools/biotools/biotools.xsd", file.path(path, "tools", "biotools"))
    file.copy("tools/biotools.R", file.path(path, "tools"))
    metadata <- matrix(c(
        "example", "A test package", "1.0", "Artistic-2.0", "UTF-8",
        "Read & write <spectra> with an annotated R package.",
        "https://example.org/", description[["Roxygen"]]),
        nrow = 1, dimnames = list(NULL, c(
            "Package", "Title", "Version", "License", "Encoding",
            "Description", "URL", "Roxygen")))
    write.dcf(metadata, file.path(path, "DESCRIPTION"))
    annotations <- c(
        "Read & write <spectra>", "", "Existing documentation is extracted.",
        "@edam_topic Proteomics (topic_0121)",
        "@edam_topic Proteomics (topic_0121)",
        "@edam_operation Data parsing (operation_1812)",
        "@edam_operation Data retrieval (operation_2422)",
        "@edam_input_data Mass spectrometry data (data_2536)",
        "@edam_input_format mzML (format_3244)",
        "@edam_input_format mzML (format_3244)",
        "@edam_input_data Mass spectrometry data (data_2536)",
        "@edam_input_format mzXML (format_3654)",
        "@edam_output_data Mass spectrometry data (data_2536)",
        "@biotools_tool_type Library",
        "@biotools_documentation User manual | https://example.org/manual",
        "@noRd")
    source <- c(paste0("#' ", annotations), "NULL")
    writeLines(source, file.path(path, "R", "example.R"))
    result <- process(parse_blocks(annotations), path)
    stopifnot(isTRUE(xml2::xml_validate(result, schema)))
    stopifnot(length(xml2::xml_find_all(result, "//b:topic", ns)) == 1L)
    stopifnot(length(xml2::xml_find_all(result, "//b:operation", ns)) == 2L)
    stopifnot(identical(texts(result, "//b:input[1]/b:format/b:term"), "mzML"))
    stopifnot(identical(texts(result, "//b:input[2]/b:format/b:term"), "mzXML"))
    stopifnot(identical(texts(result, "//b:note"),
        "Read & write <spectra> Existing documentation is extracted."))
    output <- roxygen2::roclet_output(roclet, result, path)
    stopifnot(isTRUE(xml2::xml_validate(biotools_read_xml(output), schema)))
    before <- readBin(output, "raw", n = file.info(output)$size)
    bad <- biotools_read_xml(as.character(result))
    tool_type <- xml2::xml_find_first(bad, "//b:toolType", ns)
    xml2::xml_text(tool_type) <- "Not a tool type"
    expect_error(roxygen2::roclet_output(roclet, bad, path), "Invalid bio.tools XML")
    stopifnot(identical(before, readBin(output, "raw", n = file.info(output)$size)))

    # Exercise DESCRIPTION's actual roclet expression through normal roxygenise.
    old <- setwd(path)
    on.exit(setwd(old), add = TRUE)
    roxygen2::roxygenise(".", load_code = function(path) new.env())
    generated <- biotools_read_xml(output)
    stopifnot(isTRUE(xml2::xml_validate(generated, schema)))
    stopifnot(identical(as.character(generated), as.character(result)))
    setwd(old)

    missing <- metadata
    missing[1, "URL"] <- ""
    write.dcf(missing, file.path(path, "DESCRIPTION"))
    expect_error(process(parse_blocks(annotations), path), "DESCRIPTION field URL")
})

expect_error(parse_blocks("@edam_input_data mzML (format_3244)"),
             "requires Label (data_NNNN)")
expect_error(parse_blocks("@edam_operation Data parsing (operation_181)"),
             "requires Label (operation_NNNN)")
expect_error(parse_blocks("@edam_operation operation_1812"),
             "requires Label (operation_NNNN)")
expect_error(parse_blocks("@edam_typo anything"), "is not supported")
expect_error(process(parse_blocks("@edam_topic Proteomics (topic_0121)")),
             "at least one @edam_operation")
expect_error(process(parse_blocks("@edam_input_data Mass spectrometry data (data_2536)")),
             "requires @edam_operation")
expect_error(process(parse_blocks(c(
    "@edam_operation Data parsing (operation_1812)",
    "@edam_input_format mzML (format_3244)"))), "must follow @edam_input_data")
expect_error(process(parse_blocks(c(
    "@edam_operation Data parsing (operation_1812)",
    "@biotools_documentation User manual https://example.org/"))),
    "requires Type | URL")
expect_error(process(parse_blocks(c(
    "@edam_operation Data parsing (operation_1812)",
    "@biotools_id first", "@biotools_id second"))), "Only one @biotools_id")

checked_in <- biotools_read_xml("inst/biotools.xml")
stopifnot(identical(as.character(checked_in), as.character(doc)))
cat("bio.tools metadata checks passed\n")
