# Development-only roclet: no documentation dependencies are loaded by mzR.
biotools_roclet <- function() {
    ns <- asNamespace("roxygen2")
    edam <- c("topic", "operation", "input_data", "input_format",
              "output_data", "output_format")
    metadata <- c("id", "tool_type", "operating_system", "language",
                  "link", "documentation")
    for (tag in paste0("edam_", edam))
        registerS3method("roxy_tag_parse", paste0("roxy_tag_", tag),
                         biotools_parse_edam, envir = ns)
    for (tag in paste0("biotools_", metadata))
        registerS3method("roxy_tag_parse", paste0("roxy_tag_", tag),
                         biotools_parse_text, envir = ns)
    # roxygen2 7.x tokenizes underscore tags as @edam / @biotools plus raw text.
    registerS3method("roxy_tag_parse", "roxy_tag_edam",
                     biotools_parse_underscore, envir = ns)
    registerS3method("roxy_tag_parse", "roxy_tag_biotools",
                     biotools_parse_underscore, envir = ns)
    registerS3method("roclet_process", "roclet_biotools",
                     biotools_process, envir = ns)
    registerS3method("roclet_output", "roclet_biotools",
                     biotools_output, envir = ns)
    roxygen2::roclet("biotools")
}

biotools_parse_underscore <- function(x) {
    parts <- regmatches(x$raw, regexec("^(_[a-z_]+)[[:space:]]+(.+)$",
                                      x$raw))[[1]]
    if (length(parts) != 3L)
        biotools_tag_error(x, "requires a supported underscore tag and value")
    x$tag <- paste0(x$tag, parts[2])
    x$raw <- parts[3]
    parser <- getS3method("roxy_tag_parse", paste0("roxy_tag_", x$tag),
                          optional = TRUE, envir = asNamespace("roxygen2"))
    if (is.null(parser))
        biotools_tag_error(x, "is not supported")
    parser(x)
}

biotools_tag_error <- function(x, message) {
    stop(x$file, ":", x$line, ": @", x$tag, " ", message, call. = FALSE)
}

biotools_parse_text <- function(x) {
    x$val <- trimws(gsub("[[:space:]]+", " ", x$raw))
    if (!nzchar(x$val))
        biotools_tag_error(x, "requires a value")
    x
}

biotools_parse_edam <- function(x) {
    x <- biotools_parse_text(x)
    branch <- sub("^edam_(input_|output_)?", "", x$tag)
    pattern <- paste0("^(.+) \\((", branch, "_[0-9]{4})\\)$")
    parts <- regmatches(x$val, regexec(pattern, x$val))[[1]]
    if (length(parts) != 3L)
        biotools_tag_error(x, paste0("requires Label (", branch, "_NNNN)"))
    x$val <- list(uri = paste0("http://edamontology.org/", parts[3]),
                  term = parts[2])
    x
}

biotools_values <- function(blocks, tag) {
    unique(unlist(lapply(blocks, function(block) {
        lapply(roxygen2::block_get_tags(block, tag), `[[`, "val")
    }), recursive = FALSE))
}

biotools_add_terms <- function(parent, name, terms) {
    for (term in terms) {
        node <- xml2::xml_add_child(parent, name)
        xml2::xml_add_child(node, "uri", term$uri)
        xml2::xml_add_child(node, "term", term$term)
    }
}

biotools_add_io <- function(parent, block, direction) {
    tags <- roxygen2::block_get_tags(
        block, paste0("edam_", direction, c("_data", "_format")))
    node <- NULL
    formats <- character()
    for (tag in tags) {
        if (endsWith(tag$tag, "_data")) {
            node <- xml2::xml_add_child(parent, direction)
            biotools_add_terms(node, "data", list(tag$val))
            formats <- character()
        } else {
            if (is.null(node))
                biotools_tag_error(tag, paste0(
                    "must follow @edam_", direction, "_data in the same block"))
            if (!tag$val$uri %in% formats)
                biotools_add_terms(node, "format", list(tag$val))
            formats <- c(formats, tag$val$uri)
        }
    }
}

biotools_add_links <- function(parent, name, values) {
    for (value in unique(values)) {
        parts <- trimws(strsplit(value, "|", fixed = TRUE)[[1]])
        if (length(parts) != 2L || any(!nzchar(parts)))
            stop("@biotools_", name, " requires Type | URL", call. = FALSE)
        node <- xml2::xml_add_child(parent, name)
        xml2::xml_add_child(node, "url", parts[2])
        if (name == "documentation" && parts[1] == "Package documentation") {
            xml2::xml_add_child(node, "type", "General")
            xml2::xml_add_child(node, "note", parts[1])
        } else {
            xml2::xml_add_child(node, "type", parts[1])
        }
    }
}

biotools_add_credits <- function(parent, description) {
    if (!"Authors@R" %in% names(description))
        return(invisible(NULL))
    authors <- eval(parse(text = description[["Authors@R"]]),
                    envir = list(person = utils::person), enclos = baseenv())
    roles <- c(aut = "Developer", cre = "Primary contact", ctb = "Contributor",
               cph = "Provider", fnd = "Support", rev = "Contributor",
               ths = "Contributor", trl = "Contributor")
    for (i in seq_along(authors)) {
        author <- authors[i]
        node <- xml2::xml_add_child(parent, "credit")
        xml2::xml_add_child(node, "name",
                            paste(c(author$given, author$family), collapse = " "))
        if (length(author$email))
            xml2::xml_add_child(node, "email", author$email)
        if ("ORCID" %in% names(author$comment))
            xml2::xml_add_child(node, "orcidid", paste0(
                "https://orcid.org/", author$comment[["ORCID"]]))
        xml2::xml_add_child(node, "typeEntity", "Person")
        for (role in unique(unname(roles[intersect(author$role, names(roles))])))
            xml2::xml_add_child(node, "typeRole", role)
    }
}

biotools_process <- function(x, blocks, env, base_path) {
    description <- read.dcf(file.path(base_path, "DESCRIPTION"))[1, ]
    for (field in c("Package", "Description", "URL")) {
        if (!field %in% names(description) || !nzchar(description[[field]]))
            stop("bio.tools requires DESCRIPTION field ", field, call. = FALSE)
    }
    text <- function(value) trimws(gsub("[[:space:]]+", " ", value))
    document <- xml2::xml_new_root("tools", xmlns = "biotoolsSchema")
    tool <- xml2::xml_add_child(document, "tool")
    xml2::xml_add_child(tool, "name", description[["Package"]])
    xml2::xml_add_child(tool, "description", text(description[["Description"]]))
    urls <- trimws(strsplit(description[["URL"]], "[,\n]")[[1]])
    xml2::xml_add_child(tool, "homepage", urls[1])
    ids <- biotools_values(blocks, "biotools_id")
    if (length(ids) > 1L)
        stop("Only one @biotools_id is allowed", call. = FALSE)
    if (length(ids))
        xml2::xml_add_child(tool, "biotoolsID", ids[[1]])
    if ("Version" %in% names(description))
        xml2::xml_add_child(tool, "version", description[["Version"]])
    for (value in biotools_values(blocks, "biotools_tool_type"))
        xml2::xml_add_child(tool, "toolType", value)
    biotools_add_terms(tool, "topic", biotools_values(blocks, "edam_topic"))
    for (value in biotools_values(blocks, "biotools_operating_system"))
        xml2::xml_add_child(tool, "operatingSystem", value)
    for (value in biotools_values(blocks, "biotools_language"))
        xml2::xml_add_child(tool, "language", value)
    if ("License" %in% names(description))
        xml2::xml_add_child(tool, "license", description[["License"]])
    functions <- 0L
    for (block in blocks) {
        operations <- biotools_values(list(block), "edam_operation")
        io <- roxygen2::block_get_tags(block, paste0(
            "edam_", c("input_data", "input_format", "output_data", "output_format")))
        if (!length(operations)) {
            if (length(io))
                biotools_tag_error(io[[1]], "requires @edam_operation in the same block")
            next
        }
        functions <- functions + 1L
        node <- xml2::xml_add_child(tool, "function")
        biotools_add_terms(node, "operation", operations)
        biotools_add_io(node, block, "input")
        biotools_add_io(node, block, "output")
        note <- biotools_values(list(block), c("title", "description", "details"))
        if (length(note))
            xml2::xml_add_child(node, "note", text(paste(note, collapse = " ")))
    }
    if (!functions)
        stop("bio.tools requires at least one @edam_operation", call. = FALSE)
    links <- biotools_values(blocks, "biotools_link")
    repository <- urls[grepl("^https://github.com/[^/]+/[^/]+/?$", urls)]
    if (length(repository))
        links <- c(links, paste("Repository", repository, sep = " | "))
    if ("BugReports" %in% names(description))
        links <- c(links, paste("Issue tracker", description[["BugReports"]], sep = " | "))
    biotools_add_links(tool, "link", links)
    biotools_add_links(tool, "documentation",
                        biotools_values(blocks, "biotools_documentation"))
    citation_file <- file.path(base_path, "inst", "CITATION")
    if (file.exists(citation_file)) {
        citations <- utils::readCitationFile(citation_file, meta = as.list(description))
        for (citation in unclass(citations)) {
            identifiers <- intersect(c("doi", "pmid", "pmcid"), names(citation))
            if (!length(identifiers))
                next
            node <- xml2::xml_add_child(tool, "publication")
            for (identifier in identifiers)
                xml2::xml_add_child(node, identifier, citation[[identifier]])
        }
    }
    biotools_add_credits(tool, description)
    for (node in xml2::xml_find_all(document, "//*"))
        xml2::xml_set_namespace(node, uri = "biotoolsSchema")
    document
}

biotools_read_xml <- function(path) {
    # The official schema uses a relative namespace, which libxml warns about.
    withCallingHandlers(xml2::read_xml(path, options = "NONET"), warning = function(w) {
        if (startsWith(conditionMessage(w),
                       "xmlns: URI biotoolsSchema is not absolute"))
            invokeRestart("muffleWarning")
    })
}

biotools_output <- function(x, results, base_path, ...) {
    schema <- biotools_read_xml(file.path(base_path, "tools", "biotools", "biotools.xsd"))
    valid <- xml2::xml_validate(results, schema)
    if (!isTRUE(valid)) {
        errors <- attr(valid, "errors")
        stop("Invalid bio.tools XML:\n", paste(errors, collapse = "\n"), call. = FALSE)
    }
    path <- file.path(base_path, "inst", "biotools.xml")
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    xml2::write_xml(results, path, options = "format")
    path
}
