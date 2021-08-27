#' @keywords internal
collect_results <- function(mod) {
    results <- parTable(mod$fit) %>%
        filter(label != "") %>%
        filter(label != "var_Crave_t2") %>%
        arrange(label) %>%
        mutate(pvalue = (1-pnorm(abs(est/se)))*1) %>%
        dplyr::select(parameter = label, est, se, `Pr(>|z|)` = pvalue) %>%
        mutate(
            parameter = case_when(
                parameter == "a" ~ "a",
                parameter == "b" ~ "b",
                parameter == "prop_a_b" ~ "a/b",
                parameter == "cprime" ~ "c'",
                parameter == "indirect" ~ "ab",
                parameter == "indirect_std" ~ "ab/sigma(Y)^2",
                parameter == "proportion" ~ "ab/(ab+c')",
                parameter == "ratio" ~ "ab/c'",
                parameter == "rho_m" ~ "rho_M",
                parameter == "rho_y" ~ "rho_Y",
                TRUE ~ "sigma(Y)^2"
            )
        )

    rownames(results) <- results$parameter

    results$parameter <- NULL

    class(results) <- c("anova", "data.frame")

    results <- structure(
        results,
        heading = "results:",
        class = c("anova", "data.frame")
    )

    mod$results <- results

    mod
}

#' @keywords internal
permutation_test <- function(mod) {
    parameters <- c("a", "b", "cprime", "rho_m", "rho_y", "indirect",
                    "indirect_std", "proportion", "ratio")

    not_permuted <- lapply(parameters, function(parameter) {
        tmp <- parTable(mod$fit)[parTable(mod$fit)$label == parameter, ]
        tmp <- c(tmp$est, tmp$se, (1 - pnorm(abs(tmp$est/tmp$se)))*2)
        names(tmp) <- paste0(parameter, c("_est", "_se", "_p"))
        as.data.frame(as.list(tmp))
    }) %>%
        bind_cols() %>%
        mutate(permuted = F)

    results <- replicate(mod$perm_reps, {
        warning <- ""
        result <- withCallingHandlers({
            tryCatch({
                dat_expanded <- mod$dat %>%
                    mutate(!!mod$X := sample(!!as.name(mod$X))) %>%
                    expand_data(mod$var_names, mod$lag)

                fit <- sem(mod$syntax, dat_expanded)

                lapply(parameters, function(parameter) {
                    tmp <- parTable(fit)[parTable(fit)$label == parameter,]
                    tmp <-
                        c(tmp$est, tmp$se, (1 - pnorm(abs(
                            tmp$est / tmp$se
                        ))) * 2)
                    names(tmp) <-
                        paste0(parameter, c("_est", "_se", "_p"))
                    as.data.frame(as.list(tmp))
                }) %>% bind_cols()
            },
            error = function(e) {
                data.frame(error = e$message,
                           stringsAsFactors = F)
            })
        },
        warning = function(w) {
            warning <- w$message
        })
        result$warning <- warning
        result$permuted <- T
        result
    }, simplify = F) %>%
        bind_rows()


    results <- results %>%
        bind_rows(not_permuted) %>%
        mutate(z = indirect_std_est/indirect_std_se) %>%
        summarize(
            `ab/sigma(Y)^2` = indirect_std_est[permuted==F],
            `se(ab/sigma(Y)^2)` = indirect_std_se[permuted==F],
            stat = indirect_std_est[permuted==F]/indirect_std_se[permuted==F],
            `2.5 Pct.` = quantile(z[permuted==T], 0.025, na.rm = T),
            `97.5 Pct.` = quantile(z[permuted==T], 0.975, na.rm = T),
            `Pr(>|stat|)` = mean(abs(z[permuted==T]) > abs(z[permuted==F]), na.rm = T)
        )

    class(results) <- c("anova", "data.frame")

    results <- structure(
        results,
        heading = "permutation test:",
        class = c("anova", "data.frame")
    )

    mod$perm_test <- results

    mod
}

#' @export
sccm <- function(X, M, Y, dat, lag = 1, permutation = F, perm_reps = 1000) {

    var_names <- get_var_names(X, M, Y, lag)
    dat_expanded <- expand_data(dat, var_names, lag)
    syntax <- get_model(var_names)
    fit <- sem(syntax, dat_expanded)

    mod <- list(
        fit = fit,
        var_names = var_names,
        X = X,
        M = M,
        Y = Y,
        lag = lag,
        syntax = syntax,
        dat = dat,
        dat_expanded = dat_expanded,
        perm_reps = perm_reps
    )

    class(mod) <- "sccm"

    mod <- collect_results(mod)

    if (permutation) {
        mod <- permutation_test(mod)
    }

    mod
}

#' @keywords internal
get_var_names <- function(X, M, Y, delay = 1) {
    list(
        X_t = paste0(X, "_t"),
        X_t.m.delay = paste0(X, "_t.m.", delay),
        X_t.m.delay.m.1 = paste0(X, "_t.m.", delay+1),
        M_t = paste0(M, "_t"),
        M_t.m.delay = paste0(M, "_t.m.", delay),
        M_t.m.delay.m.1 = paste0(M, "_t.m.", delay+1),
        Y_t = paste0(Y, "_t"),
        Y_t.m.1 = paste0(Y, "_t.m.1"),
        Y_t.m.delay = paste0(Y, "_t.m.", delay),
        Y_t.m.delay.m.1 = paste0(Y, "_t.m.", delay+1),
        X = X,
        M = M,
        Y = Y
    )
}

#' @keywords internal
expand_data <- function(dat, var_names, delay = 1) {
    dat %>%
        rename(
            !!var_names[["X_t.m.delay.m.1"]] := !!var_names[["X"]],
            !!var_names[["M_t.m.delay.m.1"]] := !!var_names[["M"]],
            !!var_names[["Y_t.m.delay.m.1"]] := !!var_names[["Y"]]
        ) %>%
        mutate(
            !!var_names[["X_t.m.delay"]] := lead(!!as.name(var_names[["X_t.m.delay.m.1"]]), 1),
            !!var_names[["M_t.m.delay"]] := lead(!!as.name(var_names[["M_t.m.delay.m.1"]]), 1),
            !!var_names[["Y_t.m.delay"]] := lead(!!as.name(var_names[["Y_t.m.delay.m.1"]]), 1),
            !!var_names[["X_t"]] := lead(!!as.name(var_names[["X_t.m.delay"]]), delay),
            !!var_names[["M_t"]] := lead(!!as.name(var_names[["M_t.m.delay"]]), delay),
            !!var_names[["Y_t"]] := lead(!!as.name(var_names[["Y_t.m.delay"]]), delay),
            !!var_names[["Y_t.m.1"]] := lead(!!as.name(var_names[["Y_t.m.delay"]]), delay-1)
        )
}

#' @keywords internal
get_model <- function(var_names) {
    paste(
        paste0(var_names[["M_t.m.delay"]], " ~ rho_m*", var_names[["M_t.m.delay.m.1"]], " + a*", var_names[["X_t.m.delay"]]),
        paste0(
            var_names[["Y_t"]],
            " ~ rho_y*",
            var_names[["Y_t.m.1"]],
            " + b*",
            var_names[["M_t.m.delay"]],
            " + cprime*",
            var_names[["X_t.m.delay"]]
        ),
        paste0(var_names[["Y_t"]], " ~~ var_", var_names[["Y_t"]], "*", var_names[["Y_t"]]),
        "## effects",
        "indirect := a*b",
        paste0("indirect_std := a*b/sqrt(var_", var_names[["Y_t"]], ")"),
        "proportion := a*b/(a*b+cprime)",
        "ratio := a*b/(cprime)",
        "prop_a_b := a/b",
        sep = "\n\n"
    )
}

#' @export
plot.sccm <- function(...) {
    get_plot(...)
}

#' @export
get_plot <- function(mod) {

    sccm_dag <- expand.grid(
        name.lhs = c("X", "M", "Y"),
        time.lhs = 0:mod$lag,
        name.rhs = c("X", "M", "Y"),
        time.rhs = 0:mod$lag,
        stringsAsFactors = F
    ) %>%
        as_tibble() %>%
        mutate(name.lhs.full = paste0(name.lhs, "_t", time.lhs)) %>%
        mutate(name.rhs.full = paste0(name.rhs, "_t", time.rhs)) %>%
        # remove X, M, Y if time.lhs > time.rhs
        filter(!(time.lhs > time.rhs)) %>%
        # remove autoregressive effects of X
        filter(!(name.lhs == "X" & name.rhs == "X" & time.lhs != time.rhs)) %>%
        # remove effects with wrong lag length
        filter(!(name.lhs %in% c("X", "M") & name.rhs == "Y" & time.rhs-time.lhs != mod$lag)) %>%
        # remove autoregressive effects of X, M, Y with lag > 1
        filter(!(name.lhs == name.rhs & time.lhs < time.rhs-1)) %>%
        # remove lagged effects from X to M
        filter(!(name.lhs == "X" & name.rhs == "M" & time.lhs != time.rhs)) %>%
        # remove immediate effects from X to Y
        filter(!(name.lhs == "X" & name.rhs == "Y" & time.lhs == time.rhs)) %>%
        # remove immediate effects from M to Y
        filter(!(name.lhs == "M" & name.rhs == "Y" & time.lhs == time.rhs)) %>%
        # remove any effect from M, Y to X
        filter(!(name.lhs %in% c("M", "Y") & name.rhs == "X")) %>%
        # remove any effect from Y to M
        filter(!(name.lhs == "Y" & name.rhs == "M")) %>%
        mutate(y = case_when(name.lhs == "X" ~ 2, name.lhs == "M" ~ 1, name.lhs == "Y" ~ 0)) %>%
        mutate(x = time.lhs)

    coords <- sccm_dag %>%
        filter(name.lhs == name.rhs & time.lhs == time.rhs) %>%
        select(name = name.lhs.full, x, y)

    relations <- sccm_dag %>%
        filter(!(name.lhs == name.rhs & time.lhs == time.rhs)) %>%
        group_by(1:n()) %>%
        group_map(function(x, desc) as.formula(paste0(x$name.rhs.full, " ~ ", x$name.lhs.full)))

    sccm_dag <- do.call(dagify, c(relations, list(coords = coords)))

    sccm_dag <- tidy_dagitty(sccm_dag) %>%
        mutate(lhs = sapply(name, substr, 1, 1)) %>%
        mutate(rhs = sapply(to, substr, 1, 1))

    sccm_dag <- sccm_dag %>%
        mutate(time = as.integer(sapply(name, substr, 4, length(name)))) %>%
        mutate(edge_label = case_when(
            lhs == "X" & rhs == "M" ~ "a",
            lhs == "X" & rhs == "Y" ~ "c'",
            lhs == "M" & rhs == "Y" ~ "b",
            lhs == "M" & rhs == "M" ~ "rho_M",
            lhs == "Y" & rhs == "Y" ~ "rho_Y",
            TRUE ~ NA_character_
        ))

    sccm_dag <- sccm_dag %>%
        mutate(name = paste0(lhs, "_t", ifelse(mod$lag - time == 0, "", paste0(".m.", mod$lag - time)))) %>%
        ggdag() +
        geom_dag_edges_link(
            mapping = aes(label = edge_label),
            angle_calc = "along",
            label_pos = 0.3,
            vjust = -0.2
        ) +
        theme_dag()

    sccm_dag
}

#' @export
summary.sccm <- function(mod) {
    print(mod$results)
    cat("\n\n")
    if (!is.null(mod$perm_test)) {
        (mod$perm_test)
    }
}

#' @export
#' @importFrom lavaan summary
print.sccm <- function(mod) {
    summary(mod$fit, fit.measures = T)
}

#' @export
sccm_gui <- function() {
    # Run the application
    shinyApp(ui = ui, server = server)
}
