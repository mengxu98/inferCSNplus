extract_genes <- function(
    network_list,
    regulators = NULL,
    targets = NULL,
    top_regulators_num = 10,
    top_targets_num = 10,
    common_regulators = FALSE) {

  if (common_regulators) {
    regulators_list <- purrr::map(
      network_list, function(x) {
        regulators <- x$regulator
        regulators[!duplicated(regulators)]
      }
    )
    regulators_common <- purrr::reduce(
      regulators_list,
      .f = intersect
    )
  }

  if (is.list(network_list)) {
    genes_list <- purrr::map(
      network_list, function(x) {
        regulators_list <- x$regulator
        regulators_list <- regulators_list[!duplicated(regulators_list)]

        if (common_regulators) {
          if (is.null(regulators)) {
            regulators <- regulators_common[1:min(top_regulators_num, length(regulators_common))]
          }
        } else {
          if (is.null(regulators)) {
            regulators <- regulators_list[1:min(top_regulators_num, length(regulators_list))]
          }
        }

        targets_list <- purrr::map_dfr(
          regulators, function(g) {
            x_g <- net.format(
              x,
              regulators = g,
              targets = targets,
              abs_weight = FALSE
            )
            x_g[1:top_targets_num, ]
          }
        )
        targets <- unique(targets_list$target)

        return(list(
          regulators = regulators,
          targets = targets
        ))
      }
    )
  }

  return(genes_list)
}


setdiff_genes <- function(list) {
  genes_list <- list()
  for (i in seq_len(length(list))) {
    genes1 <- list[[i]]
    nums_list <- seq_len(length(list))[-i]
    list_new <- list[nums_list]
    genes2 <- purrr::list_c(list_new)
    genes_list[[i]] <- unique(setdiff(genes1, genes2))
  }

  return(genes_list)
}


extract_genes_transition <- function(
    network_list,
    regulators = NULL,
    targets = NULL,
    transition_networks = NULL,
    top_regulators_num = 5,
    top_targets_num = 5,
    common_regulators = TRUE,
    common_targets = TRUE) {
  if (is.null(transition_networks)) {
    network_names <- sort(names(network_list))
    transition_networks <- purrr::map(
      seq_len(((length(network_names) - 1) / 2)),
      .f = function(x) {
        x <- 2 * x - 1
        network_names[x:(x + 2)]
      }
    )
  }

  transition_networks_list <- purrr::map(
    transition_networks,
    .f = function(x) {
      network_list_sub <- network_list[x]
      transition_networks_list_sub <- extract_genes(
        network_list_sub,
        common_regulators = common_regulators,
        top_regulators_num = top_regulators_num,
        top_targets_num = top_targets_num
      )
      if (common_regulators && common_targets) {
        targets_list <- purrr::map(
          transition_networks_list_sub, function(x) {
            targets <- x$target
            targets[!duplicated(targets)]
          }
        )
        targets_common <- purrr::reduce(
          targets_list,
          .f = intersect
        )

        targets_setdiff <- setdiff_genes(targets_list)
        names(targets_setdiff) <- names(targets_list)

        transition_networks_list_sub <- purrr::map2(
          transition_networks_list_sub,
          targets_setdiff,
          .f = function(x, y) {
            list(regulators = x$regulators, targets = y)
          }
        )
        transition_networks_list_sub[["intersect"]] <- list(
          regulators = transition_networks_list_sub[[1]][1],
          targets = targets_common
        )
      }

      return(transition_networks_list_sub)
    }
  )
  names(transition_networks_list) <- paste0(
    "transition_network",
    seq_len(length(transition_networks))
  )

  return(transition_networks_list)
}
