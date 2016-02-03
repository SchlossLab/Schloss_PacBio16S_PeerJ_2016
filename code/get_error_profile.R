library(rjson)

mock_error <- read.table(file="data/process/mock.error.report")
mock_error_nochim <- mock_error[mock_error$numparents==1,]
mock_error_good <- mock_error_nochim[mock_error_nochim$error <= 0.10,]

mutations <- c("AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC", "GA", "GT", "GG", "GC", "CA", "CT", "CG", "CC", "Ai", "Ti", "Gi", "Ci", "dA", "dT", "dG", "dC")

mutation_table <- mock_error_good[, which(colnames(mock_error_good)%in%mutations)]
mutation_count <- apply(mutation_table, 2, sum)

mutation_total <- sum(mutation_count)
matches <- sum(mutation_count[c("AA", "TT", "GG", "CC")])
mis_matches <- mutation_total - matches
init_error_rate <- 100 * mis_matches / mutation_total		#0.6514488


#substitutions
substitutions <- mutations[grep("[ATGC].", mutations)]
substitutions <- substitutions[grep(".[ATGC]", substitutions)]
substitutions <- substitutions[! substitutions %in% c("AA", "TT", "GG", "CC")]
substitution_counts <- mutation_count[substitutions]
substitution_total <- sum(substitution_counts)
substitution_bias <- substitution_counts / substitution_total

reference <- gsub(".(.)", "\\1", names(substitution_bias))
# aggregate(substitution_bias, by=list(reference), mean)
#
# Group.1          x
# 1       A 0.10192188
# 2       C 0.05944286
# 3       G 0.07434043
# 4       T 0.09762816

substitution_rate <- 100 * substitution_total / mis_matches

# [1] 50.90346


#deletions
deletion_counts <- mutation_count[grep("d.", mutations)]
deletion_total <- sum(deletion_counts)
deletion_bias <- 100 * deletion_counts / deletion_total

#dA        dT        dG        dC
#24.25770 17.99132 39.44415 18.30683

deletion_rate <- 100 * deletion_total / mis_matches

#[1] 17.93361

#insertions
insertion_counts <- mutation_count[grep(".i", mutations)]
insertion_total <- sum(insertion_counts)
insertion_bias <- 100 * insertion_counts / insertion_total

#Ai        Ti        Gi        Ci
#23.13509 19.63951 30.05906 27.16633

insertion_rate <- 100 * insertion_total / mis_matches

#[1] 31.16293



# Need to connect error types with the quality scores
error_quality <- read.table(file="data/process/mock.quality.report", header=T, row.names=1)

#get quantiles
matches <- rep(as.numeric(rownames(error_quality)), error_quality$matches)
m_quant <- quantile(matches, prob=c(0.025, 0.25, 0.5, 0.75, 0.975))
m_93 <- 100 * sum(matches >= 93)/length(matches)
# [1] 80.52649

subs <- rep(as.numeric(rownames(error_quality)), error_quality$substitutions)
s_quant <- quantile(subs, prob=c(0.025, 0.25, 0.5, 0.75, 0.975))
s_93 <- 100 * sum(subs >= 93)/length(subs)
# [1] 79.96451

ins <- rep(as.numeric(rownames(error_quality)), error_quality$insertions)
i_quant <- quantile(ins, prob=c(0.025, 0.25, 0.5, 0.75, 0.975))
i_93 <- 100 * sum(ins >= 93)/length(ins)
# [1] 80.35283

json <- toJSON(
	list(
		overall_error=init_error_rate, subs_rate=substitution_rate, del_rate=deletion_rate, ins_rate=insertion_rate,
		del_bias=deletion_bias,
		ins_bias=insertion_bias,
		percent_over_93 = c(m=m_93, s=s_93, i=i_93),
		qual_range = range(matches, subs, ins),
		pred_error_corr = cor(mock_error_good$error, mock_error_good$pred_error)
	)
)

write(json, "data/process/error_profile.json")
