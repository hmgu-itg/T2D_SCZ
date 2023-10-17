################################################################################
#---------------------- Child versus Adulthood BMI -----------------------------

library(ggplot2)

project_folder <- "C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/SCZ_T2D/"
setwd(project_folder)
mr <- data.table::fread("causal_inference/child_adult_bmi/univariate_child_bmi_analysis_MR.csv")
mr <- rbind(mr, data.table::fread("causal_inference/child_adult_bmi/univariate_adult_bmi_analysis_MR.csv"))
mvmr <- data.table::fread("causal_inference/child_adult_bmi/multivariate_adult_child_bmi_analysis_MR.csv")

dt <- data.table::data.table(analysis=rep(c("univariate", "univariate", "multivariate", "multivariate"), 2),
                             exposure=rep(c("childhood BMI", "adulthood BMI"), 4),
                             outcome=c(rep("type 2 diabetes",4), rep("schizophrenia",4)),
                             OR=c(mr[exposure=="child_bmi.european" & outcome=="T2D.multiancestry", `or_Inverse variance weighted`],
                                  mr[exposure=="adult_bmi.european" & outcome=="T2D.multiancestry", `or_Inverse variance weighted`],
                                  mvmr[exposure=="child_bmi.european" & outcome=="T2D.multiancestry", or],
                                  mvmr[exposure=="adult_bmi.european" & outcome=="T2D.multiancestry", or],
                                  mr[exposure=="child_bmi.european" & outcome=="SCZ.multiancestry", `or_Inverse variance weighted`],
                                  mr[exposure=="adult_bmi.european" & outcome=="SCZ.multiancestry", `or_Inverse variance weighted`],
                                  mvmr[exposure=="child_bmi.european" & outcome=="SCZ.multiancestry", or],
                                  mvmr[exposure=="adult_bmi.european" & outcome=="SCZ.multiancestry", or]),
                             lower_OR=c(mr[exposure=="child_bmi.european" & outcome=="T2D.multiancestry", `or_lci95_Inverse variance weighted`],
                                        mr[exposure=="adult_bmi.european" & outcome=="T2D.multiancestry", `or_lci95_Inverse variance weighted`],
                                        mvmr[exposure=="child_bmi.european" & outcome=="T2D.multiancestry", or_lci95],
                                        mvmr[exposure=="adult_bmi.european" & outcome=="T2D.multiancestry", or_lci95],
                                        mr[exposure=="child_bmi.european" & outcome=="SCZ.multiancestry", `or_lci95_Inverse variance weighted`],
                                        mr[exposure=="adult_bmi.european" & outcome=="SCZ.multiancestry", `or_lci95_Inverse variance weighted`],
                                        mvmr[exposure=="child_bmi.european" & outcome=="SCZ.multiancestry", or_lci95],
                                        mvmr[exposure=="adult_bmi.european" & outcome=="SCZ.multiancestry", or_lci95]),
                             upper_OR=c(mr[exposure=="child_bmi.european" & outcome=="T2D.multiancestry", `or_uci95_Inverse variance weighted`],
                                        mr[exposure=="adult_bmi.european" & outcome=="T2D.multiancestry", `or_uci95_Inverse variance weighted`],
                                        mvmr[exposure=="child_bmi.european" & outcome=="T2D.multiancestry", or_uci95],
                                        mvmr[exposure=="adult_bmi.european" & outcome=="T2D.multiancestry", or_uci95],
                                        mr[exposure=="child_bmi.european" & outcome=="SCZ.multiancestry", `or_uci95_Inverse variance weighted`],
                                        mr[exposure=="adult_bmi.european" & outcome=="SCZ.multiancestry", `or_uci95_Inverse variance weighted`],
                                        mvmr[exposure=="child_bmi.european" & outcome=="SCZ.multiancestry", or_uci95],
                                        mvmr[exposure=="adult_bmi.european" & outcome=="SCZ.multiancestry", or_uci95]))
dt <- dt[analysis=="univariate", exposure_extra:=paste(exposure, "univariate", sep="_")]
dt <- dt[analysis=="multivariate", exposure_extra:=paste(exposure, "multivariate", sep="_")]



ggplot(dt, aes(x=OR, y=exposure_extra)) +
  geom_point(aes(color=exposure, shape=analysis), size=4) + #, shape=19) +
  geom_errorbarh(aes(xmin=lower_OR, xmax=upper_OR), height=.3) +
  # coord_fixed(ratio=.3) +
  geom_vline(xintercept=1, linetype='dashed') +
  theme_bw() +
  scale_color_manual(values=RColorBrewer::brewer.pal(9, "RdPu")[c(5,8)]) +
  facet_wrap(~outcome, ncol=1) +
  theme(legend.title=element_blank(), legend.position="bottom", strip.text = element_text(size = 10), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("OR per unit increase in BMI (95% CI)") + ylab("")

ggsave("causal_inference/child_adult_bmi/forest_plot.png", width=5, height=4, dpi=600)
