"Reads all enrichment files, merge them and count them"

library(dplyr)

src = './Enrichment'
filter = 'Preserved/Enrichment'

filter = 'IV'

src = paste(src, filter, sep = '/')

files = list.files(src)

mergedAnnotations = data.frame(source = 0,
                               term_name = 0,
                               term_id = 0,
                               adjusted_p_value = 0,
                               negative_log10_of_adjusted_p_value = 0,
                               term_size = 0,
                               query_size = 0,
                               intersection_size = 0,
                               effective_domain_size = 0,
                               intersections = 0
)

for (i in files) {
  file = paste(src, i, sep = '/')
  annotatedGenes = data.frame(read.csv(file))
  
  mergedAnnotations = rbind(mergedAnnotations, annotatedGenes)
}

mergedAnnotations = mergedAnnotations[-1, ]
cnt = mergedAnnotations %>%
  group_by(source) %>%
  count(term_name, sort = TRUE) %>%
  arrange(source)

write.csv(cnt, paste(src, 'Tally.csv', sep = '/'))