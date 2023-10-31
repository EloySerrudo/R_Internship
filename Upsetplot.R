library(tidyverse)
library(ggplot2movies)
library(ComplexUpset)

movies_c <- movies

dim(movies_c)
head(movies_c)

genres <- colnames(movies_c)[18:24]
genres

movies_c[genres] <- movies_c[genres] == 1 # Turn the genre columns from 1 to TRUE
t(head(movies_c[genres], 3))

movies_c[movies_c$mpaa == '', 'mpaa'] = NA
movies_na = na.omit(movies_c)

set_size(8, 3)
upset(
  data = movies_na, 
  intersect = genres, 
  name = 'Genres', 
  width_ratio=0.1
)

upset(
  movies_na, genres, name='genre', width_ratio=0.1, 
  min_size=10, wrap=TRUE, set_sizes=FALSE
) + 
  ggtitle(
    'Without empty groups (Short dropped)'
  ) +
  upset(movies_na, genres, name='genre', width_ratio=0.1, 
        min_size=10, keep_empty_groups=TRUE, wrap=TRUE, set_sizes=FALSE
  ) + 
  ggtitle('With empty groups')
# adding plots is possible thanks to patchwork

upset(
  movies_na, genres, width_ratio=0.1,
  min_degree=3)
upset(
  movies_na, genres, width_ratio=0.1,
  n_intersections=15
)

# VENN
abc_data <- create_upset_abc_example()

abc_venn <- (
ggplot(arrange_venn(abc_data))
+ coord_fixed()
+ theme_void()
+ scale_color_venn_mix(abc_data)
)

abc_venn + geom_venn_region(data=abc_data, alpha=0.05
                            ) + 
  geom_point(aes(x=x, y=y, color=region), size=1) + 
  geom_venn_circle(abc_data) + 
  geom_venn_label_set(abc_data, aes(label=region))+ 
  geom_venn_label_region(
    abc_data, aes(label=size),
    outwards_adjust=1.75,
    position=position_nudge(y=0.2)
    ) + 
  scale_fill_venn_mix(abc_data, guide='none')

abc_upset = function(mode) upset(
  abc_data, c('A', 'B', 'C'), mode=mode, set_sizes=FALSE,
  encode_sets=FALSE,
  queries=list(upset_query(intersect=c('A', 'B'), color='orange')),
  base_annotations=list(
    'Size'=(
      intersection_size(
        mode=mode,
        mapping=aes(fill=exclusive_intersection),
        size=0,
        text=list(check_overlap=TRUE)
      ) + scale_fill_venn_mix(
        data=abc_data,
        guide='none',
        colors=c('A'='red', 'B'='blue', 'C'='green3')
      )
    )
  )
)

(abc_upset('exclusive_intersection') | abc_upset('inclusive_intersection')) / 
  (abc_upset('exclusive_union') | abc_upset('inclusive_union'))

upset(
  movies_na, genres,
  width_ratio=0.1,
  min_size=10,
  mode='inclusive_union',
  base_annotations=list('Size'=(intersection_size(counts=FALSE, mode='inclusive_union'))),
  intersections='all',
  max_degree=3
)

head(movies_na[-c(7:17)])
