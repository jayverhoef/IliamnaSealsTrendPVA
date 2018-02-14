\name{iliamna_counts}
\alias{iliamna_counts}
\docType{data}
\title{
Count data for Illiamna Seals
}
\description{
A data.frame of count data for Illiamna Seals
}
\details{
  A data frame with 1200 observations on the following 15 variables.
  \describe{
    \item{datetime}{a POSIXct}
    \item{observer}{a factor with levels \code{Kline} \code{London} \code{Mathisen/Kline} \code{Pebble Mine} \code{Small} \code{Wells} \code{Withrow} \code{Withrow/Askoak} \code{Withrow/Burns} \code{Withrow/Burns/Askaok} \code{Withrow/Burns/Askaok/Seagers} \code{Withrow/Burns/Askaok/Van Landen} \code{Withrow/Dahle} \code{Withrow/Kugo/Wynne} \code{Withrow/Yano} \code{Withrow/Yano/Gomez} \code{Withrow/Yano/London}}
    \item{polyid}{a factor with levels \code{JF15} \code{JF16} \code{JF18} \code{JF19} \code{Other}}
    \item{surveyed}{a numeric vector}
    \item{waypoint}{a factor with levels \code{JF15A} \code{JF15B} \code{JF15C} \code{JF16A} \code{JF16CD} \code{JF16EF} \code{JF16G} \code{JF16H} \code{JF16J} \code{JF16K} \code{JF16L} \code{JF16M} \code{JF18A} \code{JF19A} \code{JF19BD} \code{Other}}
    \item{nonpup_count}{a numeric vector}
    \item{pup_count}{a numeric vector}
    \item{water_count}{a numeric vector}
    \item{adult_total}{a numeric vector}
    \item{pup_total}{a numeric vector}
    \item{air_temp}{a factor with levels \code{} \code{<40} \code{<40s} \code{22} \code{3.2} \code{32} \code{37.4} \code{40s} \code{41} \code{41-50} \code{48.2} \code{50} \code{50s} \code{51-60} \code{51.8} \code{53.6} \code{55.4} \code{57.2} \code{59} \code{60.8} \code{60s} \code{61-70} \code{62.6} \code{66.2} \code{71.6}}
    \item{wind_speed}{a factor with levels \code{} \code{<5} \code{>25} \code{>5} \code{0-5} \code{10-20} \code{10.56} \code{11.81} \code{13.67} \code{14.91} \code{15-25} \code{16-25} \code{17.4} \code{18.64} \code{22.99} \code{30} \code{4.35} \code{45} \code{5-15} \code{5.59} \code{6.84} \code{8.08} \code{9.32} \code{C} \code{Calm} \code{LB}}
    \item{wind_direction}{a factor with levels \code{Calm} \code{E} \code{ESE} \code{N} \code{NNE} \code{NNW} \code{S} \code{SE} \code{SSE} \code{SSW} \code{SW} \code{Variable} \code{W} \code{WNW} \code{WSW}}
    \item{sky_condition}{a factor with levels \code{} \code{30\% Cloud} \code{Clear} \code{Cloudy} \code{Fog} \code{HO} \code{Light Rain} \code{LO} \code{Mostly Cloudy} \code{Overcast} \code{Partly Cloudy} \code{Scattered Clouds}}
    \item{precip}{a factor with levels \code{} \code{Drizzle} \code{No Precip} \code{Showers}}
  }
}
\examples{
data(iliamna_counts)
str(iliamna_counts)
plot(iliamna_counts) ...
}
\keyword{datasets}
