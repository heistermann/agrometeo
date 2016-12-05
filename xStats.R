# Arbeitsumgebung leeren
rm(list=ls())

# Paket fuer Extremwertstatistik importieren (ggf. vorher installieren)
##install.packages("extRemes", repos="https://cloud.r-project.org")
library(extRemes)

# ACHTUNG:
#    Setze Arbeitsverzeichnis via 
#    Session>Set Working Directory>To Source File Location

# ACHTUNG:
#    Setze Start- und Enddatum für die Analyse
startDate = "1961-01-01"
endDate = "1990-12-31"


# Daten lesen
tmp = read.table("precip.csv", header=TRUE, 
                 stringsAsFactors=FALSE, sep=";")

# Daten zu einheitlicher Zeitreihe (Datum - Niederschlag) umformen
date = seq(from=as.Date(startDate), to=as.Date(endDate), by="days")
df = data.frame(date=date, prec=rep(NA, length(date)))
year=strftime(date, format="X%Y")
month=as.numeric(strftime(date, format="%m"))
day=as.numeric(strftime(date, format="%d"))
for (i in 1:length(date)) {
  rowix = which((tmp$Monat==month[i]) & (tmp$Jahrestag==day[i]))
  df$prec[i] = as.numeric(tmp[rowix, year[i]])
}
##write.table(df, "precip.reformat.csv", row.names=FALSE, sep=";")

# Zeitreihe als Diagramm darstellen
plot(df$date,df$prec, type="l", xlab="Zeit", 
     ylab="Niederschlag (mm/d)")

#################################################################
# EXTREMWERTSTATISTIK FUER NIEDERSCHLAEGE
#################################################################

# Monatssummen
df$year = as.numeric(strftime(df$date,format="%Y"))
df$month = as.numeric(strftime(df$date,format="%m"))
monthly = aggregate(df$prec, 
                    by=list(year=df$year, month=df$month),
                    FUN=sum)
names(monthly) = c("year", "month", "prec")

# Monats- und Tagessummen ohne Fehlwerte
daily = df[!is.na(df$prec),]
monthly = monthly[!is.na(monthly$prec),]

# Berechne Jahresmaxima der taegl. und monatl. Summen
daymax   = blockmaxxer(x=daily, which="prec", blocks = daily$year)
monthmax = blockmaxxer(x=monthly, which="prec", blocks = monthly$year)

# Darstellung der Jahresmaxima als Histogramme
op = par(mfrow=c(2,1), mar=c(5, 4, 4, 2) + 0.1)
hist(daymax$prec, main="Annual max. of daily precipitation", xlab="Precipitation (mm/day)")
hist(monthmax$prec, main="Annual max. of monthly precipitation", xlab="Precipitation (mm/month)")
par(op)

# Anpassung der GEV-Verteilung 
#    (hier an Jahresmaxima der Monatssummen) 
fitGEV <- fevd(prec, data = monthmax)

# Auswertung der Anpassung
summary(fitGEV)
rperiods = c(2, 5, 10, 20, 50, 100)
plot(fitGEV, rperiods = rperiods)
return.level(fitGEV, return.period = rperiods, do.ci=TRUE)

# Diese Funktion gibt fuer einen Monatsniederschlag das
#   Wiederkehrintervall x wieder (Ueberschreitungsw.keit=1/x)
return.period.from.level = function(fitobj, level) {
  p = pevd(level, 
       loc=fitobj$results$par["location"],
       scale=fitobj$results$par["scale"],
       shape=fitobj$results$par["shape"],
       type=fitobj$type)
  return(1. / (1-p) )
}
return.period.from.level(fitGEV, c(80, 100, 120, 140, 160))

#################################################################
# ANALYSE VON TROCKPERIODEN
#################################################################

# Diese Funktion gibt einen data.frame von Trockenperioden zurueck
get.dry.spells = function(dates, prec, minlen, dry) {
  start = c()
  end = c()
  dry.day = prec <= dry
  i = 1
  j = i
  new.search=FALSE
  while (1) {
    # Abbruchkriterium
    if (i>=length(prec)) break
    
    if(!is.na(prec[j])) {
      # j ist kein NA-Wert
      if (dry.day[j]) {
        # j ist ein Trockentag
        j = j+1
      } else new.search=TRUE
    } else new.search=TRUE
    
    if (new.search) {
      # Auswertung nach Ende einer Suche
      if (j-i >= minlen) {
        start = c(start, i)  
        end = c(end, j-1)  
      }
      # Weitersuchen ab hier
      i = j+1
      j = i
      new.search = FALSE
    }
  }
  return(data.frame(start=dates[start],
                    end=dates[end],
                    length=end-start+1))
}

# Die detektierten Trockenperioden haengen ab von der
#   Mindestdauer (minlen) und der Definition eines
#   trockenen Tages (Niederschlag <= dry, in mm)
dry.spells = get.dry.spells(df$date, df$prec, minlen=10, dry=0.1)

# Gesamte Verteilung von Trockenperioden
hist(dry.spells$length, breaks=seq(10,35,5),
     xlab="Dauer (in Tagen)",
     ylab="Haeufigkeit",
     main="Haeufigkeit von Trockenperioden > 10 Tage"
)

# GROBE ANALYSE: Mittlere monatliche Verteilung von Trockenperioden
#   Eine Trockenperiode wird demjenigen Monat zugeteilt,
#   der den groessten Anteil an ihr hat
dry.spells$month = strftime(dry.spells$start 
                            + dry.spells$length/2,
                            format="%m")
dry.spells$month = as.numeric(dry.spells$month)

# Haeufigkeit von Trockenperioden bezogen auf Monat
plot(table(dry.spells$month), type="b",
     xlab="Monat", ylab="Haeufigkeit")
table(dry.spells$month)

# BESSER:
#   Zahl von Tagen eines Monats, die zu Trockenperioden gehoeren
#   (n.drydays) 
n.drydays = rep(0,12)
for (i in 1:nrow(dry.spells)) {
  startmonth = as.numeric( strftime(dry.spells$start[i], format="%m") )
  endmonth = as.numeric( strftime(dry.spells$end[i], format="%m") )
  if (startmonth==endmonth) {
    # Wenn Anfang und Ende in einem Monat liegen...
    n.drydays[startmonth] = n.drydays[startmonth] + dry.spells$length[i]
  } else {
    date = dry.spells$start[i]
    count = 0
    month = startmonth
    while(date <= dry.spells$end[i]) {
      if (month !=  as.numeric(strftime(date, format="%m") ) ) {
        n.drydays[month] = n.drydays[month] + count
        month = as.numeric(strftime(date, format="%m"))
        count = 0
      }
      if ( date == dry.spells$end[i] ) {
        n.drydays[month] = n.drydays[month] + count + 1
        break
      }
      count = count + 1
      date = date + 1
    }
  }  
}

# Monatsmittelwert der Zahl von Tagen in Trockenperioden
n.years = max(df$year) - min(df$year) + 1
avg.n.drydays = n.drydays / n.years
plot(1:12, avg.n.drydays, type="b", 
     xlab="Monat", ylab="Zahl der Tage",
     main="Zahl der Tage in Trockenperioden")
