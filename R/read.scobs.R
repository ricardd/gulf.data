#' @title Read Observer Data
#' 
#' @description Function to access and read at-sea and port observer data.
#' 
#' @param year Year.
#' @param file File name(s).
#' @param path Data file path.
#' @param cfvn Canadian Fishing Vessel Number(s).
#' @param type Sampling type ('sea' or 'port').
#' @param trip.number Trip number identifier.
#' @param trip.number Trip number identifier.
#' @param database Oracle database name.
#' @param username Oracle user name.
#' @param password Oracle password.
#' @param source Data source ('r', 'ascii', 'csv', or 'oracle').
#' 

# READ.SCOBS - Read snow crab observer data.
#' @export read.scobs
read.scobs <- function(source = "R", ...){
   # Parse 'source' argument:
   source <- match.arg(tolower(source), c("r", "ascii", "csv", "oracle"))

   # Read R or csv file:
   if (source == "r")      x <- read.rdata.scobs(...)
   if (source == "csv")    x <- read.csv.scobs(...) 
   if (source == "ascii")  x <- read.ascii.scobs(...)
   if (source == "oracle") x <- read.oracle.scobs(...)

   # Convert to 'scbio' object:
   x <- scobs(x)

   return(x)
}

# Read GAP database:
read.oracle.scobs <- function(x, year, zone, cfvn, type = "sea", trip.number, 
                              dsn = options()$gulf.oracle$gap$dsn, 
                              uid = options()$gulf.oracle$gap$uid, ...){

   # Parse 'year' argument:
   if (missing(year)) year <- substr(as.character(Sys.Date()), 1, 4)
   year.str <- paste0("'", paste0(year, collapse = "', '"), "'")
   
   # Parse 'type' argument:
   type <- match.arg(tolower(type), c("sea", "port"))
   type.str <- paste0("'", paste0(type, collapse = "', '"), "'")
   
   # Parse 'zone'oracle argument:
   if (missing(zone)) zone <- c('12', '18', '19', 'E', '12E', '12F', 'F')
   zone <- toupper(zone)
   if ("E" %in% zone) zone <- unique(c("12E", zone))
   if ("F" %in% zone) zone <- unique(c("12F", zone))
   if ("12E" %in% zone) zone <- unique(c("E", zone))
   if ("12F" %in% zone) zone <- unique(c("F", zone))
   zone.str <- paste0("'", paste0(zone, collapse = "', '"), "'")
   
   # Parse cfvn:
   if (!missing(cfvn)){
      cfvn <- as.numeric(cfvn)
      cfvn <- cfvn[!is.na(cfvn)]
      cfvn.str <- paste0("'", paste0(cfvn, collapse = "', '"), "'")
   }
   
   # Define default query:
   query <- "select  -- decode(a.NO_SEQ_CASIER, 99, 'P', 'S') || mod(extract(year from a.dat_deb_acti_pech), 100) || a.NO_SORTIE filename
          a.dat_deb_acti_pech date_sampled
         ,a.cod_typ_casier trap_type
         ,a.NO_VOY_OBSR || a.NO_AFFEC || a.NO_SORTIE trip_number
         ,d.NO_INDIV crab_number
         ,decode(D.COD_SEXE, 'M', '1', 'F', '2', d.COD_SEXE) SEX
         ,d.LARG_CARAP*1000 carapace_width
         ,a.VAL_DUR_IMME_CASIER soak_days
         ,a.NO_SEQ_CASIER trap_number
         ,d.HAUT_PINC*1000 chela_height_right
         ,null chela_height_left
         ,d.COD_COND_CARAP shell_condition
         ,d.VAL_DURETE_CARAP durometer
         ,F.NO_BPC_PROB cfvn
         ,null comments
         ,a.UO_VAL_LATI_DEB LL
         ,a.VAL_LONGI_DEB longitude
         ,a.VAL_LATI_DEB latitude
         ,a.VAL_PROFD_MOY DEPTH
         ,null chela_length
         ,null weight
         ,a.cod_zone_pech_gest zone
         ,decode(SIGN(a.NO_SEQ_CASIER-50), 1, 1, 2)  data_type -- ??
         ,nvl(replace(replace(replace(d.VAL_ETAT_PATTE, ' ', '*') , 'M', '1'), 'R', '2'), '**********') missing_legs
         ,v.NOM_OBSR Observer
         ,b.NOM_BPC_PROB Vessel
         ,null egg_maturity
         ,null egg_development
         ,null egg_clutch
         ,a.NB_MALE males_in_trap
         ,a.COD_QUADRIL  grid
         ,a.NB_FEM females_in_trap
         ,s.NO_JRN_BORD logbook
         ,null  Abdomen_width
         ,null MaillageC
         ,a.LARG_COTE_MAILLE*1000 MAILLAGE_C1
         ,a.LARG_COTE_MAILLE_2*1000 MAILLAGE_C2
         ,a.LARG_COTE_MAILLE_3*1000 MAILLAGE_C3
         ,null MaillageC1
         ,null MaillageC2
         ,null MaillageC3
         ,s.PDS_VIF_ESTIM_DEBAR landed_weight_kg
      from GAP.acti_peche_casier a, GAP.detail_peche_casier d, GAP.SORTIE_MER s, GAP.AFFECTATION f, GAP.VOYAGE_OBSERVATEUR v, GAP.BATEAU_PROB b
      where a.cod_esp_prob=2526
         and a.NO_VOY_OBSR=d.NO_VOY_OBSR
         and a.NO_AFFEC=d.NO_AFFEC
         and a.NO_SORTIE=d.NO_SORTIE
         and a.NO_SEQ_CASIER=d.NO_SEQ_CASIER
         and a.NO_VOY_OBSR=s.NO_VOY_OBSR
         and a.NO_AFFEC=s.NO_AFFEC
         and a.NO_SORTIE=s.NO_SORTIE
         and s.NO_VOY_OBSR=f.NO_VOY_OBSR
         and s.NO_AFFEC=f.NO_AFFEC
         and f.NO_VOY_OBSR=v.NO_VOY_OBSR
         and f.NO_BPC_PROB=b.NO_BPC_PROB
         and f.AN_EFFEC_BPC_PROB=b.AN_EFFEC_BPC_PROB
         and a.cod_zone_pech_gest in ('12E')
         and extract (year from a.dat_deb_acti_pech) in (YY99);"
   
   # Refine search using input arguments:
   query <- gsub("YY99", year.str, query)
   query <- gsub("'12E'", zone.str, query)
   if (!missing(cfvn)) query <- gsub(";", paste0(" and F.NO_BPC_PROB in ", cfvn.str, ";"), query)
   
   # Use trap number to target sea or port sampling:
   if (type == "sea")  query <- gsub(";", " and a.NO_SEQ_CASIER <= 50;", query, fixed = TRUE)
   if (type == "port") query <- gsub(";", " and a.NO_SEQ_CASIER > 50;", query, fixed = TRUE)
   
   if (!missing(trip.number)){
      trip.str <- paste0("('", paste0(trip.number, collapse = "', '"), "')")
      query <- gsub(";", paste0(" and a.NO_VOY_OBSR || a.NO_AFFEC || a.NO_SORTIE in ", trip.str, ";"), query, fixed = TRUE)
   }
   
   # Open Oracle channel:
   z <- read.oracle(query, dsn = dsn, uid = uid, ...)
   
   # Rename columns:
   names(z) <- tolower(names(z))
   names(z) <- gsub("_", ".", names(z))
   
   z$date.sampled <- substr(z$date.sampled, 1, 10)
   
   vars <- c("trap.type", "crab.number", "sex", "carapace.width", "soak.days", "trap.number",
             "chela.height.right", "chela.height.left", "durometer", "cfvn", "longitude", "latitude", "depth",
             "chela.length", "weight", "data.type", "landed.weight.kg", "males.in.trap", "females.in.trap", "abdomen.width",
             "maillagec", "maillage.c1", "maillage.c2", "maillage.c3", "maillagec1", "maillagec2", "maillagec3")
   vars <- vars[vars %in% names(z)]
   for (i in 1:length(vars)) z[, vars[i]] <- as.numeric(z[, vars[i]])
   
   z$longitude <- -abs(z$longitude)
   # z$data.type <- ifelse(z$data == 1, "port", "sea")
   z$observer <- gsub("^[ ]", "", z$observer)
   z$observer <- gsub("[ ]+", " ", z$observer)
   z$year <- as.numeric(substr(as.character(z$date.sampled), 1, 4))
   z$month <- as.numeric(substr(as.character(z$date.sampled), 6, 7))
   z$day <- as.numeric(substr(as.character(z$date.sampled), 9, 10))
   z$soak.time <- z$soak.days
   z$vessel <- gsub("^[ ]", "", z$vessel)
   z$vessel <- gsub("[ ]+", " ", z$vessel)
   z$position.type <- z$ll
   z$position.type <- "LL"
   z$male.total <- z$males.in.trap
   z$fishing.grid <- z$grid
   z$female.total <- z$females.in.trap
   z$maillageC <- z$maillagec
   z$maillage1 <- z$maillage.c1
   z$maillage2 <- z$maillage.c2
   z$maillage3 <- z$maillage.c3
   z$maillageC1 <- z$maillagec1
   z$maillageC2 <- z$maillagec2
   z$maillageC3 <- z$maillagec3
   z$weight <- z$landed.weight.kg
   z$shell.condition.mossy <- substr(z$shell.condition, 2, 2)
   z$shell.condition <- as.numeric(substr(z$shell.condition, 1, 1))
   z$species <- 2526
   z$position.type <- "LL"
   z$zone <- toupper(z$zone)
   z$zone <- gsub("12E", "E", z$zone)
   z$zone <- gsub("12F", "F", z$zone)
   
   # Extract subset of data:
   vars <- c('day', 'month', 'year', 'trap.type', 'trip.number', 'crab.number', 'sex', 'carapace.width',
             'soak.time', 'trap.number', 'chela.height.right', 'chela.height.left', 'shell.condition',
             'shell.condition.mossy', 'durometer', 'cfvn', 'comments', 'position.type', 'latitude', 'longitude',
             'depth', 'chela.length', 'weight', 'zone', 'data.type', 'missing.legs', 'observer', 'vessel',
             'egg.maturity', 'egg.development', 'egg.clutch', 'male.total', 'fishing.grid', 'female.total',
             'logbook', 'abdomen.width', 'maillageC', 'maillage1', 'maillage2', 'maillage3', 'maillageC1',
             'maillageC2', 'maillageC3', 'species')
   #print(setdiff(vars, intersect(names(z), names(y))))
   z <- z[vars]
   
   # Fix variables:
   z$carapace.width <- round(z$carapace.width)
   
   return(z)
}

# Read CSV files:
read.csv.scobs <- function(x, year, zone, type, path = options("gulf.path")[[1]]$snow.crab$observer, echo = TRUE, ...){
   # Parse 'year':
   if (missing(year) & is.numeric(x)) year <- x
   if (missing(year)) stop("'year' must be specified.")
   
   # Locate files:
   path <- paste0(path, year, "/csv")
   files <- locate(path = path, file = "*.csv")
   
   # Target file subset:
   if (!missing(type)) files <- files[grep(paste0(substr(toupper(type),1,1), year), files)]
   if (!missing(zone)) files <- files[grep(paste0("zone ", tolower(zone)), tolower(files))]
   if (length(files) == 0) return(NULL)
   
   if (length(files) > 0){
      res <- NULL
      for (i in 1:length(files)){
         if (echo) cat(paste0("Reading : '", files[i], "'\n"))
         x <- read.table(files[i], header = TRUE, sep = ",", colClasses = "character", stringsAsFactors = FALSE)
         x <- gulf.utils::recode(x)
         if (!is.null(res)){
            vars <- intersect(names(res), names(x))
            res <- res[vars]
            x <- x[vars]
         }
         res <- rbind(res, x)
      }
      return(res)
   }
}

# Read Rdata files:
read.rdata.scobs <- function(x, year, zone, type, path = options("gulf.path")[[1]]$snow.crab$observer, echo = TRUE, ...){
   # Parse 'year':
   if (!missing(x)) if (missing(year) & is.numeric(x)) year <- x
   if (missing(year)) stop("'year' must be specified.")
      
   # Locate files:
   path <- paste0(path, year, "/R")
   files <- locate(path = path, file = "*.rdata")
   
   # Target file subset:
   if (!missing(type)) files <- files[grep(paste0(substr(toupper(type),1,1), year), files)]
   if (!missing(zone)) files <- files[grep(paste0("zone ", tolower(zone)), tolower(files))]
   if (length(files) == 0) return(NULL)
   
   if (length(files) > 0){
      res <- NULL
      for (i in 1:length(files)){
         if (echo) cat(paste0("Reading : '", files[i], "'\n"))
         load(files[i])
         x <- gulf.utils::recode(x)
         if (!is.null(res)){
            vars <- intersect(names(res), names(x))
            res <- res[vars]
            x <- x[vars]
         }
         res <- rbind(res, x)
      }
      return(res)
   }
}

# Read ASCII files:
read.ascii.scobs <- function(x, ...){
   # Get file format:
   f <- fmt.scobs(x, ...)  
   
   # Define file format:
   if (!missing(x)){
      if (is.character(x)){
         files <- x
      }else{
         # Locate files:
         files <- locate.scobs(x, ...)
         if (length(files) == 0) return(NULL)
      } 
   }else{
      files <- locate.scobs(...)
   }
   
   # Read files: 
   x <- NULL
   for (i in 1:length(files)){
      cat(paste0(i, " of ", length(files), ") ", files[i], "\n"))
      tmp <- utils::read.fwf(file = files[i], 
                             widths = f$width,
                             comment.char = "", fileEncoding = "ISO-8859-1")
      x <- rbind(x, tmp)
   }
   names(x) <- f$name
   
   # Define variables which are to be converted to numeric format:
   vars <- c("day", "month", "year", "crab.number", "sex",
             "carapace.width", "soak.time", "chela.height.right",
             "chela.height.left",  "shell.condition", "durometer", "cfvn",
             "depth",  "chela.length", "weight", "data.type",
             "male.total", "female.total", "abdomen.width", "maillageC",
             "maillage1", "maillage2", "maillage3", "maillageC1", "maillageC2",
             "maillageC3")
   
   vars <- vars[vars %in% names(x)]
   
   # Remove "*" characters and convert to numeric:
   for (i in 1:length(vars)){
      x[, vars[i]] <- gsub("*", " ", x[, vars[i]], fixed = TRUE)
      x[, vars[i]] <- as.numeric(x[, vars[i]], fixed = TRUE)
   }
   
   # Fix for coordinate offset problem in earlier years:
   x$latitude <- gsub(" ", "", paste0(x$blank10, x$latitude))
   x$longitude <- gsub(" ", "", paste0(x$blank11, x$longitude))
   x$latitude[which(gsub("9", "", x$latitude) == "")] <- ""
   x$longitude[which(gsub("9", "", x$longitude) == "")] <- ""
   
   # Remove blank columns:
   x <- x[, setdiff(names(x), names(x)[grep("blank", names(x))])]
   
   # Convert coordinates:
   x$latitude  <- as.numeric(x$latitude)
   x$longitude <- as.numeric(x$longitude)

   # Clean-up 'zone' variable:
   x$zone <- gsub(" ", "", x$zone)
   x$zone[x$zone %in% c("2F", "12F")] <- "F"
   x$zone[x$zone %in% c("2E", "12E")] <- "E"
   
   # Clean-up missing leg coding:
   x$missing.legs <- toupper(x$missing.legs)
   x$missing.legs <- gsub("M", "1", x$missing.legs)
   x$missing.legs <- gsub("R", "2", x$missing.legs)
   
   # Remove leading and trailing spaces:
   x$vessel <- gsub("(^[ ]+)|([ ]+$)", "", x$vessel)
   x$trap.number <- as.numeric(x$trap.number)
   x$trip.number <- toupper(x$trip.number)
   
   # Remove blank records:
   ix <- is.na(x$sex) & is.na(x$carapace.width) & is.na(x$abdomen.width) &
         is.na(x$chela.height.right) & is.na(x$chela.height.left) & is.na(x$shell.condition)
   x <- x[!ix, ]
   
   # Remove empty columns:
   x$fishing.grid <- gsub(" ", "", x$fishing.grid)
   ix <- which(unlist(lapply(x, function(x) return(all(is.na(x) | all(x == "") | all(x == "*"))))))
   x <- x[, -ix]
   
   # Add trap number if it was removed:
   if (!("trap.number" %in% names(x))) x$trap.number <- 1
   
   # Add species code:
   x$species <- 2526   
   
   # Remove duplicate records:
   ix <- !duplicated(x[c("trip.number", "crab.number", "trap.number", "carapace.width")])
   if (sum(!ix) > 0){
      cat(paste0("Removed ", sum(!ix) ," duplicate records from ", length(unique(x$trip.number[!ix])), " trip(s).\n"))
      x <- x[ix, ]
   }

   return(x)
}
