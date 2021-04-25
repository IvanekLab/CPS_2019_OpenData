;; Created (2019) by Claire Zoellner, Renata Ivanek, Martin Wiedmann at Cornell University
;; Adapted (2019) by Genevieve Sullivan, Claire Zoellner, Renata Ivanek, Martin Wiedmann at Cornell University
;; Funding provided by American Frozen Food Institute (AFFI) and Center for Produce Safety (CPS)

;; An agent-based model that simulates Listeria spp. behavior in a food production facility and control strategies of cleaning and environmental sampling
;; each tick represents one hour

extensions [csv matrix]

globals ;define global variables that can be used and referenced anywhere in the code. Note sliders and choosers on the interface are already global so do not need to be defined here
[ production-schedule sunday monday tuesday wednesday thursday friday saturday floor-wall-juncture
  event day week
  empty-time clean-time production-time
  random-week
  zones-initial z1-initial z2-initial z3-initial p-initial e-initial
  zones-pre-prev m-zones-pre t-zones-pre w-zones-pre r-zones-pre f-zones-pre
  zones-beg-prev m-zones-beg t-zones-beg w-zones-beg r-zones-beg f-zones-beg
  zones-mid-prev m-zones-mid t-zones-mid w-zones-mid r-zones-mid f-zones-mid
  zones-post-prev m-zones-post t-zones-post w-zones-post r-zones-post f-zones-post
  z3-pre  m-z3-pre t-z3-pre w-z3-pre r-z3-pre f-z3-pre
  z1-pre  m-z1-pre t-z1-pre w-z1-pre r-z1-pre f-z1-pre
  z2-pre m-z2-pre t-z2-pre w-z2-pre r-z2-pre f-z2-pre
  p-pre  m-p-pre t-p-pre w-p-pre r-p-pre f-p-pre
  e-pre  m-e-pre t-e-pre w-e-pre r-e-pre f-e-pre
  z3-beg  m-z3-beg t-z3-beg w-z3-beg r-z3-beg f-z3-beg
  z1-beg  m-z1-beg t-z1-beg w-z1-beg r-z1-beg f-z1-beg
  z2-beg m-z2-beg t-z2-beg w-z2-beg r-z2-beg f-z2-beg
  p-beg  m-p-beg t-p-beg w-p-beg r-p-beg f-p-beg
  e-beg  m-e-beg t-e-beg w-e-beg r-e-beg f-e-beg
  z3-mid  m-z3-mid t-z3-mid w-z3-mid r-z3-mid f-z3-mid
  z1-mid  m-z1-mid t-z1-mid w-z1-mid r-z1-mid f-z1-mid
  z2-mid m-z2-mid t-z2-mid w-z2-mid r-z2-mid f-z2-mid
  p-mid  m-p-mid t-p-mid w-p-mid r-p-mid f-p-mid
  e-mid  m-e-mid t-e-mid w-e-mid r-e-mid f-e-mid
  z3-post  m-z3-post t-z3-post w-z3-post r-z3-post f-z3-post
  z1-post  m-z1-post t-z1-post w-z1-post r-z1-post f-z1-post
  z2-post m-z2-post t-z2-post w-z2-post r-z2-post f-z2-post
  p-post  m-p-post t-p-post w-p-post r-p-post f-p-post
  e-post  m-e-post t-e-post w-e-post r-e-post f-e-post
  z1-conc-pre m-z1-conc-pre t-z1-conc-pre w-z1-conc-pre r-z1-conc-pre f-z1-conc-pre
  z2-conc-pre m-z2-conc-pre t-z2-conc-pre w-z2-conc-pre r-z2-conc-pre f-z2-conc-pre
  z3-conc-pre m-z3-conc-pre t-z3-conc-pre w-z3-conc-pre r-z3-conc-pre f-z3-conc-pre
  p-conc-pre m-p-conc-pre t-p-conc-pre w-p-conc-pre r-p-conc-pre f-p-conc-pre
  e-conc-pre m-e-conc-pre t-e-conc-pre w-e-conc-pre r-e-conc-pre f-e-conc-pre
  z1-conc-mid m-z1-conc-mid t-z1-conc-mid w-z1-conc-mid r-z1-conc-mid f-z1-conc-mid
  z2-conc-mid m-z2-conc-mid t-z2-conc-mid w-z2-conc-mid r-z2-conc-mid f-z2-conc-mid
  z3-conc-mid m-z3-conc-mid t-z3-conc-mid w-z3-conc-mid r-z3-conc-mid f-z3-conc-mid
  p-conc-mid m-p-conc-mid t-p-conc-mid w-p-conc-mid r-p-conc-mid f-p-conc-mid
  e-conc-mid m-e-conc-mid t-e-conc-mid w-e-conc-mid r-e-conc-mid f-e-conc-mid
  z1-conc-post m-z1-conc-post t-z1-conc-post w-z1-conc-post r-z1-conc-post f-z1-conc-post
  z2-conc-post m-z2-conc-post t-z2-conc-post w-z2-conc-post r-z2-conc-post f-z2-conc-post
  z3-conc-post m-z3-conc-post t-z3-conc-post w-z3-conc-post r-z3-conc-post f-z3-conc-post
  p-conc-post m-p-conc-post t-p-conc-post w-p-conc-post r-p-conc-post f-p-conc-post
  e-conc-post m-e-conc-post t-e-conc-post w-e-conc-post r-e-conc-post f-e-conc-post
  zone-to-zone avg-listeria-transferred zone-to-patch condensation-drip introduction-zone4 patch-to-patch
  introduction-food sum-fruit-load sum-fruit-transfer random-event employee-FCS NFCS-FCS
  total-contam-events
  m-random-event t-random-event w-random-event r-random-event f-random-event
  collect-data?
  median-z1-time median-z2-time median-z3-time median-e-time median-p-time median-fwj-time
  median-z1-max-contam-bout median-z2-max-contam-bout median-z3-max-contam-bout median-e-max-contam-bout median-p-max-contam-bout median-fwj-max-contam-bout
  temporary-niches p-temp-niches z1-temp-niches z2-temp-niches z3-temp-niches
  z1-contacts z1-transfers z2-contacts z2-transfers z3-contacts z3-transfers
  concentration-list time-contaminated-list max-contam-bout-list contacts-list transfers-list temp-niche-list niches-estab-list
  fruit-cont-list zone4-cont-list random-cont-list employees

  ;;input parameter globals
  spread-tc randomEvent random-event-list
  m-crate-prevalence t-crate-prevalence r-crate-prevalence f-crate-prevalence fruit-conc w-crate-prevalence fruit-tc zone4-load random-load
  m-san-reduction  t-san-reduction  w-san-reduction  r-san-reduction  f-san-reduction prob-cleanable
  mu-max K p-patch-spread p-water-spread crate-prevalence-min crate-prevalence-ml crate-prevalence-max p-food-dropped p-condensation p-random-noise p-zone4-intro
  traffic-high traffic-low traffic-neg
  hoelzer-tc-matrix hoelzer-std-matrix expert-link-min-matrix expert-link-ml-matrix expert-link-max-matrix
  p-transfer tc flow-rate crate-size
  p11 p12 p13 p14 p21 p22 p23 p24 p31 p32 p33 p34 p41 p42 p43 p44
  tc11 tc12 tc13 tc14 tc21 tc22 tc23 tc24 tc31 tc32 tc33 tc34 tc41 tc42 tc43 tc44

  ;;sampling function and output globals
  a-table b-table c-table d-table a-table2 b-table2 c-table2 d-table2
  start-sampling? sample-day sample-week sample-results sample-time totalSamples takeSamples samples-taken samples-to-take-now time-to-sampling-results possible-z-sites possible-z-patches possible-p-sites zone1-sample zone2-sample zone3-sample
  z1-samples z2-samples z3-samples floors-samples floor-wall-samples
  z-sample-prev m-z-sample-prev t-z-sample-prev w-z-sample-prev r-z-sample-prev f-z-sample-prev
  daily-samples daily-prev
  z-sample-beg-prev z-sample-mid-prev z-sample-post-prev
  beg-samples mid-samples post-samples
  belt-prev cart-prev control-panel-prev door-prev drain-prev employee-prev floors-prev flume-prev frame-prev frame-in-prev ladder-prev misc-prev packing-prev squeegee-prev
  trash-gray-prev trash-white-prev trash-yellow-prev wall-prev weigher-prev a-b-prev a-flume-prev a-line-prev a-under-room-prev b-flume-prev b-line-prev coring-prev
  elevator-prev hopper-prev packing--prev scoring-prev sorting-prev spin-dryer-prev weigher-platform-prev A-prev B-prev

  ;;validation outcomes - set these up based on the facility's historical EM data
  belt-samples cart-samples control-panel-samples door-samples drain-samples employee-samples floor-samples flume-samples frame-samples frame-in-samples ladder-samples
  misc-samples packing-samples squeegee-samples trash-gray-samples trash-white-samples trash-yellow-samples wall-samples weigher-samples a-b-samples a-flume-samples
  a-line-samples a-under-room-samples b-flume-samples b-line-samples coring-samples elevator-samples hopper-samples packing--samples scoring-samples sorting-samples
  spin-dryer-samples weigher-platform-samples A-samples B-samples

  z1-prev z2-prev z3-prev
  random-floor random-ceiling random-agent

  ]

breed [ zones a-zone] ; node agents classified as zone objects
;breed [ patches a-patch]
breed [ listeria a-listeria] ; listeria agents

;; every link breed must be declared as either directed or undirected
directed-link-breed [contact-links contact-link]
undirected-link-breed [proximity-links proximity-link]

zones-own
[ z-category ; zone classification of either 1, 2, 3, 4
  z-item-name ; what is the site name
  z-height ;distance from the floor
  z-area ; surface area in sq. cm
  z-cleanable?  ;true/false to describe if item needs to be dissassembled when cleaned in order to remove listeria
  z-equipment ; if processing room is divided into "areas" relevant to activity, can assign site to equipment (i.e. hand-packing)
  z-location ; added by GBS to specify location of different agents
  z-room ; added by GBS to specify room
  z-water ; water/condensation present at that zone [dry, moist, visibly wet]
  z-listeria ;number of listeria present at zone CFU
  z-listeria-concentration ; concentration of listeria present at that zone/z-area of zone CFU/sq. cm
  ;z-sanitation-schedule ;not used currently but can be used to set frequency of dissassembly
  z-sampled?
  z-resample?
  z-time-contaminated ; keeps track of total time spent contaminated (counts ticks when z-listeria != 0)
  z-max-contam-bout ; keeps track of consecutive time contaminated
  z-contam-counter
  z-contacts ;number of times a zone is contacted by other zone with contamination
  z-transfers ;number of times a zone transfers contamination to adjacent zone
  z-fruit ; number of times contaminated from incoming product
  z-zone4 ; number of time contaminated from objects coming into the room
  z-random ; number of times contaminated from "random" event
  z-temp-niche ; number of times cleaning not executed, leaving a site contaminated
  z-time-to-niche ; number of ticks until site is not cleaned properly
  z-niche-length ; consecutive ticks an uncleanable site remains contaminated
  z-niches-established  ;number of times an uncleanable site becomes contaminated
  z-out-links ; counts out-directed links
  z-in-links ; counts in-directed links
  z-undirected-links ; counts undirected links
  ]

patches-own
[ p-water           ; the amount of water/sanitizer on this patch
  p-listeria ;number of listeria present at patch
  p-listeria-concentration  ; concentration of listeria present at patch/625cm^2
  p-traffic ;assigned traffic level from traffic map
  p-area   ;area of each patch = 25 x 25 = 625 cm^2
  p-num-positives
  p-cleanable?
  p-sampled?
  p-sanitation-schedule
  p-item-name
  p-time-contaminated
  p-max-contam-bout
  p-contam-counter
  p-zone4
  p-random
  p-temp-niche
  c-height           ;height of ceiling (20 ft)
  c-listeria  ;listeria population cfu
  c-listeria-concentration  ;listeria concentration (cfu/cm^2) on ceiling
  c-water          ;ceiling water level
  c-time-contaminated
  c-max-contam-bout
  c-contam-counter
  c-random
  m-height           ;height of mezzanine (11 ft)
  m-floor-listeria  ;listeria population on mezz floor
  m-floor-listeria-concentration  ;listeria concentration (cfu/cm^2) on mezzanine floor
  m-listeria-concentration
  m-ceiling-listeria ; listeria on underside of mezzanine
  m-ceiling-listeria-concentration ;listeria concentration (cfu/cm^2) on mezzanine underside
  m-water          ;mezz floor water
  m-traffic
  m-time-contaminated
  m-max-contam-bout
  m-contam-counter
  m-random
]

to setup
  clear-all
  file-close-all
  with-local-randomness[ ;we interrupt the global seed during setup so that we get different combinations of parameter values when sampling the parameter space
    random-seed seed-me ;this sets a new seed each iteration, even if we set a global seed from the beginning in BSpace - we do not know this seed value
  import-network
  ask zones [zones-setup]
  setup-environment
  import-matrices
  setup-submodel-params
  ask links [ set color black set thickness 0.2 ]
  create-listeria round ((initial-zone-prevalence / 100) * count zones)
  [ listeria-setup
    listeria-initial-zone-placement ]
  create-listeria round ((initial-environment-prevalence / 100) * count (patches with [p-water != 0 or p-item-name = "floor-wall-juncture"]))
  [ listeria-setup
    listeria-initial-environment-placement ]
  setup-weekly-schedule
  ;set-sampling ; uncomment if testing sampling
  ]
  if p-random-noise != 0 [set random-event-list [ ] setup-randomEvents]
  create-output-files
  reset-ticks
  ask turtle 259 [die] ;remove floor agents
  ask turtle 266 [die]
  ask turtle 268 [die]
  ask turtle 269 [die]
  ask turtle 270 [die]
  ask turtle 272 [die]
  ask turtle 278 [die]
  ask turtle 279 [die]
  ask turtle 280 [die]
  ask turtle 282 [die]
  ask turtle 283 [die]
  ask turtle 284 [die]
  ask turtle 286 [die]
  ask turtle 291 [die]

end

to setup-randomEvents ; make a list of ticks at which random events will occur
  repeat 100 [
    set random-event-list lput (round (random-exponential (1 / p-random-noise))) random-event-list] ;expert opinion (E=0.01) print (word "random event:" randomEvent)
  ;print random-event-list
  set randomEvent ticks + item 0 random-event-list
  ;print randomEvent
  set random-event-list but-first random-event-list
  ;print random-event-list
end

to create-output-files ; these output files are written to collect data when running a single iteration - mainly for testing purposes, they are not utilized in simulation experiments
  if (file-exists? "TimeSeriesZoneConcen.csv" ) ;;will refer to as File 1
    [ carefully [file-delete "TimeSeriesZoneConcen.csv" ] [ ] ]
 if (file-exists? "ZoneData.csv" ) ;;will refer to as File 2
    [ carefully [file-delete "ZoneData.csv" ] [ ] ]

  let headings-file-1 [ "hours" ] ;making col header of file "TimeSeriesZoneConcen.csv"
  let zone-item [ "," ]
  foreach sort-on [who] zones ; builds the other col headings by using who number
      [ [?1] -> ask ?1 [ set headings-file-1 lput who headings-file-1
                set zone-item lput z-item-name zone-item ] ]
    ;print headings-file-1

  file-open "TimeSeriesZoneConcen.csv"
  let headings1-converted csv:to-row headings-file-1
  let zones-converted csv:to-row zone-item
  file-print  headings1-converted
  file-print zones-converted
  file-close

  let headings-file-2 [ "who" "zone-item" "zone-category" "area" "cleanable" "out-links" "in-links" "undirected-links" ] ;making col header of file "ZoneData.csv"
  file-open "ZoneData.csv"
  let headings2-converted csv:to-row headings-file-2
  file-print  headings2-converted
  file-close
end

to import-matrices ; transfer probability and coefficient matrices
  let expert-min (csv:from-file "expert-transfer-min-matrix_produce.csv")
  set expert-link-min-matrix matrix:from-row-list expert-min

  let expert-ml (csv:from-file "expert-transfer-ml-matrix_produce.csv")
  set expert-link-ml-matrix matrix:from-row-list expert-ml

  let expert-max (csv:from-file "expert-transfer-max-matrix_produce.csv")
  set expert-link-max-matrix matrix:from-row-list expert-max

  let hoelzer-tc (csv:from-file "transfer-tc-matrix.csv")
  set hoelzer-tc-matrix matrix:from-row-list hoelzer-tc

  let hoelzer-std (csv:from-file "transfer-std-matrix.csv")
  set hoelzer-std-matrix matrix:from-row-list hoelzer-std
end

to setup-submodel-params ; here are all of the global parameter values that need to be determined for each product/facility - some have to be assumed or used from expert opinion
   let GT randomfloat-in-range 47 155 ; generation time (hr) fruit mix @ 4-8C Ziegler et al., 2018 doi: https://doi.org/10.4081/ijfs.2018.7581
   set mu-max ((ln 2) / GT) ; don't modify this
   set K (1e5) ;max growth level/carrying capacity (CFU) at <5C FDA/FSIS Lm RA (2003), can use uniform (1e4, 1e6) to represent uncertainty
   set p-patch-spread random-pert 0.03 0.25 0.65 4 ; prob spread on floor via traffic from Chambers et al., 2009 - don't need to change this
   set traffic-high 60 ; contact rate/hr on high traffic patches - observed and then calculated, but could also be a default
   set traffic-low 12 ; contact rate/hr on low traffic patches - observed and then calculated, but could also be a default
   set traffic-neg 0.2 ; contact rate/hr on negligible traffic patches - observed and then calculated, but could also be a default
   set p-water-spread randomfloat-in-range 0.01 0.05 ; prob spread on floor via water - assumed - don't need to change this
   set m-san-reduction (10 ^ (random-pert -8 -6 -1.5 4)) ; log reduction from daily cleaning and sanitation - don't need to change this
   set t-san-reduction (10 ^ (random-pert -8 -6 -1.5 4)) ; Mokhtari et al, 2018 used sanitation efficiency scenarios of 1, 2, 3, and 4 log10 ??
   set w-san-reduction (10 ^ (random-pert -8 -6 -1.5 4)) ; Zoellner et al, 2019 used sanitation reduction pert (-8, -6, -1.5, 4) from FDA/FSIS, 2013
   set r-san-reduction (10 ^ (random-pert -8 -6 -1.5 4))
   set f-san-reduction (10 ^ (random-pert -8 -6 -1.5 4))
   set p-food-dropped randomfloat-in-range 5 10 ; probability food or tools in zone 1 drop to floor - should be observed - can use this as default (or smoked salmon as 20-40)
   set p-condensation randomfloat-in-range 0 2 ; prob condensation falls - assumed - don't need to change unless no condensation or condensation is big issue
   set p-random-noise (10 ^ (random-pert -4.3 -0.9 -0.6 4.6)) ; probability of random event - taken from expert elicitation - could repeat or use these values ; GBS changed pr first shift
   set p-zone4-intro  (10 ^ (random-pert -2.3 -0.9 -0.2 4.8)) * 100 ; prob of intro from objects in zone4 - taken from expert elicitation - could repeat or use these values ; GBS changed to pz first shift
   ;set p-zone4-intro  (10 ^ (random-pert -4.3 -2.3 -0.6 4.8)) * 100 ; prob of intro from objects in zone4 - taken from expert elicitation - could repeat or use these values ; GBS changed to pz pre-op
   set spread-tc randomfloat-in-range 0 5 ; transfer coefficient for patch spread events - assumed - don't need to change
   set prob-cleanable random-in-range 95 100 ;probability that cleanable sites are truly cleanable - assumed
   set zone4-load 10 ^ random-pert 0 1.9 3.3 4.2 ; number of listeria CFU introduced during zone4-intro - expert elicitation ; GBS changed Nz
   set random-load 10 ^ random-pert 0.2 3.3 3.7 3.3 ; number of listeria CFU introduced during random-intro - expert elicitation ; GBS changed to Nr
   set crate-prevalence-min -2.3 ; min prevalence of L spp. in incoming raw materials - expert opinion - could use unless literature value or company has product testing data ; GBS changed to Rd
   set crate-prevalence-ml -0.6 ; most likely prevalence of L spp. in incoming raw materials - expert opinion - could use unless literature value or company has product testing data ; GBS changed to Rd
   set crate-prevalence-max -0.6 ; max prevalence of L spp. in incoming raw materials - Lm prev -2.22 from Johannessen et al., 2002 - strawberries in Norway ; GBS changed to Rd
   set m-crate-prevalence 10 ^ (random-pert crate-prevalence-min crate-prevalence-ml crate-prevalence-max 5.4) ; sets daily prev in raw materials - don't change
   set t-crate-prevalence 10 ^ (random-pert crate-prevalence-min crate-prevalence-ml crate-prevalence-max 5.4) ; sets daily prev in raw materials - don't change
   set w-crate-prevalence 10 ^ (random-pert crate-prevalence-min crate-prevalence-ml crate-prevalence-max 5.4) ; sets daily prev in raw materials - don't change
   set r-crate-prevalence 10 ^ (random-pert crate-prevalence-min crate-prevalence-ml crate-prevalence-max 5.4) ; sets daily prev in raw materials - don't change
   set f-crate-prevalence 10 ^ (random-pert crate-prevalence-min crate-prevalence-ml crate-prevalence-max 5.4) ; sets daily prev in raw materials - don't change
   let fruit-conc-mean 0.1 ; mean cfu/g L spp. in raw material - look in literature for appropriate value
   let fruit-conc-var 5.2 ; variance cfu/g L spp. in raw material - look in literature - var = sd^2 or range/4
   ;let fruit-conc-mean 255 ; mean cfu/g L spp. in raw material - GBS changed based on expert elicitation
   ;let fruit-conc-var 628 ; variance cfu/g L spp. in raw material - GBS changed based on expert elicitation - var = range/4 = ((10^3.4)-(10^-0.1))/4
   let fruit-conc-alpha (fruit-conc-mean ^ 2) / fruit-conc-var
   let fruit-conc-theta fruit-conc-mean / fruit-conc-var  ; calculates parameter for gamma dist - don't change
   set fruit-conc (random-gamma fruit-conc-alpha fruit-conc-theta)  ; samples CFU/g for level of L spp in contamin incoming product - don't change
   ;set fruit-conc (10 ^ (random-pert -0.1 2.2 3.4 4)) ; concentration of Listeria spp per contaminated leafy green or apple (CFU/g) changed by GBS to NR
   set fruit-tc min (list 1 (10 ^ random-normal (-0.28) (0.2))) ; transfer coefficient from product to surfaces (meat to stainless in Hoelzer 2012) - could look for better value depending on product
   set flow-rate 4.4 ; flow rate/hr of incoming materials - need to ask this from facility
   set crate-size 453592 ; grams of product introduced per unit of flow rate - need to ask this from facility
   set collect-data? false
   reset-lists
   setup-transfer-params
end

to reset-lists
   set concentration-list [] set time-contaminated-list [] set max-contam-bout-list [] set contacts-list [] set transfers-list []
   set temp-niche-list [] set niches-estab-list [] set fruit-cont-list [] set zone4-cont-list [] set random-cont-list []
end

to setup-transfer-params ; reads in transfer matrices from .csv files - don't change these
  set p-transfer [ ]
  set tc [ ]
  let min11  (matrix:get expert-link-min-matrix 0 0) ;expert elicitation probabilities
  let ml11    (matrix:get expert-link-ml-matrix 0 0)
  let max11  (matrix:get expert-link-max-matrix 0 0)
  let min12  (matrix:get expert-link-min-matrix 0 1)
  let ml12   (matrix:get expert-link-ml-matrix 0 1)
  let max12 (matrix:get expert-link-max-matrix 0 1)
  let min13  (matrix:get expert-link-min-matrix 0 2)
  let ml13   (matrix:get expert-link-ml-matrix 0 2)
  let max13 (matrix:get expert-link-max-matrix 0 2)
  let min14  (matrix:get expert-link-min-matrix 0 3)
  let ml14   (matrix:get expert-link-ml-matrix 0 3)
  let max14 (matrix:get expert-link-max-matrix 0 3)

  let min21  (matrix:get expert-link-min-matrix 1 0) ;expert elicitation probabilities
  let ml21    (matrix:get expert-link-ml-matrix 1 0)
  let max21  (matrix:get expert-link-max-matrix 1 0)
  let min22  (matrix:get expert-link-min-matrix 1 1)
  let ml22   (matrix:get expert-link-ml-matrix 1 1)
  let max22 (matrix:get expert-link-max-matrix 1 1)
  let min23  (matrix:get expert-link-min-matrix 1 2)
  let ml23   (matrix:get expert-link-ml-matrix 1 2)
  let max23 (matrix:get expert-link-max-matrix 1 2)
  let min24  (matrix:get expert-link-min-matrix 1 3)
  let ml24   (matrix:get expert-link-ml-matrix 1 3)
  let max24 (matrix:get expert-link-max-matrix 1 3)

  let min31  (matrix:get expert-link-min-matrix 2 0) ;assumed probabilities
  let ml31    (matrix:get expert-link-ml-matrix 2 0)
  let max31  (matrix:get expert-link-max-matrix 2 0)
  let min32  (matrix:get expert-link-min-matrix 2 1)
  let ml32   (matrix:get expert-link-ml-matrix 2 1)
  let max32 (matrix:get expert-link-max-matrix 2 1)
  let min33  (matrix:get expert-link-min-matrix 2 2)
  let ml33   (matrix:get expert-link-ml-matrix 2 2)
  let max33 (matrix:get expert-link-max-matrix 2 2)
  let min34  (matrix:get expert-link-min-matrix 2 3)
  let ml34   (matrix:get expert-link-ml-matrix 2 3)
  let max34 (matrix:get expert-link-max-matrix 2 3)

  let min41  (matrix:get expert-link-min-matrix 3 0) ;expert elicitation probabilities
  let ml41    (matrix:get expert-link-ml-matrix 3 0)
  let max41  (matrix:get expert-link-max-matrix 3 0)
  let min42  (matrix:get expert-link-min-matrix 3 1)
  let ml42   (matrix:get expert-link-ml-matrix 3 1)
  let max42 (matrix:get expert-link-max-matrix 3 1)
  let min43  (matrix:get expert-link-min-matrix 3 2)
  let ml43   (matrix:get expert-link-ml-matrix 3 2)
  let max43 (matrix:get expert-link-max-matrix 3 2)
  let min44  (matrix:get expert-link-min-matrix 3 3)
  let ml44   (matrix:get expert-link-ml-matrix 3 3)
  let max44 (matrix:get expert-link-max-matrix 3 3)

  set p-transfer  lput (list random-pert min11 ml11 max11 4 random-pert min12 ml12 max12 4 random-pert min13 ml13 max13 4 random-pert min14 ml14 max14 4) p-transfer
  set p-transfer  lput (list random-pert min21 ml21 max21 4 random-pert min22 ml22 max22 4 random-pert min23 ml23 max23 4 random-pert min24 ml24 max24 4) p-transfer
  set p-transfer  lput (list random-pert min31 ml31 max31 4 random-pert min32 ml32 max32 4 random-pert min33 ml33 max33 4 random-pert min34 ml34 max34 4) p-transfer
  set p-transfer  lput (list random-pert min41 ml41 max41 4 random-pert min42 ml42 max42 4 random-pert min43 ml43 max43 4 random-pert min44 ml44 max44 4) p-transfer

  set p11 (item 0 (item 0 p-transfer)) set p12 (item 1 (item 0 p-transfer)) set p13 (item 2 (item 0 p-transfer)) set p14 (item 3 (item 0 p-transfer))
  set p21 (item 0 (item 1 p-transfer)) set p22 (item 1 (item 1 p-transfer)) set p23 (item 2 (item 1 p-transfer)) set p24 (item 3 (item 1 p-transfer))
  set p31 (item 0 (item 2 p-transfer)) set p32 (item 1 (item 2 p-transfer)) set p33 (item 2 (item 2 p-transfer)) set p34 (item 3 (item 2 p-transfer))
  set p41 (item 0 (item 3 p-transfer)) set p42 (item 1 (item 3 p-transfer)) set p43 (item 2 (item 3 p-transfer)) set p44 (item 3 (item 3 p-transfer))

  let mu11  (matrix:get hoelzer-tc-matrix 0 0) ; tc estimated from Hoelzer et al., 2012a
  let std11 (matrix:get hoelzer-std-matrix 0 0) ;std estimated from Hoelzer et al., 2012a
  let mu12  (matrix:get hoelzer-tc-matrix 0 1)
  let std12 (matrix:get hoelzer-std-matrix 0 1)
  let mu13  (matrix:get hoelzer-tc-matrix 0 2)
  let std13 (matrix:get hoelzer-std-matrix 0 2)
  let mu14  (matrix:get hoelzer-tc-matrix 0 3)
  let std14 (matrix:get hoelzer-std-matrix 0 3)

  let mu21  (matrix:get hoelzer-tc-matrix 1 0) ; tc estimated from Hoelzer et al., 2012a
  let std21 (matrix:get hoelzer-std-matrix 1 0) ;std estimated from Hoelzer et al., 2012a
  let mu22  (matrix:get hoelzer-tc-matrix 1 1)
  let std22 (matrix:get hoelzer-std-matrix 1 1)
  let mu23  (matrix:get hoelzer-tc-matrix 1 2)
  let std23 (matrix:get hoelzer-std-matrix 1 2)
  let mu24  (matrix:get hoelzer-tc-matrix 1 3)
  let std24 (matrix:get hoelzer-std-matrix 1 3)

  let mu31  (matrix:get hoelzer-tc-matrix 2 0) ; tc estimated from Hoelzer et al., 2012a
  let std31 (matrix:get hoelzer-std-matrix 2 0) ;std estimated from Hoelzer et al., 2012a
  let mu32  (matrix:get hoelzer-tc-matrix 2 1)
  let std32 (matrix:get hoelzer-std-matrix 2 1)
  let mu33  (matrix:get hoelzer-tc-matrix 2 2)
  let std33 (matrix:get hoelzer-std-matrix 2 2)
  let mu34  (matrix:get hoelzer-tc-matrix 2 3)
  let std34 (matrix:get hoelzer-std-matrix 2 3)

  let mu41  (matrix:get hoelzer-tc-matrix 3 0) ; tc estimated from Hoelzer et al., 2012a
  let std41 (matrix:get hoelzer-std-matrix 3 0) ;std estimated from Hoelzer et al., 2012a
  let mu42  (matrix:get hoelzer-tc-matrix 3 1)
  let std42 (matrix:get hoelzer-std-matrix 3 1)
  let mu43  (matrix:get hoelzer-tc-matrix 3 2)
  let std43 (matrix:get hoelzer-std-matrix 3 2)
  let mu44  (matrix:get hoelzer-tc-matrix 3 3)
  let std44 (matrix:get hoelzer-std-matrix 3 3)

  set tc lput (list random-normal mu11 std11 random-normal mu12 std12 random-normal mu13 std13 random-normal mu14 std14 ) tc
  set tc lput (list random-normal mu21 std21 random-normal mu22 std22 random-normal mu23 std23 random-normal mu24 std24 ) tc
  set tc lput (list random-normal mu31 std31 random-normal mu32 std32 random-normal mu33 std33 random-normal mu34 std34 ) tc
  set tc lput (list random-normal mu41 std41 random-normal mu42 std42 random-normal mu43 std43 random-normal mu44 std44 ) tc

  set tc11 (item 0 (item 0 tc)) set tc12 (item 1 (item 0 tc)) set tc13 (item 2 (item 0 tc)) set tc14 (item 3 (item 0 tc))
  set tc21 (item 0 (item 1 tc)) set tc22 (item 1 (item 1 tc)) set tc23 (item 2 (item 1 tc)) set tc24 (item 3 (item 1 tc))
  set tc31 (item 0 (item 2 tc)) set tc32 (item 1 (item 2 tc)) set tc33 (item 2 (item 2 tc)) set tc34 (item 3 (item 2 tc))
  set tc41 (item 0 (item 3 tc)) set tc42 (item 1 (item 3 tc)) set tc43 (item 2 (item 3 tc)) set tc44 (item 3 (item 3 tc))
end

to import-network
  set-default-shape turtles "circle"
  import-attributes
  import-directed-links
  import-undirected-links
  reset-ticks
end

to import-attributes ; this creates agents with their individual characteristics that you have from a sample site list
  file-close-all
  file-open "agents_031720.txt" ; save excel file without headers as .txt - file needs to be in same folder as netlogo file
  ; each row of file is an agent
  ; each column contains agent attributes in this order:
  ; item_name xcor ycor z-category z-height z-area z-cleanable? z-equipment
let i 1
  while [i <= 294] ; change this number to match the number of rows in your excel file
 ; while [not file-at-end?]
  [
    let items split file-read-line "\t"
    let itemsA (list
      item 0 items
      read-from-string item 1 items
      read-from-string item 2 items
      read-from-string item 3 items
      read-from-string item 4 items
      read-from-string item 5 items
      read-from-string item 6 items
      item 7 items
      item 8 items
      item 9 items)
    create-zones 1 [
      set z-item-name    item 0 itemsA
      set xcor    item 1 itemsA
      set ycor    item 2 itemsA
      set z-category    item 3 itemsA
      set z-height item 4 itemsA
      set z-area item 5 itemsA
      ifelse (item 6 itemsA) = 1 [set z-cleanable? true] [set z-cleanable? false]
      set z-equipment item 7 itemsA
      set z-location item 8 itemsA
      set z-room item 9 itemsA
    ]
  set i (i + 1) ;remove if using not file-at-end?
  ]
  file-close
end

to import-directed-links ; sets up directed links
  file-close-all
  ; create excel file with 2 columns of agent who numbers, directed link will be made from who in col 1 to who in col 2
  file-open "dlinks_new.txt" ; save excel file as .txt and use name here
  while [not file-at-end?]
  [let items read-from-string (word "[" file-read-line "]")
    ask get-node (item 0 items)
    [create-contact-link-to get-node (item 1 items)]  ;directed links
  ]
  file-close
end

to import-undirected-links ; sets up undirected links
  file-close-all
  ; create excel file with 2 columns of agent who numbers, link will be made between who in col 1 and who in col 2
  file-open "links_new.txt" ; save excel file as .txt and use name here
  while [not file-at-end?]
  [ let items read-from-string (word "[" file-read-line "]")
    ask get-node (item 0 items)
    [create-proximity-link-with get-node (item 1 items)] ;undirected links
  ]
  file-close
end

to-report split [ string delim ]
  report reduce [ [?1 ?2] ->
    ifelse-value (?2 = delim)
      [ lput "" ?1 ]
      [ lput word last ?1 ?2 but-last ?1 ]
  ] fput [""] n-values (length string) [ [?1] -> substring string ?1 (?1 + 1) ]
end

to-report get-node [id]
  report one-of turtles with [who = id]
end

to zones-setup ; sets shape and color of agents depending on zone category
  set z-listeria-concentration 0
  set z-listeria 0
  set z-water one-of [ 1 2 3]
  if z-category = 1
    [ set shape "circle" set size 1 set color (16 - z-water)]
  if z-category = 2
    [ set shape "triangle" set size 1 set color (26 - z-water)]
  if z-category = 3
    [ set shape "pentagon" set size 1 set color (126 - z-water)]
  if z-item-name = "employee"
    [ set shape "person" set size 3 set color black set hidden? true]
  set z-resample? false
  set z-sampled? false
  set z-out-links count my-out-contact-links
  set z-in-links count my-in-contact-links
  set z-undirected-links count my-proximity-links
  set employees (zones with [ z-item-name = "employee" ])
end

to setup-environment ; set inital conditions and characteristics of patches
  update-water ( "no" )
  update-traffic ( "no" )
  run environment-view
  ask patches with [p-water != 0]
    [set p-listeria-concentration 0 set p-listeria 0
    set c-listeria-concentration 0  set c-listeria 0
    set c-height 20
    set p-area 900
    set p-cleanable? true
   ]
  ask patches with [p-water = 35]
     [set pcolor brown
       set p-cleanable? false
       set p-item-name "floor-wall-juncture"]
  ask patches with [(pxcor > 56 and pxcor < 83) and (pycor > 27 and pycor < 38)] ;;check where is mezzanine
    [ set m-height 8          ;height of mezzanine (8 ft)
      set m-traffic 1
      set m-floor-listeria 0 set m-floor-listeria-concentration 0
      set m-ceiling-listeria 0 set m-ceiling-listeria-concentration 0 ]
end

to update-water [level] ; reads in different water maps (if applicable) - change name of file and add other levels if needed (i.e., "high")

  if level = "high" [
    file-open "water-high_030320.txt"
  while [not file-at-end?] [
    foreach sort patches [ [?1] ->
  ask ?1 [
   set p-water file-read]
  ]]
  file-close
  ask patches with [m-height > 0]
    [set m-water one-of [1 2]]
  ]

  if level = "medium" [
    file-open "water_030320.txt"
  while [not file-at-end?] [
    foreach sort patches [ [?1] ->
  ask ?1 [
   set p-water file-read]
  ]]
  file-close
  ask patches with [m-height > 0]
    [set m-water one-of [1]]
  ]

 if level = "no" [
  file-open "A_Floorplan_key_030320.txt"
   ;read in the floor plan
  while [not file-at-end?] [
    foreach sort patches [ [?1] ->
  ask ?1 [
   set p-water file-read]
  ]]
  file-close
  ]
 water
end

;;environment-view switch requires the two following functions:
to water
  ask patches[ p-water-recolor]
end

to traffic
  ask patches [p-traffic-recolor]
end

to update-traffic [level] ; reads in different traffic maps (if applicable) - change name of file and add other levels if needed (i.e., "high")
  if level = "high" [
    file-open "traffic-high_030320.txt"
    while [not file-at-end?] [
    foreach sort patches [ [?1] ->
  ask ?1 [
   set p-traffic file-read]
  ]]
  file-close
  ask patches with [m-height > 0]
  [set m-traffic 1]
  ]

  if level = "low" [
    file-open "traffic-low_030320.txt"
    while [not file-at-end?] [
    foreach sort patches [ [?1] ->
  ask ?1 [
   set p-traffic file-read]
  ]]
  file-close
   ask patches with [m-height > 0]
    [set m-traffic 2]
  ]

  if level = "no" [
    file-open "A_Floorplan_key_030320.txt"
    while [not file-at-end?] [
    foreach sort patches [ [?1] ->
  ask ?1 [
   set p-traffic file-read]
  ]]
  file-close
  ask patches with [m-height > 0]
    [set m-traffic 1]
  ]

  if level = "production" [
    file-open "traffic_production_030320.txt"
    while [not file-at-end?] [
    foreach sort patches [ [?1] ->
  ask ?1 [
   set p-traffic file-read]
  ]]
  file-close
  ask patches with [m-height > 0]
    [set m-traffic 2]
  ]
  traffic
end

to setup-weekly-schedule ; reads in production schedule
  file-close-all
  file-open "A_hourly-events.csv" ; change name of file
  set sunday csv:from-row file-read-line
  set monday csv:from-row file-read-line
  set tuesday csv:from-row file-read-line
  set wednesday csv:from-row file-read-line
  set thursday csv:from-row file-read-line
  set friday csv:from-row file-read-line
  set saturday csv:from-row file-read-line
  file-close
end

to listeria-setup
  set color yellow
  set shape "bug"
  set size 0.5
end

to listeria-initial-zone-placement
  let load (10 ^ initial-contamination-level)
  let choice one-of zones with [ z-listeria = 0 ]
  if choice != nobody [move-to choice
  ask choice [ set z-listeria-concentration (load) set z-listeria (z-listeria-concentration * [z-area] of self)]]
end

to listeria-initial-environment-placement
  let load (10 ^ initial-contamination-level)
  let choice one-of patches with [ (p-listeria = 0) and (pcolor != black)]
  move-to choice
  ask choice [ set p-listeria-concentration (load) set p-listeria (p-listeria-concentration * [p-area] of self)]
end

to set-sampling
   set start-sampling? false
   print ( word "start-sampling:" start-sampling? )
   if sampling-time != "none" [
   ;set takeSamples 0
   set samples-taken 0
   ;set sample-day one-of [2 3 4 5 6] ;changed this to M-F for scenarios
   set sample-day one-of [4] ;changed this to W to facilitate comparisons
   set sample-week time-of-simulation
   print ( word "sample-day:" sample-day )
   ]
end

;this function is executed when "go" is clicked
to go
  if not any? turtles [ stop ]
  ;prints day of the week on the interface
  if ticks mod 168 = 0 [ output-print "Sunday"]
  if ticks mod 168 = 24 [ output-print "Monday"]
  if ticks mod 168 = 48 [ output-print "Tuesday"]
  if ticks mod 168 = 72 [ output-print "Wednesday"]
  if ticks mod 168 = 96 [ output-print "Thursday"]
  if ticks mod 168 = 120 [ output-print "Friday"]
  if ticks mod 168 = 144 [ output-print "Saturday"]

  if ticks mod 24 = 0
  [set day (day + 1)
    set-production-schedule
    print (word "day:" day)
   if day = sample-day and week = sample-week [set start-sampling? true run sampling-sites]
   print (word "possible-z-sites:" possible-z-sites)
   print (word "possible-p-sites:" possible-p-sites)
  ]

  increment-contam-time
  if ticks = randomEvent [random-noise]
  set-hourly-event
  run-event
  update-listeria
  update-niche
  grow
  run environment-view

;write-output 1 ; uncomment for a single iteration if want time series contamination data
  ;run simulation but only collect data during last week - here is where you can add outputs specific to a facility
  let start-hour 24 + (24 * (time-of-simulation - 1) * 7)
  let stop-hour (24 * (time-of-simulation) * 7)
  if ticks = start-hour [ set-initial-values set collect-data? true ]
  if ticks = stop-hour
  [ set total-contam-events (zone-to-zone + zone-to-patch + condensation-drip + introduction-zone4 + introduction-food + random-event + employee-FCS + NFCS-FCS )
    set median-z1-time median [z-time-contaminated] of zones with [z-category = 1]
    set median-z2-time median [z-time-contaminated] of zones with [z-category = 2]
    set median-z3-time median [z-time-contaminated] of zones with [z-category = 3]
    set median-e-time median [z-time-contaminated] of zones with [z-item-name = "employee" ]
    set median-p-time median [p-time-contaminated] of patches with [p-water != 0]
    set median-fwj-time median [p-time-contaminated] of patches with [p-water = 35]
    if any? zones with [z-category = 1 and z-max-contam-bout != 0]
     [set median-z1-max-contam-bout median [z-max-contam-bout] of zones with [z-category = 1 and z-max-contam-bout != 0]]
    if any? zones with [z-category = 2 and z-max-contam-bout != 0]
     [set median-z2-max-contam-bout median [z-max-contam-bout] of zones with [z-category = 2 and z-max-contam-bout != 0]]
    if any? zones with [z-category = 3 and z-max-contam-bout != 0]
     [set median-z3-max-contam-bout median [z-max-contam-bout] of zones with [z-category = 3 and z-max-contam-bout != 0]]
    if any? zones with [z-item-name = "employee" and z-max-contam-bout != 0]
     [set median-e-max-contam-bout median [z-max-contam-bout] of zones with [z-item-name = "employee" and z-max-contam-bout != 0]]
    if any? patches with [p-water != 0 and p-max-contam-bout != 0]
     [set median-p-max-contam-bout median [p-max-contam-bout] of patches with [p-water != 0 and p-max-contam-bout != 0]]
    if any? patches with [p-water = 35 and p-max-contam-bout != 0]
     [ set median-fwj-max-contam-bout median [p-max-contam-bout] of patches with [p-water = 35 and p-max-contam-bout != 0]]
    set z1-contacts mean [z-contacts] of zones with [z-category = 1] set z1-transfers mean [z-transfers] of zones with [z-category = 1]
    set z2-contacts mean [z-contacts] of zones with [z-category = 2] set z2-transfers mean [z-transfers] of zones with [z-category = 2]
    set z3-contacts mean [z-contacts] of zones with [z-category = 3] set z3-transfers mean [z-transfers] of zones with [z-category = 3]
    set p-temp-niches sum [p-temp-niche] of patches with [p-water != 0]
    set z1-temp-niches sum [z-temp-niche] of zones with [z-category = 1]
    set z2-temp-niches sum [z-temp-niche] of zones with [z-category = 2]
    set z3-temp-niches sum [z-temp-niche] of zones with [z-category = 3]
    foreach sort-on [ who ] zones
      [ [?1] -> ask ?1
        [ set concentration-list lput z-listeria-concentration concentration-list
          set time-contaminated-list lput z-time-contaminated time-contaminated-list
          set max-contam-bout-list lput z-max-contam-bout max-contam-bout-list
          set contacts-list lput z-contacts contacts-list
          set transfers-list lput z-transfers transfers-list
          set temp-niche-list lput z-temp-niche temp-niche-list
          set niches-estab-list lput z-niches-established niches-estab-list
          set fruit-cont-list lput z-fruit fruit-cont-list
          set zone4-cont-list lput z-zone4 zone4-cont-list
          set random-cont-list lput z-random random-cont-list] ]
    calculate-sample-prevalence ; uncomment if testing sampling function against hist. data - modify that function as well
    stop]

  tick
end

to set-initial-values
  let z1 0 let z2 0 let z3 0 let floors 0 let es 0
  ;if (count zones with [z-listeria > 1]) != 0 [
  set z1 (count zones with [z-category = 1 and z-listeria > 1]) / (count zones with [z-category = 1])
  set z2 (count zones with [z-category = 2 and z-listeria > 1]) / (count zones with [z-category = 2])
  set z3 (count zones with [z-category = 3 and z-listeria > 1]) / (count zones with [z-category = 3])
  set floors (count patches with [pcolor != 0 and p-listeria > 1]) / (count patches with [pcolor != 0])
  set es (count zones with [z-item-name = "employee" and z-listeria > 1]) / (count zones with [z-category = 1])
  set zones-initial ((count zones with [z-listeria > 1]) / count zones ) * 100
  set z1-initial z1 * 100
  set z2-initial  z2 * 100
  set z3-initial z3 * 100
  set p-initial ( floors ) * 100
  set e-initial es * 100 ;with [z-item-name = "employee"]
end

to   increment-contam-time
  ask zones with [z-listeria = 0 and z-contam-counter != 0] [set z-contam-counter 0]
  ask zones with [z-listeria > 0]
  [ set z-time-contaminated z-time-contaminated + 1
    set z-contam-counter z-contam-counter + 1
    set z-max-contam-bout max( list z-contam-counter z-max-contam-bout) ]

  ask patches with [p-listeria = 0 and p-contam-counter != 0] [set p-contam-counter 0]
  ask patches with [p-listeria > 0]
  [ set p-time-contaminated p-time-contaminated + 1
    set p-contam-counter p-contam-counter + 1
    set p-max-contam-bout max( list p-contam-counter p-max-contam-bout) ]
end

to set-production-schedule
  if ticks mod 168 = 0 [
    ;write-output 3
    set week (week + 1)
    set production-schedule sunday
    reset-week]
  if ticks mod 168 = 24 [
    set production-schedule monday]
  if ticks mod 168 = 48 [
    set production-schedule tuesday]
  if ticks mod 168 = 72 [
    set production-schedule wednesday]
  if ticks mod 168 = 96 [
    set production-schedule thursday]
  if ticks mod 168 = 120 [
    set production-schedule friday]
  if ticks mod 168 = 144[
    set production-schedule saturday]
  reset-day
end

to reset-day
  ;set empty-time 0  set production-time 0
  ;set daily-samples 0 set z-sample-prev 0
  ;set clean-time 0
end

to reset-week
  set day 1
  set zone-to-zone 0 set zone-to-patch 0 set avg-listeria-transferred 0 set patch-to-patch 0
  set condensation-drip 0 set introduction-zone4 0  set introduction-food 0
  set random-event 0 set employee-FCS 0  set NFCS-FCS 0 set sum-fruit-load 0 set sum-fruit-transfer 0
  set-sampling ; uncomment if testing sampling
end

to calculate-sample-prevalence ; here is where you need to define global outputs based on facility, agents, and sampling data - these are collected from sample function so need to match there
  ifelse daily-samples != 0 [set daily-prev (daily-prev / daily-samples)] [set daily-prev 0]
  ifelse beg-samples != 0 [set z-sample-beg-prev  (z-sample-beg-prev / beg-samples )] [set z-sample-beg-prev 0]
  ifelse mid-samples != 0 [set z-sample-mid-prev  (z-sample-mid-prev / mid-samples )] [set z-sample-mid-prev 0]
  ifelse post-samples != 0 [  set z-sample-post-prev  (z-sample-post-prev / post-samples )] [set z-sample-post-prev 0]
  ifelse z1-samples != 0 [  set z1-prev (z1-prev / z1-samples)] [set z1-prev 0]
  ifelse z2-samples != 0 [  set z2-prev (z2-prev / z2-samples)] [set z2-prev 0]
  ifelse z3-samples != 0 [  set z3-prev (z3-prev / z3-samples)] [set z3-prev 0]
  ifelse belt-samples != 0 [  set belt-prev (belt-prev / belt-samples )] [set belt-prev 0]
  ifelse cart-samples != 0 [  set cart-prev (cart-prev / cart-samples )] [set cart-prev 0]
  ifelse control-panel-samples != 0 [  set control-panel-prev (control-panel-prev / control-panel-samples )] [set control-panel-prev 0]
  ifelse door-samples != 0 [  set door-prev (door-prev / door-samples )] [set door-prev 0]
  ifelse drain-samples != 0 [  set drain-prev (drain-prev / drain-samples )] [set drain-prev 0]
  ifelse employee-samples != 0 [  set employee-prev (employee-prev / employee-samples )] [set employee-prev 0]
  ifelse flume-samples != 0 [  set flume-prev (flume-prev / flume-samples )] [set flume-prev 0]
  ifelse frame-samples != 0 [  set frame-prev (frame-prev / frame-samples )] [set frame-prev 0]
  ifelse frame-in-samples != 0 [  set frame-in-prev (frame-in-prev / frame-in-samples )] [set frame-in-prev 0]
  ifelse ladder-samples != 0 [  set ladder-prev (ladder-prev / ladder-samples )] [set ladder-prev 0]
  ifelse misc-samples != 0 [  set misc-prev (misc-prev / misc-samples )] [set misc-prev 0]
  ifelse packing-samples != 0 [  set packing-prev (packing-prev / packing-samples )] [set packing-prev 0]
  ifelse squeegee-samples != 0 [  set squeegee-prev (squeegee-prev / squeegee-samples )] [set squeegee-prev 0]
  ifelse trash-gray-samples != 0 [  set trash-gray-prev (trash-gray-prev / trash-gray-samples )] [set trash-gray-prev 0]
  ifelse trash-white-samples != 0 [  set trash-white-prev (trash-white-prev / trash-white-samples )] [set trash-white-prev 0]
  ifelse trash-yellow-samples != 0 [  set trash-yellow-prev (trash-yellow-prev / trash-yellow-samples )] [set trash-yellow-prev 0]
  ifelse wall-samples != 0 [  set wall-prev (wall-prev / wall-samples )] [set wall-prev 0]
  ifelse weigher-samples != 0 [  set weigher-prev (weigher-prev / weigher-samples )] [set weigher-prev 0]
  ifelse a-b-samples != 0 [  set a-b-prev (a-b-prev / a-b-samples )] [set a-b-prev 0]
  ifelse a-flume-samples != 0 [  set a-flume-prev (a-flume-prev / a-flume-samples )] [set a-flume-prev 0]
  ifelse a-line-samples != 0 [  set a-line-prev (a-line-prev / a-line-samples )] [set a-line-prev 0]
  ifelse a-under-room-samples != 0 [  set a-under-room-prev (a-under-room-prev / a-under-room-samples )] [set a-under-room-prev 0]
  ifelse b-flume-samples != 0 [  set b-flume-prev (b-flume-prev / b-flume-samples )] [set b-flume-prev 0]
  ifelse b-line-samples != 0 [  set b-line-prev (b-line-prev / b-line-samples )] [set b-line-prev 0]
  ifelse coring-samples != 0 [  set coring-prev (coring-prev / coring-samples )] [set coring-prev 0]
  ifelse elevator-samples != 0 [  set elevator-prev (elevator-prev / elevator-samples )] [set elevator-prev 0]
  ifelse hopper-samples != 0 [  set hopper-prev (hopper-prev / hopper-samples )] [set hopper-prev 0]
  ifelse packing--samples != 0 [  set packing--prev (packing--prev / packing--samples )] [set packing--prev 0]
  ifelse scoring-samples != 0 [  set scoring-prev (scoring-prev / scoring-samples )] [set scoring-prev 0]
  ifelse sorting-samples != 0 [  set sorting-prev (sorting-prev / sorting-samples )] [set sorting-prev 0]
  ifelse spin-dryer-samples != 0 [  set spin-dryer-prev (spin-dryer-prev / spin-dryer-samples )] [set spin-dryer-prev 0]
  ifelse weigher-platform-samples != 0 [  set weigher-platform-prev (weigher-platform-prev / weigher-platform-samples )] [set weigher-platform-prev 0]
  ifelse A-samples != 0 [  set A-prev (A-prev / A-samples )] [set A-prev 0]
  ifelse B-samples != 0 [  set B-prev (B-prev / B-samples )] [set B-prev 0]
  ifelse floor-samples != 0 [  set floors-prev (floors-prev / floor-samples )] [set floors-prev 0]
end

to set-hourly-event
  let hour ticks mod 24
  set event item hour production-schedule
end

to run-event ; here is where you need to bring everything together based on production schedule these events need to be terms used in the production schedule .csv file

; are there times when the room is empty?
  if event = "empty" [
    if not all? employees [hidden?] [clear-employees]
    update-traffic ("no")
    dry-water
    ]

 ; these things occur during pre-op inspection
  if event = "pre-op" [
    ;write-output 2
    ask n-of 2 employees [ ; how many employees are in the room for pre-op? - change 5 to whatever is appropriate
      set hidden? false
      ask my-proximity-links [set hidden? false]]
    update-water ("low") ; what is the water level at preop change "low" to some other level that you have defined in update-water function
    update-traffic ("low") ; what is the traffic level at preop change "low" as with above for water
    dry-water
    zone4-introduction
    patch-spread-traffic
    patch-spread-water
    ;if start-sampling? [sample]
    if collect-data? [collect-data]
    clear-employees
    set clean-time 0
    set production-time 0
    print (word "pre-op")
    print (word "production-time:" production-time)
    ]

  ; these things occur during production
  if event = "production" [
    if all? employees [hidden?] [
      ask n-of (random-in-range 26 34) employees [ ; how many employees are in the room during production?
        set hidden? false
        ask my-proximity-links [set hidden? false]]]
    update-water ("high") ; update water level to high during all hours of production
    ifelse production-time = 0 or production-time = 8 or production-time = 17 ; update traffic level at shift change?
       ;[update-traffic ("high")]
       [update-traffic ("high") zone4-introduction]
       [update-traffic ("production")]
    zone-spread
    patch-spread-traffic
    patch-spread-water
    zone-to-patchBelow
    condensation
    food-introduction
    zone4-introduction
    ;dry-water
    if production-time = 0 and sampling-time = "first-hour" and start-sampling? [
      print ( "to-sample-first-hour" )
      sample
      print ( "sampled-first-hour" )
    ]

    if production-time = 3 and sampling-time = "mid-production" and start-sampling? [sample]
    if production-time = 9 and sampling-time = "second-shift" and start-sampling? [sample]
    ;if start-sampling? [sample]
    if collect-data? [collect-data]
    if production-time = 8 or production-time = 17 [clear-employees] ;clear room at the end of each shift
    set production-time (production-time + 1)
    print (word "production-time:" production-time)
    ]

  ; these things occur during cleaning and sanitation
  if event = "clean" [

    ifelse day = 7
    [update-traffic ("low")
     ask zones [set z-water 3]
     if all? employees [hidden?] [
      ask n-of 15 employees [ ; how many employees are in the room for rinse?
        set hidden? false
        ask my-proximity-links [set hidden? false]]]
    zone4-introduction
    condensation
    patch-spread-water
    patch-spread-traffic]

    [if clean-time = 0 [ ;update water and traffic levels as appropriate
      update-traffic ("low")
      ;update-water ("low")
      ask zones [set z-water 3]]
    if all? employees [hidden?] [
      ask n-of 6 employees [ ; how many employees are in the room for cleaning?
        set hidden? false
        ask my-proximity-links [set hidden? false]]]
  ;  zone4-introduction
    condensation
    patch-spread-water
    patch-spread-traffic
    ; how long is the cleaning and san shift?  we assume sanitation log reduction only occurs in last hour when sanitizer is applied
    if day = 2 and clean-time = 10 [ ;Monday
        clean
        ask patches with [p-water >= 1 and p-water < 5] [set p-water 4]
        ask zones with [z-cleanable?] [set z-water 4]
        clear-employees]
    if day > 2 and day < 7 and clean-time = 4 [ ;Tues-Friday
        clean
        ask patches with [p-water >= 1 and p-water < 5] [set p-water 4]
        ask zones with [z-cleanable?] [set z-water 4]
        clear-employees]
    set clean-time (clean-time + 1)]
    ]
end

to clear-employees
  ask employees [
      set hidden? true
      ask my-proximity-links [set hidden? true]
      ask listeria-here [ die ] ; assume that listeria are removed from employees when they are not in the room
      set z-listeria 0
      set z-listeria-concentration 0]
end

to collect-data ; this is where you need to collect data for global outputs defined at the beginning - these are specific to a facility and agents

  if event = "pre-op" [
  set zones-pre-prev ((count zones with [z-listeria > 0]) / count zones ) * 100
  set z1-pre ((count zones with [z-category = 1 and z-listeria > 0 ]) / count zones with [z-category = 1]) * 100
  set z2-pre ((count zones with [z-category = 2 and z-listeria > 0 ] )/ count zones with [z-category = 2] ) * 100
  set z3-pre ((count zones with [z-category = 3 and z-listeria > 0]) / count zones with [z-category = 3]) * 100
  set p-pre ((count patches with [pcolor != 0 and p-listeria > 0 ]) / count patches with [pcolor != 0])  * 100
  set e-pre ((count zones with [z-item-name = "employee" and z-listeria > 0 ]) / count zones with [z-category = 1]) * 100
  ]
  if (event = "production" and production-time = 0) [
  set zones-beg-prev ((count zones with [z-listeria > 0])/ count zones ) * 100
  set z1-beg ((count zones with [z-category = 1 and z-listeria > 0 ]) / count zones with [z-category = 1]) * 100
  set z2-beg ((count zones with [z-category = 2 and z-listeria > 0 ] )/ count zones with [z-category = 2] ) * 100
  set z3-beg ((count zones with [z-category = 3 and z-listeria > 0]) / count zones with [z-category = 3]) * 100
  set p-beg ((count patches with [pcolor != 0 and p-listeria > 0 ]) / count patches with [pcolor != 0])  * 100
  set e-beg ((count zones with [z-item-name = "employee" and z-listeria > 0 ]) / count zones with [z-category = 1]) * 100
  ]
  if (event = "production" and production-time = 8) [
  set zones-mid-prev ((count zones with [z-listeria > 0])/ count zones ) * 100
  set z1-mid ((count zones with [z-category = 1 and z-listeria > 0 ]) / count zones with [z-category = 1]) * 100
  set z2-mid ((count zones with [z-category = 2 and z-listeria > 0 ] )/ count zones with [z-category = 2] ) * 100
  set z3-mid ((count zones with [z-category = 3 and z-listeria > 0]) / count zones with [z-category = 3]) * 100
  set p-mid ((count patches with [pcolor != 0 and p-listeria > 0 ]) / count patches with [pcolor != 0])  * 100
  set e-mid ((count zones with [z-item-name = "employee" and z-listeria > 0 ]) / count zones with [z-category = 1]) * 100
  ]
  if (event = "production" and production-time = 17) [
  set zones-post-prev ((count zones with [z-listeria > 0])/ count zones ) * 100
  set z1-post ((count zones with [z-category = 1 and z-listeria > 0 ]) / count zones with [z-category = 1]) * 100
  set z2-post ((count zones with [z-category = 2 and z-listeria > 0 ] )/ count zones with [z-category = 2] ) * 100
  set z3-post ((count zones with [z-category = 3 and z-listeria > 0]) / count zones with [z-category = 3]) * 100
  set p-post ((count patches with [pcolor != 0 and p-listeria > 0 ]) / count patches with [pcolor != 0])  * 100
  set e-post ((count zones with [z-item-name = "employee" and z-listeria > 0 ]) / count zones with [z-category = 1]) * 100
  ]

  if day = 2 [
    set m-z-sample-prev z-sample-prev  set m-zones-pre zones-pre-prev set m-zones-beg zones-beg-prev set m-zones-post zones-post-prev set m-zones-mid zones-mid-prev
    set m-z1-pre z1-pre set m-z2-pre z2-pre set m-p-pre p-pre set m-e-pre e-pre set m-z3-pre z3-pre
    set m-z1-beg z1-beg set m-z2-beg z2-beg set m-p-beg p-beg set m-e-beg e-beg set m-z3-beg z3-beg
    set m-z1-mid z1-mid set m-z2-mid z2-mid set m-p-mid p-mid set m-e-mid e-mid set m-z3-mid z3-mid
    set m-z1-post z1-post set m-z2-post z2-post set m-p-post p-post set m-e-post e-post set m-z3-post z3-post]
  if day = 3 [
    set t-z-sample-prev z-sample-prev set t-zones-pre zones-pre-prev set t-zones-beg zones-beg-prev set t-zones-post zones-post-prev set t-zones-mid zones-mid-prev
    set t-z1-pre z1-pre set t-z2-pre z2-pre set t-p-pre p-pre set t-e-pre e-pre set t-z3-pre z3-pre
    set t-z1-beg z1-beg set t-z2-beg z2-beg set t-p-beg p-beg set t-e-beg e-beg set t-z3-beg z3-beg
    set t-z1-mid z1-mid set t-z2-mid z2-mid set t-p-mid p-mid set t-e-mid e-mid set t-z3-mid z3-mid
    set t-z1-post z1-post set t-z2-post z2-post set t-p-post p-post set t-e-post e-post set t-z3-post z3-post]
  if day = 4 [
    set w-z-sample-prev z-sample-prev set w-zones-pre zones-pre-prev set w-zones-beg zones-beg-prev set w-zones-post zones-post-prev set w-zones-mid zones-mid-prev
    set w-z1-pre z1-pre set w-z2-pre z2-pre set w-p-pre p-pre set w-e-pre e-pre set w-z3-pre z3-pre
    set w-z1-beg z1-beg set w-z2-beg z2-beg set w-p-beg p-beg set w-e-beg e-beg set w-z3-beg z3-beg
    set w-z1-mid z1-mid set w-z2-mid z2-mid set w-p-mid p-mid set w-e-mid e-mid set w-z3-mid z3-mid
    set w-z1-post z1-post set w-z2-post z2-post set w-p-post p-post set w-e-post e-post set w-z3-post z3-post]
  if day = 5 [
    set r-z-sample-prev z-sample-prev set r-zones-pre zones-pre-prev set r-zones-beg zones-beg-prev set r-zones-post zones-post-prev set r-zones-mid zones-mid-prev
    set r-z1-pre z1-pre set r-z2-pre z2-pre set r-p-pre p-pre set r-e-pre e-pre set r-z3-pre z3-pre
    set r-z1-beg z1-beg set r-z2-beg z2-beg set r-p-beg p-beg set r-e-beg e-beg set r-z3-beg z3-beg
    set r-z1-mid z1-mid set r-z2-mid z2-mid set r-p-mid p-mid set r-e-mid e-mid set r-z3-mid z3-mid
    set r-z1-post z1-post set r-z2-post z2-post set r-p-post p-post set r-e-post e-post set r-z3-post z3-post]
  if day = 6 [
    set f-z-sample-prev z-sample-prev set f-zones-pre zones-pre-prev set f-zones-beg zones-beg-prev set f-zones-post zones-post-prev set f-zones-mid zones-mid-prev
    set f-z1-pre z1-pre set f-z2-pre z2-pre set f-p-pre p-pre set f-e-pre e-pre set f-z3-pre z3-pre
    set f-z1-beg z1-beg set f-z2-beg z2-beg set f-p-beg p-beg set f-e-beg e-beg set f-z3-beg z3-beg
    set f-z1-mid z1-mid set f-z2-mid z2-mid set f-p-mid p-mid set f-e-mid e-mid set f-z3-mid z3-mid
    set f-z1-post z1-post set f-z2-post z2-post set f-p-post p-post set f-e-post e-post set f-z3-post z3-post]
end

to write-output [file-number]

  if file-number = 1 [
  ;OUTPUT FILE #1 DATA
  let zone-data [  ] ; defines an empty list each tick (as a local var.) that holds concentration for each zone
  set zone-data fput ticks zone-data
  foreach sort-on [ who ] zones
   [ [?1] -> ask ?1
     [ set zone-data lput z-listeria-concentration zone-data ] ]
  file-open  "TimeSeriesZoneConcen.csv"
  let data-to-csv csv:to-row zone-data
  file-print data-to-csv
  file-close]
end

; Verhulst discrete growth function for nodes and patches - don't need to change
; only growth when water is present (2 or 3) and same growth rate
to grow
  ask zones with [  (z-water > 1 and z-water <= 3) ]
    [let N (z-listeria-concentration)
     set z-listeria-concentration ( ( K * N * e ^ mu-max)/ (K + N * ((e ^ mu-max) - 1 )))
     set z-listeria round (z-listeria-concentration * [z-area] of self)]

  ask patches with [ (p-water > 1 and p-water <= 3)]
    [let Np (p-listeria-concentration)
     set p-listeria-concentration ( ( K * Np * e ^ mu-max)/ (K + Np * ((e ^ mu-max) - 1 )))
     set p-listeria round (p-listeria-concentration * [p-area] of self)]

  ask patches with [ (c-water > 1 and c-water <= 3)]
    [let Nc (c-listeria-concentration)
     set c-listeria-concentration ( ( K * Nc * e ^ mu-max)/ (K + Nc * ((e ^ mu-max) - 1 )))
     set c-listeria round (c-listeria-concentration * [p-area] of self)]

  ask patches with [ (m-water > 1 and m-water <= 3)]
    [let Nm (m-floor-listeria-concentration)
     set m-floor-listeria-concentration ( ( K * Nm * e ^ mu-max)/ (K + Nm * ((e ^ mu-max) - 1 )))
     set m-floor-listeria round (m-floor-listeria-concentration * [p-area] of self)]
end

to update-niche
  ask zones with [not z-cleanable?]
  [if z-listeria = 0
    [ if z-niches-established = 0 [set z-time-to-niche z-time-to-niche + 1
      set z-niche-length 0]]
  if z-listeria > 0
    [ set z-niche-length z-niche-length + 1]]
end

to dry-water
  ask zones with [z-water > 1 and z-water <= 3]
    [set z-water (z-water - 0.2)]
  ask zones with [(z-water = 4) and z-cleanable?] [set z-water 1]
  ask zones with [(z-water = 4) and not z-cleanable?] [set z-water 3]

  ask patches with [(p-water > 1) and (p-water <= 3)]
    [set p-water (p-water - 0.2)]
  ask patches with [p-water = 4] [set p-water 1 ]
  ask patches with [c-water > 1 and c-water <= 3]
    [set c-water (c-water - 0.2)]
  ask patches with [m-water > 1 and m-height > 0]
    [set m-water (m-water - 0.2)]
end

to update-listeria
; ask listeria-on patches with [(p-listeria < 1) and (c-listeria < 1) and (m-floor-listeria < 0) and (m-ceiling-listeria < 0)]
;   [if (not any? zones-here)
;      [ask self [ die ]]
;   ]

ask patches with [(p-listeria = 0) and (c-listeria = 0) and (m-floor-listeria = 0) and (m-ceiling-listeria = 0)]
  [if (not any? zones-here)
    [ask listeria-here [die]]
  ]

 ask zones with [ (z-listeria < 1) ]
  [ let z count listeria-here
    if z >= 1
    [ask patch-here
    [ if (p-listeria = 0) and (c-listeria = 0) and (m-floor-listeria = 0) and (m-ceiling-listeria = 0)
      [ask listeria-here [ die ] ]
   ] ]]

 ask zones with [ (z-listeria >= 1) ]
  [ let z count listeria-here
    if z > 1
    [ask listeria-here [ die ]
     hatch-listeria 1 [ listeria-setup setxy xcor ycor]]
   ]
end

to zone-spread ;TAKEN FROM FDA DELI MODEL - may need modification if using zone 4 agents - matrices will be 5x5 instead of 4x4 as here
  ask zones with [z-listeria > 0]
  [ let N1 z-listeria
    let me [who] of self
    let num ([z-category] of self) - 1
    let NO (count zones with [z-category = (num + 1)])
    if z-item-name = "employee" [set num 3]
    let neighbor-zones (turtle-set (link-neighbors) (out-link-neighbors) )

    if (neighbor-zones != nobody)
    [ask neighbor-zones

       [let num2 ([z-category] of self) - 1
        let NO2 (count zones with [z-category = (num2 + 1)])
        if z-item-name = "employee" [set num2 3]

        let prob-12  (item num2 (item num p-transfer))
        let prob-21  (item num (item num2 p-transfer))
        if not out-link-neighbor? a-zone me [set prob-21 0]

        let trans-12  (item num2 (item num tc)) ;we keep track of this during variability simulations
        let trans-21  (item num (item num2 tc))
        let N2 z-listeria
        let z count (listeria-here)
        let T12 0
        let T21 0

        ifelse random-float 100 < prob-12
         [ set T12 min (list 1 (10 ^ trans-12)) ;(10 ^ random-normal (tc12) (std12)) ;
           ask myself [set z-transfers z-transfers + 1]
           ask self [set z-contacts z-contacts + 1]
           ifelse (num = 3 and num2 = 0) [set employee-FCS (employee-FCS + 1)]
             [ifelse (num = 1  and num2 = 0) [set NFCS-FCS (NFCS-FCS + 1)] [set zone-to-zone (zone-to-zone + 1)]]]
         [ set T12 0 ]
        ifelse random-float 100 < prob-21
         [ set T21 min (list 1 (10 ^ trans-21)) ;(10 ^ random-normal (tc21) (std21)) ;
           ask self [set z-transfers z-transfers + 1]
           ask myself [set z-contacts z-contacts + 1]
           ifelse (num2 = 3 and num = 0) [set employee-FCS (employee-FCS + 1)]
             [ifelse (num2 = 1 and num = 0) [set NFCS-FCS (NFCS-FCS + 1)] [set zone-to-zone (zone-to-zone + 1)]]]
         [ set T21 0 ]
        ;print (word "T12:" T12)
        let x11 (random-binomial (N1) (1 - T12))
        let x21 (random-binomial (N2) (T21))
        if zone-to-zone != 0 [set avg-listeria-transferred (avg-listeria-transferred + (N1 - x11 + x21) / (zone-to-zone + employee-FCS + NFCS-FCS))]
        set z-listeria (N1 + N2 - (x11 + x21))
        set z-listeria-concentration (z-listeria / [z-area] of self)

        set N1 (x11 + x21)
        ;print (word "N1:" N1)

        if (([z-listeria] of self) > 0) and (z = 0)
         [ hatch-listeria 1 [ listeria-setup setxy xcor ycor]
           if ([not z-cleanable?] of self) [set z-niches-established z-niches-established + 1]]

       ] ]
    set z-listeria N1
    set z-listeria-concentration (z-listeria / [z-area] of self)

  ]
end

to patch-spread-traffic
  let p-i (1 / count patches with [(p-traffic != 0) and (p-traffic != 5)])

  ask patches with [p-listeria-concentration > 1 and p-traffic > 0]
  [ let x p-traffic
    let z 0
    let M count neighbors with [p-water != 0 and p-water != 35]
    let possible-patches ( patch-set neighbors with [p-traffic >= x])
    let num-patches count possible-patches
    ifelse x = 3 [set z (p-patch-spread * p-i * traffic-high * (num-patches / M))]
    [ifelse x = 2 [set z (p-patch-spread * p-i * traffic-low * (num-patches / M))] [set z (p-patch-spread * p-i * traffic-neg * (num-patches / M))]]
    let N p-listeria ;assume traffic hits entire patch
    let trans (spread-tc / 100)
    ;print (word "z:" z)
    let p 1 - exp (- z)
    set p p * 100
    ;print (word "p:" p)
    if (possible-patches != nobody) and (random-float 100 <  p )
    [ set p-listeria (p-listeria - (N * trans ))
      set p-listeria-concentration (p-listeria / [p-area] of self)
      ask one-of possible-patches
       [ let y count (listeria-here)
         set patch-to-patch (patch-to-patch + 1)
         set p-listeria (p-listeria + (N * trans ))
         set p-listeria-concentration (p-listeria / [p-area] of self)
         if (([p-listeria] of self) > 0) and (y = 0)
           [ sprout-listeria 1 [ listeria-setup setxy xcor ycor]]
       ]] ]

  ;;include this if including the mezzanine in model
    ask patches with [m-floor-listeria > 0 and m-traffic > 0]
  [ let x m-traffic
    let z 0
    let M count neighbors with [m-height != 0]
    let possible-patches ( patch-set neighbors with [m-traffic >= x])
    let num-patches count possible-patches
    ifelse x = 3 [set z (p-patch-spread * p-i * traffic-high * (num-patches / M))]
    [ifelse x = 2 [set z (p-patch-spread * p-i * traffic-low * (num-patches / M))] [set z (p-patch-spread * p-i * traffic-neg * (num-patches / M))]]
    let N m-floor-listeria ;assume traffic hits entire patch
    let trans (spread-tc / 100)
    ;print (word "z:" z)
    let p 1 - exp (- z)
    set p p * 100
    ;print (word "p:" p)
    if (possible-patches != nobody) and (random-float 100 <  p )
    [ set m-floor-listeria (m-floor-listeria - (N * trans ))
      set m-floor-listeria-concentration (m-floor-listeria / [p-area] of self)
      ask one-of possible-patches
       [ let y count (listeria-here)
         set patch-to-patch (patch-to-patch + 1)
         set m-floor-listeria (m-floor-listeria + (N * trans ))
         set m-floor-listeria-concentration (m-floor-listeria / [p-area] of self)
         if (([m-floor-listeria] of self) > 0) and (y = 0)
           [ sprout-listeria 1 [ listeria-setup setxy xcor ycor]]
       ]]  ]

end

to patch-spread-water
 ask patches with [p-listeria-concentration > 0 and p-water = 3]
  [ let N round (p-listeria-concentration) ;assume 1 ml moves
    let z p-water-spread
    let trans (spread-tc / 100)
    let possible-patches ( patch-set neighbors with [p-water = 3 or p-item-name = "floor-wall-juncture"])
    let num-patches count possible-patches
    set z z * 100
    if (num-patches != 0) and ((random-float 100 < z) or (N > (1E8)))
    [let listeria-transfer-per-patch ( (trans * N ) / num-patches )
    set p-listeria (p-listeria - (trans * N))
    set p-listeria-concentration (p-listeria / [p-area] of self)
    ask possible-patches
       [let y count (listeria-here)
        set patch-to-patch (patch-to-patch + 1)
        set p-listeria (p-listeria + listeria-transfer-per-patch)
        set p-listeria-concentration (p-listeria / [p-area] of self)
        if (([p-listeria] of self) > 0) and (y = 0)
        [sprout-listeria 1 [ listeria-setup setxy pxcor pycor]]
      ]] ]

 ;;include this if including the mezzanine in model
 ask patches with [m-floor-listeria > 0 and m-water = 3]
  [ let N round (m-floor-listeria-concentration) ;assume 1 ml moves
    let z p-water-spread
    let trans (spread-tc / 100)
    let possible-patches ( patch-set neighbors with [m-water = 3])
    let num-patches count possible-patches
    set z z * 100
    if (num-patches != 0) and ((random-float 100 < z) or (N > K))
    [let listeria-transfer-per-patch ( (trans * N ) / num-patches )
    set m-floor-listeria (m-floor-listeria - (trans * N))
    set m-floor-listeria-concentration (m-floor-listeria / [p-area] of self)
    ask possible-patches
       [let y count (listeria-here)
        set patch-to-patch (patch-to-patch + 1)
        set m-floor-listeria (m-floor-listeria + listeria-transfer-per-patch)
        set m-floor-listeria-concentration (m-floor-listeria / [p-area] of self)
        if (([m-floor-listeria] of self) > 0) and (y = 0)
        [sprout-listeria 1 [ listeria-setup setxy pxcor pycor]]
      ]]   ]

end


to zone-to-patchBelow
   ask zones with [z-listeria > 1 and z-category = 1]
 [ let N (z-listeria-concentration * 5)
   let trans (spread-tc / 100)
   ;;include this section if modeling the mezzanine
   ifelse z-height > 12 ;;change this number depending on height of mezzanine
   [let possible-patch one-of neighbors with [m-height != 0]
   if (possible-patch != nobody) and (random-float 100 < p-food-dropped)
    [ set z-listeria (z-listeria - (N * trans))
      set z-listeria-concentration  (z-listeria / [z-area] of self)
      set zone-to-patch (zone-to-patch + 1)
      ask possible-patch
      [ifelse (p-listeria = 0 and c-listeria = 0 and m-floor-listeria = 0 and m-ceiling-listeria = 0)
         [set m-floor-listeria (N * trans)
          set m-floor-listeria-concentration (m-floor-listeria / [p-area] of self)
          sprout-listeria 1 [ listeria-setup setxy xcor ycor]]
         [set m-floor-listeria (m-floor-listeria + (trans * N ))
          set m-floor-listeria-concentration (m-floor-listeria / [p-area] of self)]]
    ]]
   [;; until here
     let possible-patch one-of neighbors with [p-water != 0]
   if (possible-patch != nobody) and (random-float 100 < p-food-dropped)
    [ set z-listeria (z-listeria - (N * trans))
      set z-listeria-concentration  (z-listeria / [z-area] of self)
      set zone-to-patch (zone-to-patch + 1)
      ask possible-patch
      [ifelse (p-listeria = 0 and c-listeria = 0 and m-floor-listeria = 0 and m-ceiling-listeria = 0)
         [set p-listeria (N * trans)
          set p-listeria-concentration (p-listeria / [p-area] of self)
          sprout-listeria 1 [ listeria-setup setxy pxcor pycor]]
         [set p-listeria (p-listeria + (trans * N ))
          set p-listeria-concentration (p-listeria / [p-area] of self)]]
    ]
   ] ;; and remove this bracket too if not using mezzanine

 ]
end

to condensation
ask patches with [(c-water = 3 or m-water = 3) and (pxcor > 0 and pxcor < 93 and pycor > 0 and pycor < 40)] ;;change coordinates to match the ceiling area of the production room
[let prob p-condensation
   let possible-ceiling count patches with [p-water != 0 and p-water != 35 and p-water != 5 and p-water != 6 and pxcor > 0 and pxcor < 93 and pycor > 0 and pycor < 40]
   let possible-mezzanine count patches with [m-height != 0]
   ifelse [m-water = 3] of self
     [ if random-float 100 < (prob * 1 / possible-mezzanine)
     [set condensation-drip (condensation-drip + 1)
     let N m-ceiling-listeria-concentration * 2 ; 2ml droplet
     set m-ceiling-listeria (m-ceiling-listeria - N)
     set m-ceiling-listeria-concentration (m-ceiling-listeria / [p-area] of self)
     ifelse (any? zones-here) ; if there is a zone below, condensation falls there
      [ask zones-here with-max [z-height]
       [ifelse (z-listeria = 0)
         [ set z-listeria (N)
           hatch-listeria 1 [listeria-setup setxy xcor ycor]
           set z-water 2
           if not z-cleanable? [set z-niches-established z-niches-established + 1]]
         [set z-listeria (z-listeria + N)
           set z-listeria-concentration (z-listeria / [z-area] of self)
           set z-water 2]]]
      [set p-listeria (p-listeria + N)
       set p-listeria-concentration (p-listeria / [p-area] of self)
       if p-water < 2 [set p-water 2 ] ; if no zone below, condensation falls to floor
      ]]]

   [if random-float 100 < (prob * 1 / possible-ceiling)   ;probability that condensation falls from specific patch
    [set condensation-drip (condensation-drip + 1)
     let N c-listeria-concentration * 2 ; 2ml droplet
     set c-listeria (c-listeria - N)
     set c-listeria-concentration (c-listeria / [p-area] of self)
     ifelse (any? zones-here) ; if there is a zone below, condensation falls there
      [ask zones-here with-max [z-height]
       [ifelse (z-listeria = 0)
         [ set z-listeria (N)
           hatch-listeria 1 [listeria-setup setxy xcor ycor]
           set z-water 2
           if not z-cleanable? [set z-niches-established z-niches-established + 1]]
         [set z-listeria (z-listeria + N)
           set z-listeria-concentration (z-listeria / [z-area] of self)
           set z-water 2]]]
      [set p-listeria (p-listeria + N)
       set p-listeria-concentration (p-listeria / [p-area] of self)
       if p-water < 2 [set p-water 2 ] ; if no zone below, condensation falls to floor
      ]]]
  ]
end


to food-introduction ; just make sure the global variable names are changed here as well
  let conc fruit-conc
  let prevalence 0
  if day = 2 [set prevalence m-crate-prevalence]
  if day = 3 [set prevalence t-crate-prevalence]
  if day = 4 [set prevalence w-crate-prevalence]
  if day = 5 [set prevalence r-crate-prevalence]
  if day = 6 [set prevalence f-crate-prevalence]

  let bulk-crates round (random-binomial flow-rate prevalence )
  let surfaces (zones with [ z-location = "scoring" or z-location = "coring"]) ; change this to define the agents where raw material first enters production room

  repeat bulk-crates [
  ask one-of surfaces
     [let x count listeria-here
       let load (conc * crate-size) ;CFU subject to transfer
       let y (random-binomial load fruit-tc)
       set introduction-food (introduction-food + 1)
       set sum-fruit-load (sum-fruit-load + load)
       set sum-fruit-transfer (sum-fruit-transfer + y)
       ask self [set z-listeria (z-listeria + y) set z-listeria-concentration (z-listeria / [z-area] of self)   set z-fruit z-fruit + 1]
       if (([z-listeria] of self) > 0) and (x = 0)
        [hatch-listeria 1 [ listeria-setup setxy xcor ycor]
        if ([not z-cleanable?] of self) [set z-niches-established z-niches-established + 1]]]
  ]
end


to zone4-introduction
  let p  p-zone4-intro

  ask patches with [p-water = 5] [
    if random-float 100 < p [
    let sites [zones in-radius 5] of self
    let spots [(patches with [p-water != 0 and p-item-name != "floor-wall-juncture"]) in-radius 5] of self
    let n count sites
    let q count spots

    ifelse random 100 < 50 ;assume equal chance of introduction going to floor or equipment

    [if (n != 0) [
    ask one-of sites [
         let x count listeria-here
         set introduction-zone4 (introduction-zone4 + 1)
           ask self [set z-listeria (z-listeria + zone4-load) set z-listeria-concentration (z-listeria / [z-area] of self) set z-zone4 z-zone4 + 1]
           if (([z-listeria] of self) > 0) and (x = 0)
             [hatch-listeria 1 [ listeria-setup setxy xcor ycor]
             if ([not z-cleanable?] of self) [set z-niches-established z-niches-established + 1]]
    ]]]

    [if (q != 0) [
    ask one-of spots [
         let y count listeria-here
         set introduction-zone4 (introduction-zone4 + 1)
           ask self [set p-listeria (p-listeria + zone4-load) set p-listeria-concentration (p-listeria / [p-area] of self) set p-zone4 p-zone4 + 1]
           if (([p-listeria] of self) > 0) and (y = 0)
              [sprout-listeria 1 [ listeria-setup setxy pxcor pycor]]
    ]]]
    ]]
end

to random-noise
  set random-event (random-event + 1)
  if day = 2 [set m-random-event (m-random-event + 1)]
  if day = 3 [set t-random-event (t-random-event + 1)]
  if day = 4 [set w-random-event (w-random-event + 1)]
  if day = 5 [set r-random-event (r-random-event + 1)]
  if day = 6 [set f-random-event (f-random-event + 1)]


;Sets up probability of random introduction to floor (75%), ceiling (5%), and agents (20%)
ifelse random 100 < 75
 [set random-floor random-floor + 1
   ask one-of patches with [ p-water != 0 ]
     [let z count listeria-here
       ask self [set p-listeria (p-listeria + random-load) set p-listeria-concentration (p-listeria / [p-area] of self) set p-random p-random + 1]
       if (([p-listeria] of self) > 0) and (z = 0)
         [sprout-listeria 1 [ listeria-setup setxy pxcor pycor]]]]
[ifelse random 100 > 94
  [set random-ceiling random-ceiling + 1
    ask one-of patches with [ p-water != 0]
     [let z count listeria-here
       ask self [set c-listeria (c-listeria + random-load) set c-listeria-concentration (c-listeria / [p-area] of self) set p-random p-random + 1]
       if (([c-listeria] of self) > 0) and (z = 0)
         [sprout-listeria 1 [ listeria-setup setxy pxcor pycor]]]]
     [ ask one-of zones
     [set random-agent random-agent + 1
       let x count listeria-here
       ask self [set z-listeria (z-listeria + random-load) set z-listeria-concentration (z-listeria / [z-area] of self) set z-random z-random + 1]
       if (([z-listeria] of self) > 0) and (x = 0)
        [hatch-listeria 1 [ listeria-setup setxy xcor ycor]
          if ([not z-cleanable?] of self) [set z-niches-established z-niches-established + 1]]]]]


  if p-random-noise != 0 [set randomEvent ticks + (item 0 random-event-list)
    set random-event-list but-first random-event-list
    if random-event-list = [] [setup-randomEvents]]

end

to clean
  print("clean")
 let reduction 0
  if day = 2 [set reduction (m-san-reduction)]
  if day = 3 [set reduction (t-san-reduction)]
  if day = 4 [set reduction (w-san-reduction)]
  if day = 5 [set reduction (r-san-reduction)]
  if day = 6 [set reduction (f-san-reduction)]

 ask zones with [z-listeria > 0 and z-cleanable?]
    [ ifelse random 100 < prob-cleanable
      [set z-listeria (random-binomial (z-listeria) (reduction))
      set z-listeria-concentration (z-listeria / [z-area] of self)]
      [set temporary-niches temporary-niches + 1
       set z-temp-niche z-temp-niche + 1]]
 ask patches with [p-listeria > 0 and p-cleanable?]
    [ ifelse random 100 < prob-cleanable
      [set p-listeria (random-binomial (p-listeria) (reduction))
      set p-listeria-concentration (p-listeria / [p-area] of self)]
      [set temporary-niches temporary-niches + 1
      set p-temp-niche p-temp-niche + 1]]
 ;;include for mezzanine
 ask patches with [m-floor-listeria > 0 ]
    [ ifelse random 100 < prob-cleanable
      [set m-floor-listeria (random-binomial (m-floor-listeria) (reduction))
      set m-floor-listeria-concentration (m-floor-listeria / [p-area] of self)]
      [set temporary-niches temporary-niches + 1
      set p-temp-niche p-temp-niche + 1]]
end

to sample ; will need modification depending on historical sampling practice and data to be collected
  set sample-time ticks
  ask zones [set z-sampled? false]
   ask possible-z-sites  [
        print who
        if z-category = 1 [set z1-samples z1-samples + 1]
        if z-category = 2 [set z2-samples z2-samples + 1]
        if z-category = 3 [set z3-samples z3-samples + 1]
        set daily-samples daily-samples + 1
        if production-time = 0 [set beg-samples beg-samples + 1]
        if production-time = 3 [set mid-samples mid-samples + 1]
        if production-time > 5 [set post-samples post-samples + 1]
        if z-equipment = "belt" [set belt-samples belt-samples + 1]
        if z-equipment = "cart" [set cart-samples cart-samples + 1]
        if z-equipment = "control-panel" [set control-panel-samples control-panel-samples + 1]
        if z-equipment = "door" [set door-samples door-samples + 1]
        if z-equipment = "drain" [set drain-samples drain-samples + 1]
        if z-equipment = "employee" [set employee-samples employee-samples + 1]
        ;if z-equipment = "floors" [set floors-samples floors-samples + 1]
        if z-equipment = "flume" [set flume-samples flume-samples + 1]
        if z-equipment = "frame" [set frame-samples frame-samples + 1]
        if z-equipment = "frame-in" [set frame-in-samples frame-in-samples + 1]
        if z-equipment = "ladder" [set ladder-samples ladder-samples + 1]
        if z-equipment = "misc" [set misc-samples misc-samples + 1]
        if z-equipment = "packing" [set packing-samples packing-samples + 1]
        if z-equipment = "squeegee" [set squeegee-samples squeegee-samples + 1]
        if z-equipment = "trash-gray" [set trash-gray-samples trash-gray-samples + 1]
        if z-equipment = "trash-white" [set trash-white-samples trash-white-samples + 1]
        if z-equipment = "trash-yellow" [set trash-yellow-samples trash-yellow-samples + 1]
        if z-equipment = "wall" [set wall-samples wall-samples + 1]
        if z-equipment = "weigher" [set weigher-samples weigher-samples + 1]
        if z-location = "a-b" [set a-b-samples a-b-samples + 1]
        if z-location = "a-flume" [set a-flume-samples a-flume-samples + 1]
        if z-location = "a-line" [set a-line-samples a-line-samples + 1]
        if z-location = "a-under-room" [set a-under-room-samples a-under-room-samples + 1]
        if z-location = "b-flume" [set b-flume-samples b-flume-samples + 1]
        if z-location = "b-line" [set b-line-samples b-line-samples + 1]
        if z-location = "coring" [set coring-samples coring-samples + 1]
        if z-location = "elevator" [set elevator-samples elevator-samples + 1]
        if z-location = "hopper" [set hopper-samples hopper-samples + 1]
        if z-location = "packing-" [set packing--samples packing--samples + 1]
        if z-location = "scoring" [set scoring-samples scoring-samples + 1]
        if z-location = "sorting" [set sorting-samples sorting-samples + 1]
        if z-location = "spin-dryer" [set spin-dryer-samples spin-dryer-samples + 1]
        if z-location = "weigher-platform" [set weigher-platform-samples weigher-platform-samples + 1]
        if z-room = "A" [set A-samples A-samples + 1]
        if z-room = "B" [set B-samples B-samples + 1]
        set z-sampled? true set color lime
        let result 0

        if (z-listeria-concentration > 0 and z-listeria-concentration <= 10)
           [ ifelse random-float 100 < 10 ; when concentration between 0-10 CFU/cm^2, 10% FN rate
             [set result -1][set result 1]]
        if (z-listeria-concentration > 10 and z-listeria-concentration <= 100)
            [ ifelse random-float 100 < 1 ; when concentration between 10-100 CFU/cm^2, 1% FN rate
                [set result -1][set result 1]]
        if (z-listeria-concentration > 100)
                [set result 1]

        ;second contingency table
        if result = 1 [set a-table2 a-table2 + 1]
        print (word "a2:" a-table2)
        set b-table2 0
        print (word "b2:" b-table2)
        if result = -1 [set c-table2 c-table2 + 1]
        print (word "c2:" c-table2)
        if result = 0 [set d-table2 d-table2 + 1]
        print (word "d2:" d-table2)


        if result = 1
          [set z-sample-prev z-sample-prev + 1
           if z-category = 1 [set z1-prev z1-prev + 1]
           if z-category = 2 [set z2-prev z2-prev + 1]
           if z-category = 3 [set z3-prev z3-prev + 1]
           set daily-prev daily-prev + 1
           if production-time = 0 [set z-sample-beg-prev z-sample-beg-prev + 1]
           if production-time = 3 [set z-sample-mid-prev z-sample-mid-prev + 1]
           if production-time > 5 [set z-sample-post-prev z-sample-post-prev + 1]
           if z-equipment = "belt" [set belt-prev belt-prev + 1]
           if z-equipment = "cart" [set cart-prev cart-prev + 1]
           if z-equipment = "control-panel" [set control-panel-prev control-panel-prev + 1]
           if z-equipment = "door" [set door-prev door-prev + 1]
           if z-equipment = "drain" [set drain-prev drain-prev + 1]
           if z-equipment = "employee" [set employee-prev employee-prev + 1]
           ;if z-equipment = "floors" [set floors-prev floors-prev + 1]
           if z-equipment = "flume" [set flume-prev flume-prev + 1]
           if z-equipment = "frame" [set frame-prev frame-prev + 1]
           if z-equipment = "frame-in" [set frame-in-prev frame-in-prev + 1]
           if z-equipment = "ladder" [set ladder-prev ladder-prev + 1]
           if z-equipment = "misc" [set misc-prev misc-prev + 1]
           if z-equipment = "packing" [set packing-prev packing-prev + 1]
           if z-equipment = "squeegee" [set squeegee-prev squeegee-prev + 1]
           if z-equipment = "trash-gray" [set trash-gray-prev trash-gray-prev + 1]
           if z-equipment = "trash-white" [set trash-white-prev trash-white-prev + 1]
           if z-equipment = "trash-yellow" [set trash-yellow-prev trash-yellow-prev + 1]
           if z-equipment = "wall" [set wall-prev wall-prev + 1]
           if z-equipment = "weigher" [set weigher-prev weigher-prev + 1]
           if z-location = "a-b" [set a-b-prev a-b-prev + 1]
           if z-location = "a-flume" [set a-flume-prev a-flume-prev + 1]
           if z-location = "a-line" [set a-line-prev a-line-prev + 1]
           if z-location = "a-under-room" [set a-under-room-prev a-under-room-prev + 1]
           if z-location = "b-flume" [set b-flume-prev b-flume-prev + 1]
           if z-location = "b-line" [set b-line-prev b-line-prev + 1]
           if z-location = "coring" [set coring-prev coring-prev + 1]
           if z-location = "elevator" [set elevator-prev elevator-prev + 1]
           if z-location = "hopper" [set hopper-prev hopper-prev + 1]
           if z-location = "packing-" [set packing--prev packing--prev + 1]
           if z-location = "scoring" [set scoring-prev scoring-prev + 1]
           if z-location = "sorting" [set sorting-prev sorting-prev + 1]
           if z-location = "spin-dryer" [set spin-dryer-prev spin-dryer-prev + 1]
           if z-location = "weigher-platform" [set weigher-platform-prev weigher-platform-prev + 1]
           if z-room = "A" [set A-prev A-prev + 1]
           if z-room = "B" [set B-prev B-prev + 1]
          ]
    set possible-z-sites possible-z-sites with [self != myself]
    ;ask possible-z-sites [print self]
  ]


 ; ask possible-p-sites  [
 ;    	   set floor-samples floor-samples + 1
 ;          set p-sampled? true
 ;          let result 0
;	ifelse m-height = 0 [
 ;          if (p-listeria-concentration > 0 and p-listeria-concentration <= 10)
 ;          [ ifelse random-float 100 < 10 ; when concentration between 0-10 CFU/cm^2, 10% FN rate
 ;            [set result -1][set result 1]]
 ;          if (p-listeria-concentration > 10 and p-listeria-concentration <= 100)
 ;           [ ifelse random-float 100 < 1 ; when concentration between 10-100 CFU/cm^2, 1% FN rate
 ;               [set result -1][set result 1]]
 ;          if (p-listeria-concentration > 100)
 ;               [set result 1]]
;	[  if (m-listeria-concentration > 0 and m-listeria-concentration <= 10)
 ;          [ ifelse random-float 100 < 10 ; when concentration between 0-10 CFU/cm^2, 10% FN rate
 ;            [set result -1][set result 1]]
 ;          if (m-listeria-concentration > 10 and m-listeria-concentration <= 100)
 ;           [ ifelse random-float 100 < 1 ; when concentration between 10-100 CFU/cm^2, 1% FN rate
 ;               [set result -1][set result 1]]
 ;          if (m-listeria-concentration > 100)
 ;               [set result 1]]

 ;       if result = 1
 ;        [set z-sample-prev z-sample-prev + 1
 ;         set floors-prev floors-prev + 1]

 ;   set possible-p-sites possible-p-sites with [self != myself]
 ;   ]

  ;first contingency table
  set start-sampling? false
  set a-table (count zones with [z-sampled? and z-listeria-concentration > 0])
  print (word "a:" a-table)
  set b-table (count zones with [z-sampled? and z-listeria-concentration = 0])
  print (word "b:" b-table)
  set c-table (count zones with [not z-sampled? and z-listeria-concentration > 0])
  print (word "c:" c-table)
  set d-table (count zones with [not z-sampled? and z-listeria-concentration = 0])
  print (word "d:" d-table)
end

;these are options from the chooser on interface for sampling-sites - still experimenting here but may not be relevant to your experiments for corrective actions

to Facility-G-sites
  let b random-in-range 14 21
  ;print b
  set possible-z-sites (turtle-set a-zone 34 a-zone 35 a-zone 67 a-zone 77 a-zone 93 a-zone 106 a-zone 110 a-zone 118 a-zone 119 a-zone 120 a-zone 121 a-zone 122
    a-zone 123 a-zone 125 a-zone 127 a-zone 129 a-zone 139 a-zone 155 a-zone 157 a-zone 159 a-zone 161 a-zone 163 a-zone 165 a-zone 167 a-zone 180 a-zone 194 a-zone 196
    a-zone 210 a-zone 242 a-zone 260 a-zone 261 a-zone 262 a-zone 263 a-zone 264 a-zone 265 a-zone 267 a-zone 271 a-zone 273 a-zone 274 a-zone 277 a-zone 281 a-zone 285
    a-zone 287 a-zone 289 a-zone 290 )
  set possible-z-sites n-of b possible-z-sites

  let f random-in-range 4 10
  set possible-z-patches (patch-set patch 53 26 patch 14 39 patch 20 32 patch 20 32 patch 1	10 patch 54	27 patch 92	34 patch 70	1 patch 1	26 patch 52	39 patch 51	2 patch 63 39 patch 69 27 patch 54 30)
  set possible-z-patches n-of f possible-z-patches
end

to Scenario-1
  set zone1-sample 0
  set zone2-sample 10
  set zone3-sample 20 ; change to 20 if not sampling floor
  set possible-z-sites (turtle-set (n-of zone1-sample zones with [ z-category = 1])
  (n-of zone2-sample zones with [ z-category = 2]) (n-of zone3-sample zones with [z-category = 3]))
  ;set possible-p-sites n-of 5 patches with [pcolor != 0 ] ;comment out this line if not modeling floor
end

to Scenario-2
  set zone1-sample 10
  set zone2-sample 10
  set zone3-sample 10 ; change to 10 if not sampling floor
  set possible-z-sites (turtle-set (n-of zone1-sample zones with [ z-category = 1])
  (n-of zone2-sample zones with [ z-category = 2]) (n-of zone3-sample zones with [z-category = 3]))
  ;set possible-p-sites n-of 5 patches with [pcolor != 0 ] ;comment out this line if not modeling floor
end

to Scenario-3
  set possible-z-sites (turtle-set (n-of 30 zones))
  ;set possible-p-sites n-of 5 patches with [pcolor != 0]
end

to Scenario-4
  set zone1-sample 0
  set zone2-sample 0
  set zone3-sample 30
  set possible-z-sites (turtle-set (n-of zone1-sample zones with [ z-category = 1])
  (n-of zone2-sample zones with [ z-category = 2]) (n-of zone3-sample zones with [z-category = 3]) )
  ;set possible-p-sites n-of 5 patches with [pcolor != 0 ]
end




;;; these are reporter functions used throughout the code - mainly probability distributions - don't change
to-report random-binomial [n p]
     report sum n-values n [ifelse-value (p > random-float 1) [1] [0]]
end

to-report random-in-range [low high]
  report low + random (high - low + 1)
end

to-report randomfloat-in-range [low high]
  report low + random-float (high - low)
end

to-report random-pert [minval likeval maxval lambda]  ;;taken from internet at: http://stackoverflow.com/questions/30807377/netlogo-sampling-from-a-beta-pert-distribtuion
  ;use pert params to draw from a beta distribution
  if not (minval <= likeval and likeval <= maxval) [error "wrong argument ranking"]
  if (minval = likeval and likeval = maxval) [report minval] ;;handle trivial inputs
  let pert-var 1 / 36
  let pert-mean (maxval + lambda * likeval - 5 * minval) / (6 * (maxval - minval))
  let temp pert-mean * (1 - pert-mean) / pert-var
  let alpha1 pert-mean * (temp - 1)
  let alpha2 (1 - pert-mean) * (temp - 1)
  let x1 random-gamma alpha1 1
  let x2 random-gamma alpha2 1
  report (x1 / (x1 + x2)) * (maxval - minval) + minval
end

;;; these are functions to update colors on interface - can change if want a different color scheme
to p-water-recolor
  if (p-water != 0) and (p-water != 5) and (p-water != 35) [set pcolor (blue + 4.9 - p-water)]
  if (p-water <= 1) [set pcolor white]
  if (p-water = 0) [set pcolor black]
  if (p-water = 5) [set pcolor gray]
  if (p-water = 6) [set pcolor gray]
  if (p-water = 35) [set pcolor brown]
end

to p-traffic-recolor
  if (p-traffic != 0) and (p-water != 5) and (p-water != 35) [set pcolor (green + 4.9 - p-traffic)]
  if (p-traffic = 0) [set pcolor black]
  if (p-traffic = 1) [set pcolor white]
  if (p-traffic = 5) [set pcolor gray]
  if (p-traffic = 6) [set pcolor gray]
  if (p-traffic = 35) [set pcolor brown]
end

to zone-recolor
  if z-category = 1 [set color (red)]
  if z-category = 2 [set color (orange)]
  if z-category = 3 [ set color (magenta)]
  if z-item-name = "employee" [set color (black)]
end
@#$#@#$#@
GRAPHICS-WINDOW
274
51
1231
474
-1
-1
10.1
1
10
1
1
1
0
0
0
1
0
93
0
40
1
1
1
ticks
30.0

SLIDER
25
92
257
125
initial-zone-prevalence
initial-zone-prevalence
0
10
0.0
1
1
%
HORIZONTAL

BUTTON
30
10
104
49
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
123
11
195
50
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
18
398
75
443
zone 1
count zones with [(z-category = 1) and (z-listeria-concentration > 0)]
17
1
11

MONITOR
78
399
134
444
zone 2
count zones with [(z-category = 2 ) and (z-listeria-concentration > 0)]
0
1
11

MONITOR
139
399
195
444
zone 3
count zones with [(z-category = 3) and (z-listeria-concentration > 0)]
0
1
11

SLIDER
24
130
257
163
initial-environment-prevalence
initial-environment-prevalence
0
10
0.0
1
1
%
HORIZONTAL

PLOT
1234
285
1592
548
Facility Listeria Population
time (h)
listeria-density [log(cfu)]
0.0
336.0
0.0
3.0
true
true
"" ""
PENS
"zone-1 " 1.0 0 -10899396 true "" "plot log(mean [z-listeria-concentration + 1e-1] of zones with [z-category = 1]) 10"
"zone-2" 1.0 0 -955883 true "" "plot log (mean [z-listeria-concentration + 1e-1] of zones with [z-category = 2]) 10"
"zone-3" 1.0 0 -5825686 true "" "plot log ( mean [z-listeria-concentration + 1e-1] of zones with [z-category = 3]) 10"
"patches" 1.0 0 -13345367 true "" "plot log( mean [p-listeria-concentration + 1e-1] of patches ) 10"

SLIDER
24
169
258
202
initial-contamination-level
initial-contamination-level
0.3
8
1.0
.1
1
log CFU/cm^2
HORIZONTAL

CHOOSER
131
207
256
252
environment-view
environment-view
"water" "traffic"
1

MONITOR
1359
552
1482
593
mean patch CFU/sq. cm
(mean [p-listeria-concentration] of patches with [p-water != 0])
2
1
10

MONITOR
1342
596
1452
641
cleanable + contam.
count zones with [(z-cleanable?) and (z-listeria-concentration > 0)]
0
1
11

MONITOR
1456
595
1587
640
not cleanable +  contam.
count zones with [(not z-cleanable?) and (z-listeria-concentration > 0)]
0
1
11

OUTPUT
914
11
1057
43
14

MONITOR
1236
596
1338
641
contam. sites
count zones with [z-listeria-concentration > 0 ]
0
1
11

MONITOR
1235
552
1352
593
mean zone CFU/sq. cm
mean [z-listeria-concentration] of zones with [z-listeria-concentration > 0]
2
1
10

PLOT
1233
10
1587
281
Prevalence
Time (hr)
Contaminated Sites
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"zones" 1.0 1 -16777216 true "" "plot (count zones with [z-listeria-concentration > 0 ] / count zones)"
"patches" 1.0 1 -2674135 true "" "plot (count patches with [p-listeria-concentration > 0 ] / ((count patches with [p-water != 0]) + 1))"

MONITOR
601
606
681
651
zone-to-zone
zone-to-zone
17
1
11

MONITOR
343
557
406
602
food falls
zone-to-patch
0
1
11

MONITOR
410
556
508
601
condensation-drip
condensation-drip
0
1
11

MONITOR
698
557
767
602
cont. fruit
introduction-food
0
1
11

MONITOR
436
605
506
650
introduction
introduction-zone4
0
1
11

MONITOR
4
448
74
493
employees
count zones with [z-item-name = \"employee\" and z-listeria-concentration > 0]
17
1
11

MONITOR
68
350
125
395
Z(dry)
count zones with [z-water <= 1]
17
1
11

MONITOR
133
349
191
394
Z(moist)
count zones with [z-water > 1 and z-water <= 2]
17
1
11

MONITOR
198
349
255
394
Z(wet)
count zones with [z-water > 2 and z-water <= 3]
17
1
11

MONITOR
511
556
597
601
random events
random-event
0
1
11

MONITOR
1080
10
1219
55
current-event
event
17
1
11

MONITOR
137
496
194
541
drains
count zones with [(z-equipment = \"drain\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
197
448
268
493
Room A
count patches with [ (pxcor < 52) and (p-listeria-concentration > 0 )]
17
1
11

MONITOR
71
546
147
591
sorting
count zones with [ (z-equipment = \"Sorting\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
10
496
68
541
scoring
count zones with [ (z-equipment = \"scoring\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
152
546
226
591
packing
count zones with [ (z-equipment = \"Packing\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
5
594
68
639
A Flumes
count zones with [ (z-equipment = \"a-flume\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
72
595
134
640
B Flumes
count zones with [ (z-equipment = \"b-flume\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
143
643
205
688
mezzanine
count patches with [(pxcor > 56 and pxcor < 83 and pycor > 27 and pycor < 38) and (m-floor-listeria > 0)]
0
1
11

MONITOR
230
546
317
591
QA Employee
count zones with [ (z-equipment = \"QA\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
208
643
314
688
weigher platform
count zones with [ (z-equipment = \"weigher-platform\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
72
643
139
688
pack floor
count patches with [(pxcor > 54 and pxcor < 92 and pycor < 27) and (p-listeria > 0)]
0
1
11

MONITOR
78
448
132
493
ceiling 
count patches with [c-listeria-concentration > 0]
0
1
11

MONITOR
72
496
132
541
coring
count zones with [ (z-equipment = \"coring\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
134
448
194
493
Room B
count patches with [(pxcor > 52) and (p-listeria > 0)]
0
1
11

MONITOR
769
557
860
602
cont. fruit CFUs
sum-fruit-load
0
1
11

MONITOR
864
556
972
601
fruit CFU transfered
sum-fruit-transfer
0
1
11

MONITOR
684
606
752
651
z1-contacts
mean [z-contacts] of zones with [z-category = 1]
2
1
11

MONITOR
755
606
828
651
z1-transfers
mean [z-transfers] of zones with [z-category = 1]
2
1
11

MONITOR
830
606
920
651
z1-zone4 intros
mean [z-zone4] of zones with [z-category = 1]
2
1
11

MONITOR
601
557
695
602
z1 -random intro
mean [z-random] of zones with [z-category = 1]
2
1
11

MONITOR
510
606
598
651
NIL
patch-to-patch
0
1
11

SLIDER
25
53
257
86
time-of-simulation
time-of-simulation
1
12
2.0
1
1
weeks
HORIZONTAL

MONITOR
1484
552
1607
593
mean ceiling CFU/sq. cm
(mean [c-listeria-concentration] of patches with [p-water != 0 and p-water != 35])
2
1
10

CHOOSER
24
207
127
252
detection-limit
detection-limit
1 10 100
0

MONITOR
376
606
433
651
z1-fruit
mean [z-fruit] of zones with [z-category = 1]
2
1
11

SLIDER
25
255
255
288
numSamples
numSamples
0
100
10.0
1
1
samples/wk
HORIZONTAL

CHOOSER
24
296
148
341
sampling-time
sampling-time
"none" "first-hour" "mid-production" "second-shift"
3

CHOOSER
151
295
255
340
sampling-sites
sampling-sites
"Facility-G-sites" "Scenario-1" "Scenario-2" "Scenario-3" "Scenario-4"
1

MONITOR
138
595
203
640
hoppers
count zones with [ (z-equipment = \"hopper\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
206
595
273
640
elevators
count zones with [(z-equipment = \"elevator\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
198
495
299
540
cutting platform
count zones with [(z-equipment = \"cutting-platform\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
3
546
66
591
hand tools
count zones with [(z-item-name = \"hammer\" or z-item-name = \"scoop\" or z-item-name = \"white-squeegee\") and (z-listeria-concentration > 0)]
0
1
11

MONITOR
5
643
67
688
AB Wall
count zones with [(z-equipment = \"a-b\") and (z-listeria-concentration > 0)]
0
1
11

SLIDER
1034
553
1206
586
seed-me
seed-me
0
100
1000.0
1
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="NT validation" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>spread-tc</metric>
    <metric>m-crate-prevalence</metric>
    <metric>t-crate-prevalence</metric>
    <metric>w-crate-prevalence</metric>
    <metric>r-crate-prevalence</metric>
    <metric>f-crate-prevalence</metric>
    <metric>fruit-conc</metric>
    <metric>fruit-tc</metric>
    <metric>zone4-load</metric>
    <metric>random-load</metric>
    <metric>m-san-reduction</metric>
    <metric>t-san-reduction</metric>
    <metric>w-san-reduction</metric>
    <metric>r-san-reduction</metric>
    <metric>f-san-reduction</metric>
    <metric>prob-cleanable</metric>
    <metric>mu-max</metric>
    <metric>K</metric>
    <metric>p-patch-spread</metric>
    <metric>p-water-spread</metric>
    <metric>p-food-dropped</metric>
    <metric>p-condensation</metric>
    <metric>p-random-noise</metric>
    <metric>p-zone4-intro</metric>
    <metric>p11</metric>
    <metric>p12</metric>
    <metric>p13</metric>
    <metric>p14</metric>
    <metric>p21</metric>
    <metric>p22</metric>
    <metric>p23</metric>
    <metric>p24</metric>
    <metric>p31</metric>
    <metric>p32</metric>
    <metric>p33</metric>
    <metric>p34</metric>
    <metric>p41</metric>
    <metric>p42</metric>
    <metric>p43</metric>
    <metric>p44</metric>
    <metric>tc11</metric>
    <metric>tc12</metric>
    <metric>tc13</metric>
    <metric>tc14</metric>
    <metric>tc21</metric>
    <metric>tc22</metric>
    <metric>tc23</metric>
    <metric>tc24</metric>
    <metric>tc31</metric>
    <metric>tc32</metric>
    <metric>tc33</metric>
    <metric>tc34</metric>
    <metric>tc41</metric>
    <metric>tc42</metric>
    <metric>tc43</metric>
    <metric>tc44</metric>
    <metric>z1-initial</metric>
    <metric>z2-initial</metric>
    <metric>z3-initial</metric>
    <metric>p-initial</metric>
    <metric>e-initial</metric>
    <metric>m-zones-pre</metric>
    <metric>t-zones-pre</metric>
    <metric>w-zones-pre</metric>
    <metric>r-zones-pre</metric>
    <metric>f-zones-pre</metric>
    <metric>m-zones-beg</metric>
    <metric>t-zones-beg</metric>
    <metric>w-zones-beg</metric>
    <metric>r-zones-beg</metric>
    <metric>f-zones-beg</metric>
    <metric>m-zones-mid</metric>
    <metric>t-zones-mid</metric>
    <metric>w-zones-mid</metric>
    <metric>r-zones-mid</metric>
    <metric>f-zones-mid</metric>
    <metric>m-zones-post</metric>
    <metric>t-zones-post</metric>
    <metric>w-zones-post</metric>
    <metric>r-zones-post</metric>
    <metric>f-zones-post</metric>
    <metric>zone-to-zone</metric>
    <metric>avg-listeria-transferred</metric>
    <metric>zone-to-patch</metric>
    <metric>condensation-drip</metric>
    <metric>introduction-zone4</metric>
    <metric>patch-to-patch</metric>
    <metric>introduction-food</metric>
    <metric>sum-fruit-load</metric>
    <metric>sum-fruit-transfer</metric>
    <metric>random-event</metric>
    <metric>employee-FCS</metric>
    <metric>NFCS-FCS</metric>
    <metric>total-contam-events</metric>
    <metric>m-random-event</metric>
    <metric>t-random-event</metric>
    <metric>w-random-event</metric>
    <metric>r-random-event</metric>
    <metric>f-random-event</metric>
    <metric>median-z1-time</metric>
    <metric>median-z2-time</metric>
    <metric>median-z3-time</metric>
    <metric>median-e-time</metric>
    <metric>median-p-time</metric>
    <metric>median-fwj-time</metric>
    <metric>median-z1-max-contam-bout</metric>
    <metric>median-z2-max-contam-bout</metric>
    <metric>median-z3-max-contam-bout</metric>
    <metric>median-e-max-contam-bout</metric>
    <metric>median-p-max-contam-bout</metric>
    <metric>median-fwj-max-contam-bout</metric>
    <metric>temporary-niches</metric>
    <metric>p-temp-niches</metric>
    <metric>z1-temp-niches</metric>
    <metric>z2-temp-niches</metric>
    <metric>z3-temp-niches</metric>
    <metric>z1-contacts</metric>
    <metric>z1-transfers</metric>
    <metric>z2-contacts</metric>
    <metric>z2-transfers</metric>
    <metric>z3-contacts</metric>
    <metric>z3-transfers</metric>
    <metric>z1-prev</metric>
    <metric>z2-prev</metric>
    <metric>z3-prev</metric>
    <metric>concentration-list</metric>
    <metric>time-contaminated-list</metric>
    <metric>max-contam-bout-list</metric>
    <metric>contacts-list</metric>
    <metric>transfers-list</metric>
    <metric>temp-niche-list</metric>
    <metric>niches-estab-list</metric>
    <metric>fruit-cont-list</metric>
    <metric>random-cont-list</metric>
    <metric>belt-prev</metric>
    <metric>cart-prev</metric>
    <metric>control-panel-prev</metric>
    <metric>door-prev</metric>
    <metric>drain-prev</metric>
    <metric>employee-prev</metric>
    <metric>flume-prev</metric>
    <metric>frame-prev</metric>
    <metric>frame-in-prev</metric>
    <metric>ladder-prev</metric>
    <metric>misc-prev</metric>
    <metric>packing-prev</metric>
    <metric>squeegee-prev</metric>
    <metric>trash-gray-prev</metric>
    <metric>trash-white-prev</metric>
    <metric>trash-yellow-prev</metric>
    <metric>wall-prev</metric>
    <metric>weigher-prev</metric>
    <metric>a-b-prev</metric>
    <metric>a-flume-prev</metric>
    <metric>a-line-prev</metric>
    <metric>a-under-room-prev</metric>
    <metric>b-flume-prev</metric>
    <metric>b-line-prev</metric>
    <metric>coring-prev</metric>
    <metric>elevator-prev</metric>
    <metric>hopper-prev</metric>
    <metric>packing--prev</metric>
    <metric>scoring-prev</metric>
    <metric>sorting-prev</metric>
    <metric>spin-dryer-prev</metric>
    <metric>weigher-platform-prev</metric>
    <metric>A-prev</metric>
    <metric>B-prev</metric>
    <metric>floors-prev</metric>
    <metric>random-floor</metric>
    <metric>random-ceiling</metric>
    <metric>random-agent</metric>
    <enumeratedValueSet variable="random-seed">
      <value value="4519"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="site_specific_validation" repetitions="10000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>spread-tc</metric>
    <metric>m-crate-prevalence</metric>
    <metric>t-crate-prevalence</metric>
    <metric>w-crate-prevalence</metric>
    <metric>r-crate-prevalence</metric>
    <metric>f-crate-prevalence</metric>
    <metric>fruit-conc</metric>
    <metric>fruit-tc</metric>
    <metric>zone4-load</metric>
    <metric>random-load</metric>
    <metric>m-san-reduction</metric>
    <metric>t-san-reduction</metric>
    <metric>w-san-reduction</metric>
    <metric>r-san-reduction</metric>
    <metric>f-san-reduction</metric>
    <metric>prob-cleanable</metric>
    <metric>mu-max</metric>
    <metric>K</metric>
    <metric>p-patch-spread</metric>
    <metric>p-water-spread</metric>
    <metric>p-food-dropped</metric>
    <metric>p-condensation</metric>
    <metric>p-random-noise</metric>
    <metric>p-zone4-intro</metric>
    <metric>p11</metric>
    <metric>p12</metric>
    <metric>p13</metric>
    <metric>p14</metric>
    <metric>p21</metric>
    <metric>p22</metric>
    <metric>p23</metric>
    <metric>p24</metric>
    <metric>p31</metric>
    <metric>p32</metric>
    <metric>p33</metric>
    <metric>p34</metric>
    <metric>p41</metric>
    <metric>p42</metric>
    <metric>p43</metric>
    <metric>p44</metric>
    <metric>tc11</metric>
    <metric>tc12</metric>
    <metric>tc13</metric>
    <metric>tc14</metric>
    <metric>tc21</metric>
    <metric>tc22</metric>
    <metric>tc23</metric>
    <metric>tc24</metric>
    <metric>tc31</metric>
    <metric>tc32</metric>
    <metric>tc33</metric>
    <metric>tc34</metric>
    <metric>tc41</metric>
    <metric>tc42</metric>
    <metric>tc43</metric>
    <metric>tc44</metric>
    <metric>z1-initial</metric>
    <metric>z2-initial</metric>
    <metric>z3-initial</metric>
    <metric>p-initial</metric>
    <metric>e-initial</metric>
    <metric>m-zones-pre</metric>
    <metric>t-zones-pre</metric>
    <metric>w-zones-pre</metric>
    <metric>r-zones-pre</metric>
    <metric>f-zones-pre</metric>
    <metric>m-zones-beg</metric>
    <metric>t-zones-beg</metric>
    <metric>w-zones-beg</metric>
    <metric>r-zones-beg</metric>
    <metric>f-zones-beg</metric>
    <metric>m-zones-mid</metric>
    <metric>t-zones-mid</metric>
    <metric>w-zones-mid</metric>
    <metric>r-zones-mid</metric>
    <metric>f-zones-mid</metric>
    <metric>m-zones-post</metric>
    <metric>t-zones-post</metric>
    <metric>w-zones-post</metric>
    <metric>r-zones-post</metric>
    <metric>f-zones-post</metric>
    <metric>zone-to-zone</metric>
    <metric>avg-listeria-transferred</metric>
    <metric>zone-to-patch</metric>
    <metric>condensation-drip</metric>
    <metric>introduction-zone4</metric>
    <metric>patch-to-patch</metric>
    <metric>introduction-food</metric>
    <metric>sum-fruit-load</metric>
    <metric>sum-fruit-transfer</metric>
    <metric>random-event</metric>
    <metric>employee-FCS</metric>
    <metric>NFCS-FCS</metric>
    <metric>total-contam-events</metric>
    <metric>m-random-event</metric>
    <metric>t-random-event</metric>
    <metric>w-random-event</metric>
    <metric>r-random-event</metric>
    <metric>f-random-event</metric>
    <metric>median-z1-time</metric>
    <metric>median-z2-time</metric>
    <metric>median-z3-time</metric>
    <metric>median-e-time</metric>
    <metric>median-p-time</metric>
    <metric>median-fwj-time</metric>
    <metric>median-z1-max-contam-bout</metric>
    <metric>median-z2-max-contam-bout</metric>
    <metric>median-z3-max-contam-bout</metric>
    <metric>median-e-max-contam-bout</metric>
    <metric>median-p-max-contam-bout</metric>
    <metric>median-fwj-max-contam-bout</metric>
    <metric>temporary-niches</metric>
    <metric>p-temp-niches</metric>
    <metric>z1-temp-niches</metric>
    <metric>z2-temp-niches</metric>
    <metric>z3-temp-niches</metric>
    <metric>z1-contacts</metric>
    <metric>z1-transfers</metric>
    <metric>z2-contacts</metric>
    <metric>z2-transfers</metric>
    <metric>z3-contacts</metric>
    <metric>z3-transfers</metric>
    <metric>z1-prev</metric>
    <metric>z2-prev</metric>
    <metric>z3-prev</metric>
    <metric>concentration-list</metric>
    <metric>time-contaminated-list</metric>
    <metric>max-contam-bout-list</metric>
    <metric>contacts-list</metric>
    <metric>transfers-list</metric>
    <metric>temp-niche-list</metric>
    <metric>niches-estab-list</metric>
    <metric>fruit-cont-list</metric>
    <metric>random-cont-list</metric>
    <metric>belt-prev</metric>
    <metric>cart-prev</metric>
    <metric>control-panel-prev</metric>
    <metric>door-prev</metric>
    <metric>drain-prev</metric>
    <metric>employee-prev</metric>
    <metric>flume-prev</metric>
    <metric>frame-prev</metric>
    <metric>frame-in-prev</metric>
    <metric>ladder-prev</metric>
    <metric>misc-prev</metric>
    <metric>packing-prev</metric>
    <metric>squeegee-prev</metric>
    <metric>trash-gray-prev</metric>
    <metric>trash-white-prev</metric>
    <metric>trash-yellow-prev</metric>
    <metric>wall-prev</metric>
    <metric>weigher-prev</metric>
    <metric>a-b-prev</metric>
    <metric>a-flume-prev</metric>
    <metric>a-line-prev</metric>
    <metric>a-under-room-prev</metric>
    <metric>b-flume-prev</metric>
    <metric>b-line-prev</metric>
    <metric>coring-prev</metric>
    <metric>elevator-prev</metric>
    <metric>hopper-prev</metric>
    <metric>packing--prev</metric>
    <metric>scoring-prev</metric>
    <metric>sorting-prev</metric>
    <metric>spin-dryer-prev</metric>
    <metric>weigher-platform-prev</metric>
    <metric>A-prev</metric>
    <metric>B-prev</metric>
    <metric>floors-prev</metric>
    <metric>random-floor</metric>
    <metric>random-ceiling</metric>
    <metric>random-agent</metric>
    <metric>m-zone-34b</metric>
    <metric>m-zone-35b</metric>
    <metric>m-zone-67b</metric>
    <metric>m-zone-69b</metric>
    <metric>m-zone-71b</metric>
    <metric>m-zone-73b</metric>
    <metric>m-zone-75b</metric>
    <metric>m-zone-77b</metric>
    <metric>m-zone-79b</metric>
    <metric>m-zone-81b</metric>
    <metric>m-zone-83b</metric>
    <metric>m-zone-85b</metric>
    <metric>m-zone-86b</metric>
    <metric>m-zone-87b</metric>
    <metric>m-zone-88b</metric>
    <metric>m-zone-89b</metric>
    <metric>m-zone-90b</metric>
    <metric>m-zone-91b</metric>
    <metric>m-zone-92b</metric>
    <metric>m-zone-93b</metric>
    <metric>m-zone-94b</metric>
    <metric>m-zone-95b</metric>
    <metric>m-zone-96b</metric>
    <metric>m-zone-97b</metric>
    <metric>m-zone-98b</metric>
    <metric>m-zone-99b</metric>
    <metric>m-zone-106b</metric>
    <metric>m-zone-107b</metric>
    <metric>m-zone-108b</metric>
    <metric>m-zone-109b</metric>
    <metric>m-zone-110b</metric>
    <metric>m-zone-111b</metric>
    <metric>m-zone-112b</metric>
    <metric>m-zone-113b</metric>
    <metric>m-zone-118b</metric>
    <metric>m-zone-119b</metric>
    <metric>m-zone-120b</metric>
    <metric>m-zone-121b</metric>
    <metric>m-zone-122b</metric>
    <metric>m-zone-123b</metric>
    <metric>m-zone-125b</metric>
    <metric>m-zone-127b</metric>
    <metric>m-zone-129b</metric>
    <metric>m-zone-139b</metric>
    <metric>m-zone-155b</metric>
    <metric>m-zone-157b</metric>
    <metric>m-zone-159b</metric>
    <metric>m-zone-161b</metric>
    <metric>m-zone-163b</metric>
    <metric>m-zone-165b</metric>
    <metric>m-zone-167b</metric>
    <metric>m-zone-180b</metric>
    <metric>m-zone-194b</metric>
    <metric>m-zone-196b</metric>
    <metric>m-zone-210b</metric>
    <metric>m-zone-242b</metric>
    <metric>m-zone-259b</metric>
    <metric>m-zone-260b</metric>
    <metric>m-zone-261b</metric>
    <metric>m-zone-262b</metric>
    <metric>m-zone-263b</metric>
    <metric>m-zone-264b</metric>
    <metric>m-zone-265b</metric>
    <metric>m-zone-266b</metric>
    <metric>m-zone-267b</metric>
    <metric>m-zone-268b</metric>
    <metric>m-zone-269b</metric>
    <metric>m-zone-270b</metric>
    <metric>m-zone-271b</metric>
    <metric>m-zone-272b</metric>
    <metric>m-zone-273b</metric>
    <metric>m-zone-274b</metric>
    <metric>m-zone-277b</metric>
    <metric>m-zone-278b</metric>
    <metric>m-zone-279b</metric>
    <metric>m-zone-280b</metric>
    <metric>m-zone-281b</metric>
    <metric>m-zone-282b</metric>
    <metric>m-zone-283b</metric>
    <metric>m-zone-284b</metric>
    <metric>m-zone-285b</metric>
    <metric>m-zone-286b</metric>
    <metric>m-zone-287b</metric>
    <metric>m-zone-289b</metric>
    <metric>m-zone-290b</metric>
    <metric>m-zone-291b</metric>
    <metric>t-zone-34b</metric>
    <metric>t-zone-35b</metric>
    <metric>t-zone-67b</metric>
    <metric>t-zone-69b</metric>
    <metric>t-zone-71b</metric>
    <metric>t-zone-73b</metric>
    <metric>t-zone-75b</metric>
    <metric>t-zone-77b</metric>
    <metric>t-zone-79b</metric>
    <metric>t-zone-81b</metric>
    <metric>t-zone-83b</metric>
    <metric>t-zone-85b</metric>
    <metric>t-zone-86b</metric>
    <metric>t-zone-87b</metric>
    <metric>t-zone-88b</metric>
    <metric>t-zone-89b</metric>
    <metric>t-zone-90b</metric>
    <metric>t-zone-91b</metric>
    <metric>t-zone-92b</metric>
    <metric>t-zone-93b</metric>
    <metric>t-zone-94b</metric>
    <metric>t-zone-95b</metric>
    <metric>t-zone-96b</metric>
    <metric>t-zone-97b</metric>
    <metric>t-zone-98b</metric>
    <metric>t-zone-99b</metric>
    <metric>t-zone-106b</metric>
    <metric>t-zone-107b</metric>
    <metric>t-zone-108b</metric>
    <metric>t-zone-109b</metric>
    <metric>t-zone-110b</metric>
    <metric>t-zone-111b</metric>
    <metric>t-zone-112b</metric>
    <metric>t-zone-113b</metric>
    <metric>t-zone-118b</metric>
    <metric>t-zone-119b</metric>
    <metric>t-zone-120b</metric>
    <metric>t-zone-121b</metric>
    <metric>t-zone-122b</metric>
    <metric>t-zone-123b</metric>
    <metric>t-zone-125b</metric>
    <metric>t-zone-127b</metric>
    <metric>t-zone-129b</metric>
    <metric>t-zone-139b</metric>
    <metric>t-zone-155b</metric>
    <metric>t-zone-157b</metric>
    <metric>t-zone-159b</metric>
    <metric>t-zone-161b</metric>
    <metric>t-zone-163b</metric>
    <metric>t-zone-165b</metric>
    <metric>t-zone-167b</metric>
    <metric>t-zone-180b</metric>
    <metric>t-zone-194b</metric>
    <metric>t-zone-196b</metric>
    <metric>t-zone-210b</metric>
    <metric>t-zone-242b</metric>
    <metric>t-zone-259b</metric>
    <metric>t-zone-260b</metric>
    <metric>t-zone-261b</metric>
    <metric>t-zone-262b</metric>
    <metric>t-zone-263b</metric>
    <metric>t-zone-264b</metric>
    <metric>t-zone-265b</metric>
    <metric>t-zone-266b</metric>
    <metric>t-zone-267b</metric>
    <metric>t-zone-268b</metric>
    <metric>t-zone-269b</metric>
    <metric>t-zone-270b</metric>
    <metric>t-zone-271b</metric>
    <metric>t-zone-272b</metric>
    <metric>t-zone-273b</metric>
    <metric>t-zone-274b</metric>
    <metric>t-zone-277b</metric>
    <metric>t-zone-278b</metric>
    <metric>t-zone-279b</metric>
    <metric>t-zone-280b</metric>
    <metric>t-zone-281b</metric>
    <metric>t-zone-282b</metric>
    <metric>t-zone-283b</metric>
    <metric>t-zone-284b</metric>
    <metric>t-zone-285b</metric>
    <metric>t-zone-286b</metric>
    <metric>t-zone-287b</metric>
    <metric>t-zone-289b</metric>
    <metric>t-zone-290b</metric>
    <metric>t-zone-291b</metric>
    <metric>w-zone-34b</metric>
    <metric>w-zone-35b</metric>
    <metric>w-zone-67b</metric>
    <metric>w-zone-69b</metric>
    <metric>w-zone-71b</metric>
    <metric>w-zone-73b</metric>
    <metric>w-zone-75b</metric>
    <metric>w-zone-77b</metric>
    <metric>w-zone-79b</metric>
    <metric>w-zone-81b</metric>
    <metric>w-zone-83b</metric>
    <metric>w-zone-85b</metric>
    <metric>w-zone-86b</metric>
    <metric>w-zone-87b</metric>
    <metric>w-zone-88b</metric>
    <metric>w-zone-89b</metric>
    <metric>w-zone-90b</metric>
    <metric>w-zone-91b</metric>
    <metric>w-zone-92b</metric>
    <metric>w-zone-93b</metric>
    <metric>w-zone-94b</metric>
    <metric>w-zone-95b</metric>
    <metric>w-zone-96b</metric>
    <metric>w-zone-97b</metric>
    <metric>w-zone-98b</metric>
    <metric>w-zone-99b</metric>
    <metric>w-zone-106b</metric>
    <metric>w-zone-107b</metric>
    <metric>w-zone-108b</metric>
    <metric>w-zone-109b</metric>
    <metric>w-zone-110b</metric>
    <metric>w-zone-111b</metric>
    <metric>w-zone-112b</metric>
    <metric>w-zone-113b</metric>
    <metric>w-zone-118b</metric>
    <metric>w-zone-119b</metric>
    <metric>w-zone-120b</metric>
    <metric>w-zone-121b</metric>
    <metric>w-zone-122b</metric>
    <metric>w-zone-123b</metric>
    <metric>w-zone-125b</metric>
    <metric>w-zone-127b</metric>
    <metric>w-zone-129b</metric>
    <metric>w-zone-139b</metric>
    <metric>w-zone-155b</metric>
    <metric>w-zone-157b</metric>
    <metric>w-zone-159b</metric>
    <metric>w-zone-161b</metric>
    <metric>w-zone-163b</metric>
    <metric>w-zone-165b</metric>
    <metric>w-zone-167b</metric>
    <metric>w-zone-180b</metric>
    <metric>w-zone-194b</metric>
    <metric>w-zone-196b</metric>
    <metric>w-zone-210b</metric>
    <metric>w-zone-242b</metric>
    <metric>w-zone-259b</metric>
    <metric>w-zone-260b</metric>
    <metric>w-zone-261b</metric>
    <metric>w-zone-262b</metric>
    <metric>w-zone-263b</metric>
    <metric>w-zone-264b</metric>
    <metric>w-zone-265b</metric>
    <metric>w-zone-266b</metric>
    <metric>w-zone-267b</metric>
    <metric>w-zone-268b</metric>
    <metric>w-zone-269b</metric>
    <metric>w-zone-270b</metric>
    <metric>w-zone-271b</metric>
    <metric>w-zone-272b</metric>
    <metric>w-zone-273b</metric>
    <metric>w-zone-274b</metric>
    <metric>w-zone-277b</metric>
    <metric>w-zone-278b</metric>
    <metric>w-zone-279b</metric>
    <metric>w-zone-280b</metric>
    <metric>w-zone-281b</metric>
    <metric>w-zone-282b</metric>
    <metric>w-zone-283b</metric>
    <metric>w-zone-284b</metric>
    <metric>w-zone-285b</metric>
    <metric>w-zone-286b</metric>
    <metric>w-zone-287b</metric>
    <metric>w-zone-289b</metric>
    <metric>w-zone-290b</metric>
    <metric>w-zone-291b</metric>
    <metric>r-zone-34b</metric>
    <metric>r-zone-35b</metric>
    <metric>r-zone-67b</metric>
    <metric>r-zone-69b</metric>
    <metric>r-zone-71b</metric>
    <metric>r-zone-73b</metric>
    <metric>r-zone-75b</metric>
    <metric>r-zone-77b</metric>
    <metric>r-zone-79b</metric>
    <metric>r-zone-81b</metric>
    <metric>r-zone-83b</metric>
    <metric>r-zone-85b</metric>
    <metric>r-zone-86b</metric>
    <metric>r-zone-87b</metric>
    <metric>r-zone-88b</metric>
    <metric>r-zone-89b</metric>
    <metric>r-zone-90b</metric>
    <metric>r-zone-91b</metric>
    <metric>r-zone-92b</metric>
    <metric>r-zone-93b</metric>
    <metric>r-zone-94b</metric>
    <metric>r-zone-95b</metric>
    <metric>r-zone-96b</metric>
    <metric>r-zone-97b</metric>
    <metric>r-zone-98b</metric>
    <metric>r-zone-99b</metric>
    <metric>r-zone-106b</metric>
    <metric>r-zone-107b</metric>
    <metric>r-zone-108b</metric>
    <metric>r-zone-109b</metric>
    <metric>r-zone-110b</metric>
    <metric>r-zone-111b</metric>
    <metric>r-zone-112b</metric>
    <metric>r-zone-113b</metric>
    <metric>r-zone-118b</metric>
    <metric>r-zone-119b</metric>
    <metric>r-zone-120b</metric>
    <metric>r-zone-121b</metric>
    <metric>r-zone-122b</metric>
    <metric>r-zone-123b</metric>
    <metric>r-zone-125b</metric>
    <metric>r-zone-127b</metric>
    <metric>r-zone-129b</metric>
    <metric>r-zone-139b</metric>
    <metric>r-zone-155b</metric>
    <metric>r-zone-157b</metric>
    <metric>r-zone-159b</metric>
    <metric>r-zone-161b</metric>
    <metric>r-zone-163b</metric>
    <metric>r-zone-165b</metric>
    <metric>r-zone-167b</metric>
    <metric>r-zone-180b</metric>
    <metric>r-zone-194b</metric>
    <metric>r-zone-196b</metric>
    <metric>r-zone-210b</metric>
    <metric>r-zone-242b</metric>
    <metric>r-zone-259b</metric>
    <metric>r-zone-260b</metric>
    <metric>r-zone-261b</metric>
    <metric>r-zone-262b</metric>
    <metric>r-zone-263b</metric>
    <metric>r-zone-264b</metric>
    <metric>r-zone-265b</metric>
    <metric>r-zone-266b</metric>
    <metric>r-zone-267b</metric>
    <metric>r-zone-268b</metric>
    <metric>r-zone-269b</metric>
    <metric>r-zone-270b</metric>
    <metric>r-zone-271b</metric>
    <metric>r-zone-272b</metric>
    <metric>r-zone-273b</metric>
    <metric>r-zone-274b</metric>
    <metric>r-zone-277b</metric>
    <metric>r-zone-278b</metric>
    <metric>r-zone-279b</metric>
    <metric>r-zone-280b</metric>
    <metric>r-zone-281b</metric>
    <metric>r-zone-282b</metric>
    <metric>r-zone-283b</metric>
    <metric>r-zone-284b</metric>
    <metric>r-zone-285b</metric>
    <metric>r-zone-286b</metric>
    <metric>r-zone-287b</metric>
    <metric>r-zone-289b</metric>
    <metric>r-zone-290b</metric>
    <metric>r-zone-291b</metric>
    <metric>m-patch-1b</metric>
    <metric>m-patch-2b</metric>
    <metric>m-patch-3b</metric>
    <metric>m-patch-4b</metric>
    <metric>m-patch-5b</metric>
    <metric>m-patch-6b</metric>
    <metric>m-patch-7b</metric>
    <metric>m-patch-8b</metric>
    <metric>m-patch-9b</metric>
    <metric>m-patch-10b</metric>
    <metric>m-patch-11b</metric>
    <metric>m-patch-12b</metric>
    <metric>m-patch-13b</metric>
    <metric>m-patch-14b</metric>
    <metric>t-patch-1b</metric>
    <metric>t-patch-2b</metric>
    <metric>t-patch-3b</metric>
    <metric>t-patch-4b</metric>
    <metric>t-patch-5b</metric>
    <metric>t-patch-6b</metric>
    <metric>t-patch-7b</metric>
    <metric>t-patch-8b</metric>
    <metric>t-patch-9b</metric>
    <metric>t-patch-10b</metric>
    <metric>t-patch-11b</metric>
    <metric>t-patch-12b</metric>
    <metric>t-patch-13b</metric>
    <metric>t-patch-14b</metric>
    <metric>w-patch-1b</metric>
    <metric>w-patch-2b</metric>
    <metric>w-patch-3b</metric>
    <metric>w-patch-4b</metric>
    <metric>w-patch-5b</metric>
    <metric>w-patch-6b</metric>
    <metric>w-patch-7b</metric>
    <metric>w-patch-8b</metric>
    <metric>w-patch-9b</metric>
    <metric>w-patch-10b</metric>
    <metric>w-patch-11b</metric>
    <metric>w-patch-12b</metric>
    <metric>w-patch-13b</metric>
    <metric>w-patch-14b</metric>
    <metric>r-patch-1b</metric>
    <metric>r-patch-2b</metric>
    <metric>r-patch-3b</metric>
    <metric>r-patch-4b</metric>
    <metric>r-patch-5b</metric>
    <metric>r-patch-6b</metric>
    <metric>r-patch-7b</metric>
    <metric>r-patch-8b</metric>
    <metric>r-patch-9b</metric>
    <metric>r-patch-10b</metric>
    <metric>r-patch-11b</metric>
    <metric>r-patch-12b</metric>
    <metric>r-patch-13b</metric>
    <metric>r-patch-14b</metric>
    <enumeratedValueSet variable="random-seed">
      <value value="4519"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sensitivity_new_reporters" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>p11</metric>
    <metric>p12</metric>
    <metric>p13</metric>
    <metric>p14</metric>
    <metric>p21</metric>
    <metric>p22</metric>
    <metric>p23</metric>
    <metric>p24</metric>
    <metric>p31</metric>
    <metric>p32</metric>
    <metric>p33</metric>
    <metric>p34</metric>
    <metric>p41</metric>
    <metric>p42</metric>
    <metric>p43</metric>
    <metric>p44</metric>
    <metric>tc11</metric>
    <metric>tc12</metric>
    <metric>tc13</metric>
    <metric>tc14</metric>
    <metric>tc21</metric>
    <metric>tc22</metric>
    <metric>tc23</metric>
    <metric>tc24</metric>
    <metric>tc31</metric>
    <metric>tc32</metric>
    <metric>tc33</metric>
    <metric>tc34</metric>
    <metric>tc41</metric>
    <metric>tc42</metric>
    <metric>tc43</metric>
    <metric>tc44</metric>
    <metric>mu-max</metric>
    <metric>p-patch-spread</metric>
    <metric>p-water-spread</metric>
    <metric>p-food-dropped</metric>
    <metric>p-condensation</metric>
    <metric>p-random-noise</metric>
    <metric>p-zone4-intro</metric>
    <metric>spread-tc</metric>
    <metric>m-crate-prevalence</metric>
    <metric>t-crate-prevalence</metric>
    <metric>w-crate-prevalence</metric>
    <metric>r-crate-prevalence</metric>
    <metric>f-crate-prevalence</metric>
    <metric>fruit-conc</metric>
    <metric>fruit-tc</metric>
    <metric>zone4-load</metric>
    <metric>random-load</metric>
    <metric>m-san-reduction</metric>
    <metric>t-san-reduction</metric>
    <metric>w-san-reduction</metric>
    <metric>r-san-reduction</metric>
    <metric>f-san-reduction</metric>
    <metric>zones-initial</metric>
    <metric>z1-initial</metric>
    <metric>z2-initial</metric>
    <metric>z3-initial</metric>
    <metric>p-initial</metric>
    <metric>e-initial</metric>
    <metric>m-zones-pre</metric>
    <metric>t-zones-pre</metric>
    <metric>w-zones-pre</metric>
    <metric>r-zones-pre</metric>
    <metric>f-zones-pre</metric>
    <metric>m-zones-beg</metric>
    <metric>t-zones-beg</metric>
    <metric>w-zones-beg</metric>
    <metric>r-zones-beg</metric>
    <metric>f-zones-beg</metric>
    <metric>m-zones-mid</metric>
    <metric>t-zones-mid</metric>
    <metric>w-zones-mid</metric>
    <metric>r-zones-mid</metric>
    <metric>f-zones-mid</metric>
    <metric>m-zones-post</metric>
    <metric>t-zones-post</metric>
    <metric>w-zones-post</metric>
    <metric>r-zones-post</metric>
    <metric>f-zones-post</metric>
    <metric>zone-to-zone</metric>
    <metric>avg-listeria-transferred</metric>
    <metric>zone-to-patch</metric>
    <metric>condensation-drip</metric>
    <metric>introduction-zone4</metric>
    <metric>patch-to-patch</metric>
    <metric>introduction-food</metric>
    <metric>sum-fruit-load</metric>
    <metric>sum-fruit-transfer</metric>
    <metric>random-event</metric>
    <metric>employee-FCS</metric>
    <metric>NFCS-FCS</metric>
    <metric>total-contam-events</metric>
    <metric>m-random-event</metric>
    <metric>t-random-event</metric>
    <metric>w-random-event</metric>
    <metric>r-random-event</metric>
    <metric>f-random-event</metric>
    <metric>median-z1-time</metric>
    <metric>median-z2-time</metric>
    <metric>median-z3-time</metric>
    <metric>median-e-time</metric>
    <metric>median-p-time</metric>
    <metric>median-fwj-time</metric>
    <metric>median-z1-max-contam-bout</metric>
    <metric>median-z2-max-contam-bout</metric>
    <metric>median-z3-max-contam-bout</metric>
    <metric>median-e-max-contam-bout</metric>
    <metric>median-p-max-contam-bout</metric>
    <metric>median-fwj-max-contam-bout</metric>
    <metric>temporary-niches</metric>
    <metric>p-temp-niches</metric>
    <metric>z1-temp-niches</metric>
    <metric>z2-temp-niches</metric>
    <metric>z3-temp-niches</metric>
    <metric>z1-contacts</metric>
    <metric>z1-transfers</metric>
    <metric>z2-contacts</metric>
    <metric>z2-transfers</metric>
    <metric>z3-contacts</metric>
    <metric>z3-transfers</metric>
    <metric>z1-prev</metric>
    <metric>z2-prev</metric>
    <metric>z3-prev</metric>
    <metric>concentration-list</metric>
    <metric>time-contaminated-list</metric>
    <metric>max-contam-bout-list</metric>
    <metric>contacts-list</metric>
    <metric>transfers-list</metric>
    <metric>temp-niche-list</metric>
    <metric>niches-estab-list</metric>
    <metric>fruit-cont-list</metric>
    <metric>random-cont-list</metric>
    <metric>random-floor</metric>
    <metric>random-ceiling</metric>
    <metric>random-agent</metric>
    <metric>K</metric>
    <metric>prob-cleanable</metric>
    <metric>belt-prev</metric>
    <metric>cart-prev</metric>
    <metric>control-panel-prev</metric>
    <metric>door-prev</metric>
    <metric>drain-prev</metric>
    <metric>employee-prev</metric>
    <metric>flume-prev</metric>
    <metric>frame-prev</metric>
    <metric>frame-in-prev</metric>
    <metric>ladder-prev</metric>
    <metric>misc-prev</metric>
    <metric>packing-prev</metric>
    <metric>squeegee-prev</metric>
    <metric>trash-gray-prev</metric>
    <metric>trash-white-prev</metric>
    <metric>trash-yellow-prev</metric>
    <metric>wall-prev</metric>
    <metric>weigher-prev</metric>
    <metric>a-b-prev</metric>
    <metric>a-flume-prev</metric>
    <metric>a-line-prev</metric>
    <metric>a-under-room-prev</metric>
    <metric>b-flume-prev</metric>
    <metric>b-line-prev</metric>
    <metric>coring-prev</metric>
    <metric>elevator-prev</metric>
    <metric>hopper-prev</metric>
    <metric>packing--prev</metric>
    <metric>scoring-prev</metric>
    <metric>sorting-prev</metric>
    <metric>spin-dryer-prev</metric>
    <metric>weigher-platform-prev</metric>
    <metric>A-prev</metric>
    <metric>B-prev</metric>
    <metric>floors-prev</metric>
    <metric>belt-samples</metric>
    <metric>cart-samples</metric>
    <metric>control-panel-samples</metric>
    <metric>door-samples</metric>
    <metric>drain-samples</metric>
    <metric>employee-samples</metric>
    <metric>floor-samples</metric>
    <metric>flume-samples</metric>
    <metric>frame-samples</metric>
    <metric>frame-in-samples</metric>
    <metric>ladder-samples</metric>
    <metric>misc-samples</metric>
    <metric>packing-samples</metric>
    <metric>squeegee-samples</metric>
    <metric>trash-gray-samples</metric>
    <metric>trash-white-samples</metric>
    <metric>trash-yellow-samples</metric>
    <metric>wall-samples</metric>
    <metric>weigher-samples</metric>
    <metric>a-b-samples</metric>
    <metric>a-flume-samples</metric>
    <metric>a-line-samples</metric>
    <metric>a-under-room-samples</metric>
    <metric>b-flume-samples</metric>
    <metric>b-line-samples</metric>
    <metric>coring-samples</metric>
    <metric>elevator-samples</metric>
    <metric>hopper-samples</metric>
    <metric>packing--samples</metric>
    <metric>scoring-samples</metric>
    <metric>sorting-samples</metric>
    <metric>spin-dryer-samples</metric>
    <metric>weigher-platform-samples</metric>
    <metric>A-samples</metric>
    <metric>B-samples</metric>
    <enumeratedValueSet variable="random-seed">
      <value value="4519"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sensitivity_prevalence" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>p11</metric>
    <metric>p12</metric>
    <metric>p13</metric>
    <metric>p14</metric>
    <metric>p21</metric>
    <metric>p22</metric>
    <metric>p23</metric>
    <metric>p24</metric>
    <metric>p31</metric>
    <metric>p32</metric>
    <metric>p33</metric>
    <metric>p34</metric>
    <metric>p41</metric>
    <metric>p42</metric>
    <metric>p43</metric>
    <metric>p44</metric>
    <metric>tc11</metric>
    <metric>tc12</metric>
    <metric>tc13</metric>
    <metric>tc14</metric>
    <metric>tc21</metric>
    <metric>tc22</metric>
    <metric>tc23</metric>
    <metric>tc24</metric>
    <metric>tc31</metric>
    <metric>tc32</metric>
    <metric>tc33</metric>
    <metric>tc34</metric>
    <metric>tc41</metric>
    <metric>tc42</metric>
    <metric>tc43</metric>
    <metric>tc44</metric>
    <metric>mu-max</metric>
    <metric>p-patch-spread</metric>
    <metric>p-water-spread</metric>
    <metric>p-food-dropped</metric>
    <metric>p-condensation</metric>
    <metric>p-random-noise</metric>
    <metric>p-zone4-intro</metric>
    <metric>spread-tc</metric>
    <metric>m-crate-prevalence</metric>
    <metric>t-crate-prevalence</metric>
    <metric>w-crate-prevalence</metric>
    <metric>r-crate-prevalence</metric>
    <metric>f-crate-prevalence</metric>
    <metric>fruit-conc</metric>
    <metric>fruit-tc</metric>
    <metric>zone4-load</metric>
    <metric>random-load</metric>
    <metric>m-san-reduction</metric>
    <metric>t-san-reduction</metric>
    <metric>w-san-reduction</metric>
    <metric>r-san-reduction</metric>
    <metric>f-san-reduction</metric>
    <metric>zones-initial</metric>
    <metric>z1-initial</metric>
    <metric>z2-initial</metric>
    <metric>z3-initial</metric>
    <metric>p-initial</metric>
    <metric>e-initial</metric>
    <metric>m-zones-pre</metric>
    <metric>t-zones-pre</metric>
    <metric>w-zones-pre</metric>
    <metric>r-zones-pre</metric>
    <metric>f-zones-pre</metric>
    <metric>m-zones-beg</metric>
    <metric>t-zones-beg</metric>
    <metric>w-zones-beg</metric>
    <metric>r-zones-beg</metric>
    <metric>f-zones-beg</metric>
    <metric>m-zones-mid</metric>
    <metric>t-zones-mid</metric>
    <metric>w-zones-mid</metric>
    <metric>r-zones-mid</metric>
    <metric>f-zones-mid</metric>
    <metric>m-zones-post</metric>
    <metric>t-zones-post</metric>
    <metric>w-zones-post</metric>
    <metric>r-zones-post</metric>
    <metric>f-zones-post</metric>
    <metric>zone-to-zone</metric>
    <metric>avg-listeria-transferred</metric>
    <metric>zone-to-patch</metric>
    <metric>condensation-drip</metric>
    <metric>introduction-zone4</metric>
    <metric>patch-to-patch</metric>
    <metric>introduction-food</metric>
    <metric>sum-fruit-load</metric>
    <metric>sum-fruit-transfer</metric>
    <metric>random-event</metric>
    <metric>employee-FCS</metric>
    <metric>NFCS-FCS</metric>
    <metric>total-contam-events</metric>
    <metric>m-random-event</metric>
    <metric>t-random-event</metric>
    <metric>w-random-event</metric>
    <metric>r-random-event</metric>
    <metric>f-random-event</metric>
    <metric>median-z1-time</metric>
    <metric>median-z2-time</metric>
    <metric>median-z3-time</metric>
    <metric>median-e-time</metric>
    <metric>median-p-time</metric>
    <metric>median-fwj-time</metric>
    <metric>median-z1-max-contam-bout</metric>
    <metric>median-z2-max-contam-bout</metric>
    <metric>median-z3-max-contam-bout</metric>
    <metric>median-e-max-contam-bout</metric>
    <metric>median-p-max-contam-bout</metric>
    <metric>median-fwj-max-contam-bout</metric>
    <metric>temporary-niches</metric>
    <metric>p-temp-niches</metric>
    <metric>z1-temp-niches</metric>
    <metric>z2-temp-niches</metric>
    <metric>z3-temp-niches</metric>
    <metric>z1-contacts</metric>
    <metric>z1-transfers</metric>
    <metric>z2-contacts</metric>
    <metric>z2-transfers</metric>
    <metric>z3-contacts</metric>
    <metric>z3-transfers</metric>
    <metric>z1-prev</metric>
    <metric>z2-prev</metric>
    <metric>z3-prev</metric>
    <metric>concentration-list</metric>
    <metric>time-contaminated-list</metric>
    <metric>max-contam-bout-list</metric>
    <metric>contacts-list</metric>
    <metric>transfers-list</metric>
    <metric>temp-niche-list</metric>
    <metric>niches-estab-list</metric>
    <metric>fruit-cont-list</metric>
    <metric>random-cont-list</metric>
    <metric>random-floor</metric>
    <metric>random-ceiling</metric>
    <metric>random-agent</metric>
    <metric>K</metric>
    <metric>prob-cleanable</metric>
    <metric>belt-prev</metric>
    <metric>cart-prev</metric>
    <metric>control-panel-prev</metric>
    <metric>door-prev</metric>
    <metric>drain-prev</metric>
    <metric>employee-prev</metric>
    <metric>flume-prev</metric>
    <metric>frame-prev</metric>
    <metric>frame-in-prev</metric>
    <metric>ladder-prev</metric>
    <metric>misc-prev</metric>
    <metric>packing-prev</metric>
    <metric>squeegee-prev</metric>
    <metric>trash-gray-prev</metric>
    <metric>trash-white-prev</metric>
    <metric>trash-yellow-prev</metric>
    <metric>wall-prev</metric>
    <metric>weigher-prev</metric>
    <metric>a-b-prev</metric>
    <metric>a-flume-prev</metric>
    <metric>a-line-prev</metric>
    <metric>a-under-room-prev</metric>
    <metric>b-flume-prev</metric>
    <metric>b-line-prev</metric>
    <metric>coring-prev</metric>
    <metric>elevator-prev</metric>
    <metric>hopper-prev</metric>
    <metric>packing--prev</metric>
    <metric>scoring-prev</metric>
    <metric>sorting-prev</metric>
    <metric>spin-dryer-prev</metric>
    <metric>weigher-platform-prev</metric>
    <metric>A-prev</metric>
    <metric>B-prev</metric>
    <metric>floors-prev</metric>
    <metric>belt-samples</metric>
    <metric>cart-samples</metric>
    <metric>control-panel-samples</metric>
    <metric>door-samples</metric>
    <metric>drain-samples</metric>
    <metric>employee-samples</metric>
    <metric>floor-samples</metric>
    <metric>flume-samples</metric>
    <metric>frame-samples</metric>
    <metric>frame-in-samples</metric>
    <metric>ladder-samples</metric>
    <metric>misc-samples</metric>
    <metric>packing-samples</metric>
    <metric>squeegee-samples</metric>
    <metric>trash-gray-samples</metric>
    <metric>trash-white-samples</metric>
    <metric>trash-yellow-samples</metric>
    <metric>wall-samples</metric>
    <metric>weigher-samples</metric>
    <metric>a-b-samples</metric>
    <metric>a-flume-samples</metric>
    <metric>a-line-samples</metric>
    <metric>a-under-room-samples</metric>
    <metric>b-flume-samples</metric>
    <metric>b-line-samples</metric>
    <metric>coring-samples</metric>
    <metric>elevator-samples</metric>
    <metric>hopper-samples</metric>
    <metric>packing--samples</metric>
    <metric>scoring-samples</metric>
    <metric>sorting-samples</metric>
    <metric>spin-dryer-samples</metric>
    <metric>weigher-platform-samples</metric>
    <metric>A-samples</metric>
    <metric>B-samples</metric>
    <enumeratedValueSet variable="random-seed">
      <value value="4519"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="scenario-trials-claire" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>p11</metric>
    <metric>p12</metric>
    <metric>p13</metric>
    <metric>p14</metric>
    <metric>p21</metric>
    <metric>p22</metric>
    <metric>p23</metric>
    <metric>p24</metric>
    <metric>p31</metric>
    <metric>p32</metric>
    <metric>p33</metric>
    <metric>p34</metric>
    <metric>p41</metric>
    <metric>p42</metric>
    <metric>p43</metric>
    <metric>p44</metric>
    <metric>tc11</metric>
    <metric>tc12</metric>
    <metric>tc13</metric>
    <metric>tc14</metric>
    <metric>tc21</metric>
    <metric>tc22</metric>
    <metric>tc23</metric>
    <metric>tc24</metric>
    <metric>tc31</metric>
    <metric>tc32</metric>
    <metric>tc33</metric>
    <metric>tc34</metric>
    <metric>tc41</metric>
    <metric>tc42</metric>
    <metric>tc43</metric>
    <metric>tc44</metric>
    <metric>mu-max</metric>
    <metric>p-patch-spread</metric>
    <metric>p-water-spread</metric>
    <metric>p-food-dropped</metric>
    <metric>p-condensation</metric>
    <metric>p-random-noise</metric>
    <metric>p-zone4-intro</metric>
    <metric>spread-tc</metric>
    <metric>m-crate-prevalence</metric>
    <metric>t-crate-prevalence</metric>
    <metric>w-crate-prevalence</metric>
    <metric>r-crate-prevalence</metric>
    <metric>f-crate-prevalence</metric>
    <metric>fruit-conc</metric>
    <metric>fruit-tc</metric>
    <metric>zone4-load</metric>
    <metric>random-load</metric>
    <metric>m-san-reduction</metric>
    <metric>t-san-reduction</metric>
    <metric>w-san-reduction</metric>
    <metric>r-san-reduction</metric>
    <metric>f-san-reduction</metric>
    <metric>zones-initial</metric>
    <metric>z1-initial</metric>
    <metric>z2-initial</metric>
    <metric>z3-initial</metric>
    <metric>p-initial</metric>
    <metric>e-initial</metric>
    <metric>m-zones-pre</metric>
    <metric>t-zones-pre</metric>
    <metric>w-zones-pre</metric>
    <metric>r-zones-pre</metric>
    <metric>f-zones-pre</metric>
    <metric>m-zones-beg</metric>
    <metric>t-zones-beg</metric>
    <metric>w-zones-beg</metric>
    <metric>r-zones-beg</metric>
    <metric>f-zones-beg</metric>
    <metric>m-zones-mid</metric>
    <metric>t-zones-mid</metric>
    <metric>w-zones-mid</metric>
    <metric>r-zones-mid</metric>
    <metric>f-zones-mid</metric>
    <metric>m-zones-post</metric>
    <metric>t-zones-post</metric>
    <metric>w-zones-post</metric>
    <metric>r-zones-post</metric>
    <metric>f-zones-post</metric>
    <metric>zone-to-zone</metric>
    <metric>avg-listeria-transferred</metric>
    <metric>zone-to-patch</metric>
    <metric>condensation-drip</metric>
    <metric>introduction-zone4</metric>
    <metric>patch-to-patch</metric>
    <metric>introduction-food</metric>
    <metric>sum-fruit-load</metric>
    <metric>sum-fruit-transfer</metric>
    <metric>random-event</metric>
    <metric>employee-FCS</metric>
    <metric>NFCS-FCS</metric>
    <metric>total-contam-events</metric>
    <metric>m-random-event</metric>
    <metric>t-random-event</metric>
    <metric>w-random-event</metric>
    <metric>r-random-event</metric>
    <metric>f-random-event</metric>
    <metric>median-z1-time</metric>
    <metric>median-z2-time</metric>
    <metric>median-z3-time</metric>
    <metric>median-e-time</metric>
    <metric>median-p-time</metric>
    <metric>median-fwj-time</metric>
    <metric>median-z1-max-contam-bout</metric>
    <metric>median-z2-max-contam-bout</metric>
    <metric>median-z3-max-contam-bout</metric>
    <metric>median-e-max-contam-bout</metric>
    <metric>median-p-max-contam-bout</metric>
    <metric>median-fwj-max-contam-bout</metric>
    <metric>temporary-niches</metric>
    <metric>p-temp-niches</metric>
    <metric>z1-temp-niches</metric>
    <metric>z2-temp-niches</metric>
    <metric>z3-temp-niches</metric>
    <metric>z1-contacts</metric>
    <metric>z1-transfers</metric>
    <metric>z2-contacts</metric>
    <metric>z2-transfers</metric>
    <metric>z3-contacts</metric>
    <metric>z3-transfers</metric>
    <metric>z1-prev</metric>
    <metric>z2-prev</metric>
    <metric>z3-prev</metric>
    <metric>belt-prev</metric>
    <metric>cart-prev</metric>
    <metric>control-panel-prev</metric>
    <metric>door-prev</metric>
    <metric>drain-prev</metric>
    <metric>employee-prev</metric>
    <metric>flume-prev</metric>
    <metric>frame-prev</metric>
    <metric>frame-in-prev</metric>
    <metric>ladder-prev</metric>
    <metric>misc-prev</metric>
    <metric>packing-prev</metric>
    <metric>squeegee-prev</metric>
    <metric>trash-gray-prev</metric>
    <metric>trash-white-prev</metric>
    <metric>trash-yellow-prev</metric>
    <metric>wall-prev</metric>
    <metric>weigher-prev</metric>
    <metric>a-b-prev</metric>
    <metric>a-flume-prev</metric>
    <metric>a-line-prev</metric>
    <metric>a-under-room-prev</metric>
    <metric>b-flume-prev</metric>
    <metric>b-line-prev</metric>
    <metric>coring-prev</metric>
    <metric>elevator-prev</metric>
    <metric>hopper-prev</metric>
    <metric>packing--prev</metric>
    <metric>scoring-prev</metric>
    <metric>sorting-prev</metric>
    <metric>spin-dryer-prev</metric>
    <metric>weigher-platform-prev</metric>
    <metric>A-prev</metric>
    <metric>B-prev</metric>
    <metric>floors-prev</metric>
    <metric>belt-samples</metric>
    <metric>cart-samples</metric>
    <metric>control-panel-samples</metric>
    <metric>door-samples</metric>
    <metric>drain-samples</metric>
    <metric>employee-samples</metric>
    <metric>floor-samples</metric>
    <metric>flume-samples</metric>
    <metric>frame-samples</metric>
    <metric>frame-in-samples</metric>
    <metric>ladder-samples</metric>
    <metric>misc-samples</metric>
    <metric>packing-samples</metric>
    <metric>squeegee-samples</metric>
    <metric>trash-gray-samples</metric>
    <metric>trash-white-samples</metric>
    <metric>trash-yellow-samples</metric>
    <metric>wall-samples</metric>
    <metric>weigher-samples</metric>
    <metric>a-b-samples</metric>
    <metric>a-flume-samples</metric>
    <metric>a-line-samples</metric>
    <metric>a-under-room-samples</metric>
    <metric>b-flume-samples</metric>
    <metric>b-line-samples</metric>
    <metric>coring-samples</metric>
    <metric>elevator-samples</metric>
    <metric>hopper-samples</metric>
    <metric>packing--samples</metric>
    <metric>scoring-samples</metric>
    <metric>sorting-samples</metric>
    <metric>spin-dryer-samples</metric>
    <metric>weigher-platform-samples</metric>
    <metric>A-samples</metric>
    <metric>B-samples</metric>
    <metric>concentration-list</metric>
    <metric>time-contaminated-list</metric>
    <metric>max-contam-bout-list</metric>
    <metric>contacts-list</metric>
    <metric>transfers-list</metric>
    <metric>temp-niche-list</metric>
    <metric>niches-estab-list</metric>
    <metric>fruit-cont-list</metric>
    <metric>random-cont-list</metric>
    <metric>random-floor</metric>
    <metric>random-ceiling</metric>
    <metric>random-agent</metric>
    <metric>a-table</metric>
    <metric>b-table</metric>
    <metric>c-table</metric>
    <metric>d-table</metric>
    <metric>a-table2</metric>
    <metric>b-table2</metric>
    <metric>c-table2</metric>
    <metric>d-table2</metric>
    <enumeratedValueSet variable="random-seed">
      <value value="4519"/>
    </enumeratedValueSet>
    <steppedValueSet variable="seed-me" first="1" step="1" last="1000"/>
    <enumeratedValueSet variable="sampling-time">
      <value value="&quot;first-hour&quot;"/>
      <value value="&quot;mid-production&quot;"/>
      <value value="&quot;second-shift&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sampling-sites">
      <value value="&quot;Scenario-1&quot;"/>
      <value value="&quot;Scenario-2&quot;"/>
      <value value="&quot;Scenario-3&quot;"/>
      <value value="&quot;Scenario-4&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
