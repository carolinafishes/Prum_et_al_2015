##begin by reading in internal functions from PI_functions file###########################
###########################
###########################
###########################
source('~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/PI_functions.R', chdir = TRUE)


##read in the tree
read.tree("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/ALLTREE-lambda-0.1.phy")->tree

##for PI visualization let's drop the outgroup crocodiles
c("I0460_HERP_14722_Crocodylidae_Crocodylus_porosus","I0461_HERP_15451_Alligatoridae_Caiman_croccodilus")->dropping
drop.tip(tree,dropping)->tree2
####get the rate vectors, i was just pasting these in
##set working directory
setwd("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness")


###Next read in all the rate vector files 
source("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/loci1_30Allrates_file.allrates")
source("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/loci31_62Allrates_file.allrates")
source("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/loci63_100Allrates_file.allrates")
source("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/loci101_141Allrates_file.allrates")
source("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/loci142_197Allrates_file.allrates")
source("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/loci200_240Allrates_file.allrates")
source("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/242_280Allrates_file.allrates")
source("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/281_314Allrates_file.allrates")
source("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/314_349Allrates_file.allrates")

##note on the above, we encountered memory limitations getting all the individual rate vectors, so we broke the dataset into stages to save memory. With a better machine you could just estimates all the site rates at once. 


##here we are generating PI profiles
L1<-ratesA[1:1594]
L2<-ratesA[1595:3354]
L3<-ratesA[3355:4990]
L4<-ratesA[4991:6483]
L5<-ratesA[6484:8117]
L6<-ratesA[8118:9799]
L7<-ratesA[9800:10831]
L8<-ratesA[10832:12320]
L9<-ratesA[12321:14008]

ratesApart1<-ratesA[1:14008]
as.matrix(ratesApart1)-> ratesAApart1
upperA_1<-c(1594,3354, 4990, 6483, 8117, 9799, 10831, 12320, 14008)
lowerA_1<-c(0, 1595, 3355, 4991, 6484, 8118, 9800, 10832, 12321 )
cbind(lowerA_1, upperA_1)->breaksA1
defined.multi.profile(ratesAApart1,tree, breaksA1)->ll
write.table(ll, file="profilebylocus", quote=FALSE)


L10<-ratesA[14009:15201]
L11<-ratesA[15202:16983]
L12<-ratesA[16984:18640]
L21<-ratesA[18641:19530]
L22<-ratesA[19531:21236]
L24<-ratesA[21237:22925]
L27<-ratesA[22926:23933]
L28<-ratesA[23934:25539]
L30<-ratesA[25540:27222]
ratesApart2<-ratesA[14009:27222]
as.matrix(ratesApart2)->ratesAApart2
length(ratesAApart2)->subtract
upperA_2<-c(15201-subtract, 16983-subtract, 18640-subtract, 19530-subtract, 21236-subtract, 22925-subtract, 23933-subtract, 25539-subtract, 27222-subtract)
lowerA_2<-c(14009-subtract, 15202-subtract, 16984-subtract, 18641-subtract, 19531-subtract, 21237-subtract, 22926-subtract, 23934-subtract, 25540-subtract)
cbind(lowerA_2, upperA_2)->breaksA2
defined.multi.profile(ratesAApart2,tree, breaksA2)->ll2
write.table(ll2, file="profilebylocus", quote=FALSE, append=TRUE)


### For the purpose of illustration we split the above two subsets out.   Let's write the rest of the profiles to file, ignore the warnings telling you that you are appending. These are still broken into subsets to make tracking easier. 

L31<-ratesB[1:1743]
L32<-ratesB[1744:3155]
L34<-ratesB[3156:4772]
L35<-ratesB[4773:6490]
L36<-ratesB[6491:7839]
L37<-ratesB[7840:9461]
L38<-ratesB[9462:11134]
L39<-ratesB[11135:12773]
L41<-ratesB[12774:14451]
L42<-ratesB[14452:15980]
L43<-ratesB[15981:17635]
L44<-ratesB[17636:19274]
L48<-ratesB[19275:20892]
L49<-ratesB[20893:22501]
L50<-ratesB[22502:24067]
L51<-ratesB[24068:24676]
L53<-ratesB[24677:26336]
L54<-ratesB[26337:27837]
L57<-ratesB[27838:28703]
L58<-ratesB[28704:30437]
L59<-ratesB[30438:32254]
L61<-ratesB[32255:33594]
L62<-ratesB[33595:35320]
as.matrix(ratesB)-> ratesBm
upperA_1<-c(1743, 3155, 4772, 6490, 7839, 9461, 11134, 12773, 14451, 15980, 17635, 19274, 20892, 22501, 24067, 24676, 26336, 27837, 28703, 30437, 32254, 33594, 35320)
lowerA_1<-c(0, 1744, 3156, 4773, 6491, 7840, 9462, 11135, 12774, 14452, 15981, 17636, 19275, 20893, 22502, 24068, 24677, 26337, 27838, 28704, 30438, 32255, 33595)
cbind(lowerA_1, upperA_1)->breaksB
defined.multi.profile(ratesBm,tree, breaksB)->ll3
write.table(ll3, file="profilebylocus", quote=FALSE, append=TRUE)


L63<-ratesC[1:1639] 
L64<-ratesC[1640:3343] 
L65<-ratesC[3344:4969] 
L67<-ratesC[4970:6692] 
L72<-ratesC[6693:8126] 
L73<-ratesC[8127:9466] 
L74<-ratesC[9467:10561] 
L75<-ratesC[10562:12104] 
L76<-ratesC[12105:13564] 
L77<-ratesC[13565:14900] 
L78<-ratesC[14901:16176] 
L79<-ratesC[16177:17743] 
L80<-ratesC[17744:19041] 
L81<-ratesC[19042:20639] 
L82<-ratesC[20640:22331] 
L85<-ratesC[22332:24059] 
L86<-ratesC[24060:25809] 
L87<-ratesC[25810:27611] 
L88<-ratesC[27612:29278] 
L89<-ratesC[29279:30703] 
L91<-ratesC[30704:32444] 
L94<-ratesC[32445:33290] 
L95<-ratesC[33291:34744] 
L98<-ratesC[34745:36509] 
L99<-ratesC[36510:38221] 
L100<-ratesC[38222:39639]
as.matrix(ratesC)-> ratesCm
upperA_1<-c(1639, 3343, 4969, 6692, 8126, 9466, 10561, 12104, 13564, 14900, 16176, 17743, 19041, 20639, 22331, 24059, 25809, 27611, 29278, 30703, 32444, 33290, 34744, 36509, 38221, 39639)
lowerA_1<-c(0, 1640, 3344, 4970, 6693, 8127, 9467, 10562, 12105, 13565, 14901, 16177, 17744, 19042, 20640, 22332, 24060, 25810, 27612, 29279, 30704, 32445, 33291, 34745, 36510, 38222 )
cbind(lowerA_1, upperA_1)->breaksC
defined.multi.profile(ratesCm,tree, breaksC)->ll4
write.table(ll4, file="profilebylocus", quote=FALSE, append=TRUE)


L102<-ratesD[1:1351] 
L103<-ratesD[1352:3040] 
L104<-ratesD[3041:4772]
L105<-ratesD[4773:6235] 
L107<-ratesD[6236:7937] 
L108<-ratesD[7938:9637] 
L110<-ratesD[9638:11009] 
L111<-ratesD[11010:12347] 
L113<-ratesD[12348:14304] 
L115<-ratesD[14305:15037] 
L116<-ratesD[15038:16416] 
L118<-ratesD[16417:17430] 
L119<-ratesD[17431:19016] 
L120<-ratesD[19017:20594] 
L121<-ratesD[20595:22531] 
L122<-ratesD[22532:24477] 
L125<-ratesD[24478:26123] 
L126<-ratesD[26124:27793] 
L127<-ratesD[27794:29519] 
L128<-ratesD[29520:30838] 
L129<-ratesD[30839:32483] 
L130<-ratesD[32484:34253] 
L131<-ratesD[34254:35386] 
L133<-ratesD[35387:36951] 
L134<-ratesD[36952:38522] 
L135<-ratesD[38523:40150] 
L136<-ratesD[40151:42208] 
L137<-ratesD[42209:43745] 
L141<-ratesD[43746:44950]

as.matrix(ratesD)-> ratesDm
upperA_1<-c(1351,3040, 4772, 6235, 7937, 9637, 11009, 12347, 14304, 15037, 16416, 17430, 19016, 20594, 22531, 24477, 26123, 27793, 29519, 30838, 32483, 34253, 35386, 36951, 38522, 40150, 42208, 43745, 44950)
lowerA_1<-c(0, 1352, 3041, 4773, 6236, 7938, 9638, 11010, 12348, 14305, 15038, 16417, 17431, 19017, 20595, 22532, 24478, 26124, 27794, 29520, 30839, 32484, 34254, 35387, 36952, 38523, 40151, 42209, 43746)
cbind(lowerA_1, upperA_1)->breaksD
defined.multi.profile(ratesDm,tree, breaksD)->ll5
write.table(ll5, file="profilebylocus", quote=FALSE, append=TRUE)



L142<-ratesE[1:1728] 
L143<-ratesE[1729:3316] 
L145<-ratesE[3317:5021] 
L146<-ratesE[5022:6809] 
L147<-ratesE[6810:8445] 
L148<-ratesE[8446:9774] 
L149<-ratesE[9775:11581] 
L150<-ratesE[11582:13436] 
L151<-ratesE[13437:15147] 
L152<-ratesE[15148:16584] 
L153<-ratesE[16585:18687] 
L154<-ratesE[18688:20374] 
L155<-ratesE[20375:22047] 
L157<-ratesE[22048:23581] 
L158<-ratesE[23582:25368] 
L159<-ratesE[25369:26648] 
L161<-ratesE[26649:28612] 
L162<-ratesE[28613:30146] 
L163<-ratesE[30147:31220] 
L165<-ratesE[31221:32704] 
L169<-ratesE[32705:34487] 
L170<-ratesE[34488:36262] 
L171<-ratesE[36263:38252] 
L174<-ratesE[38253:39153] 
L175<-ratesE[39154:40882] 
L176<-ratesE[40883:42794] 
L177<-ratesE[42795:44516] 
L178<-ratesE[44517:46226] 
L180<-ratesE[46227:47912] 
L181<-ratesE[47913:49716] 
L182<-ratesE[49717:51708] 
L183<-ratesE[51709:53072] 
L184<-ratesE[53073:54660] 
L186<-ratesE[54661:56124] 
L187<-ratesE[56125:57793] 
L188<-ratesE[57794:59329] 
L189<-ratesE[59330:61222] 
L190<-ratesE[61223:62883] 
L192<-ratesE[62884:64243] 
L193<-ratesE[64244:64715] 
L194<-ratesE[64716:66490] 
L195<-ratesE[66491:68210] 
L196<-ratesE[68211:69909] 
L197<-ratesE[69910:71170]

as.matrix(ratesE)-> ratesEm
upperA_1<-c(1728,3316, 5021, 6809, 8445, 9774, 11581, 13436, 15147, 16584, 18687, 20374, 22047, 23581, 25368, 26648, 28612, 30146, 31220, 32704, 34487, 36262, 38252, 39153, 40882, 42794, 44516, 46226, 47912, 49716, 51708, 53072, 54660, 56124, 57793, 59329, 61222, 62883, 64243, 64715, 66490, 68210, 69909, 71170)
lowerA_1<-c(0,1729, 3317, 5022, 6810, 8446, 9775,11582, 13437, 15148, 16585, 18688, 20375, 22048, 23582, 25369, 26649, 28613, 30147, 31221, 32705, 34488, 36263, 38253, 39154, 40883, 42795, 44517, 46227, 47913, 49717, 51709, 53073, 54661, 56125, 57794, 59330, 61223, 62884, 64244, 64716, 66491, 68211, 69910)
cbind(lowerA_1, upperA_1)->breaksE
defined.multi.profile(ratesEm,tree, breaksE)->ll6
write.table(ll6, file="profilebylocus", quote=FALSE, append=TRUE)




L201<-ratesF[1:1060] 
L202<-ratesF[1061:2758] 
L203<-ratesF[2759:4671] 
L204<-ratesF[4672:6700] 
L205<-ratesF[6701:8387] 
L206<-ratesF[8388:9997] 
L207<-ratesF[9998:11414] 
L208<-ratesF[11415:11864] 
L209<-ratesF[11865:12949] 
L210<-ratesF[12950:14353] 
L211<-ratesF[14354:16035] 
L213<-ratesF[16036:18351] 
L214<-ratesF[18352:19979] 
L215<-ratesF[19980:21151] 
L216<-ratesF[21152:22828] 
L217<-ratesF[22829:23765] 
L219<-ratesF[23766:24489] 
L221<-ratesF[24490:26260] 
L222<-ratesF[26261:27927] 
L223<-ratesF[27928:28371] 
L224<-ratesF[28372:29642] 
L225<-ratesF[29643:30613] 
L226<-ratesF[30614:32428] 
L228<-ratesF[32429:34260] 
L229<-ratesF[34261:35864] 
L230<-ratesF[35865:37529] 
L231<-ratesF[37530:39039] 
L232<-ratesF[39040:40783] 
L233<-ratesF[40784:41629] 
L234<-ratesF[41630:43314] 
L235<-ratesF[43315:45010] 
L237<-ratesF[45011:45831] 
L239<-ratesF[45832:46680] 
L240<-ratesF[46681:48146]

as.matrix(ratesF)-> ratesFm
upperA_1<-c(1060,2758, 4671, 6700, 8387, 9997, 11414,11864, 12949, 14353, 16035, 18351, 19979, 21151, 22828, 23765, 24489, 26260, 27927, 28371, 29642, 30613, 32428, 34260, 35864, 37529, 39039, 40783, 41629, 43314, 45010, 45831, 46680, 48146)
lowerA_1<-c(0,1061, 2759, 4672, 6701, 8388, 9998, 11415, 11865, 12950, 14354, 16036, 18352, 19980, 21152, 22829, 23766, 24490, 26261, 27928, 28372, 29643, 30614, 32429, 34261, 35865, 37530, 39040, 40784, 41630, 43315, 45011, 45832, 46681)
cbind(lowerA_1, upperA_1)->breaksF
defined.multi.profile(ratesFm,tree, breaksF)->ll7
write.table(ll7, file="profilebylocus", quote=FALSE, append=TRUE)



L241<-ratesG[1:1737] 
L242<-ratesG[1738:3137] 
L243<-ratesG[3138:4101] 
L244<-ratesG[4102:5203] 
L245<-ratesG[5204:6982] 
L246<-ratesG[6983:8699] 
L247<-ratesG[8700:10305] 
L248<-ratesG[10306:11853] 
L249<-ratesG[11854:12840] 
L250<-ratesG[12841:14367] 
L251<-ratesG[14368:16098] 
L252<-ratesG[16099:17753] 
L253<-ratesG[17754:19097] 
L255<-ratesG[19098:20790] 
L256<-ratesG[20791:22514] 
L257<-ratesG[22515:24232] 
L258<-ratesG[24233:26099] 
L260<-ratesG[26100:27700] 
L261<-ratesG[27701:29447] 
L262<-ratesG[29448:31029] 
L264<-ratesG[31030:32822] 
L265<-ratesG[32823:33361] 
L266<-ratesG[33362:35547] 
L267<-ratesG[35548:37206] 
L269<-ratesG[37207:38604] 
L270<-ratesG[38605:40123] 
L271<-ratesG[40124:41740] 
L272<-ratesG[41741:42304] 
L275<-ratesG[42305:43515] 
L276<-ratesG[43516:45199] 
L277<-ratesG[45200:46347] 
L279<-ratesG[46348:47213] 
L280<-ratesG[47214:48742]

as.matrix(ratesG)-> ratesGm
upperA_1<-c(1737,3137, 4101, 5203, 6982, 8699, 10305,11853, 12840, 14367, 16098, 17753, 19097, 20790, 22514, 24232, 26099, 27700, 29447, 31029, 32822, 33361, 35547, 37206, 38604, 40123, 41740, 42304, 43515, 45199, 46347, 47213, 48742)
lowerA_1<-c(0,1738, 3138, 4102, 5204, 6983, 8700,10306, 11854, 12841, 14368, 16099, 17754, 19098, 20791, 22515, 24233, 26100, 27701, 29448, 31030, 32823, 33362, 35548, 37207, 38605, 40124, 41741, 42305, 43516, 45200, 46348, 47214)
cbind(lowerA_1, upperA_1)->breaksG
defined.multi.profile(ratesGm,tree, breaksG)->ll8
write.table(ll8, file="profilebylocus", quote=FALSE, append=TRUE)


L281<-ratesH[1:1726] 
L282<-ratesH[1727:3670] 
L288<-ratesH[3671:5341] 
L289<-ratesH[5342:7002] 
L290<-ratesH[7003:8661] 
L292<-ratesH[8662:10362] 
L293<-ratesH[10363:12103] 
L294<-ratesH[12104:13576] 
L295<-ratesH[13577:15304] 
L296<-ratesH[15305:16994] 
L297<-ratesH[16995:18550] 
L298<-ratesH[18551:20402] 
L301<-ratesH[20403:21124] 
L302<-ratesH[21125:21990] 
L303<-ratesH[21991:23590] 
L305<-ratesH[23591:24342] 
L306<-ratesH[24343:26134] 
L307<-ratesH[26135:27276] 
L310<-ratesH[27277:28956] 
L311<-ratesH[28957:30731] 
L312<-ratesH[30732:32141] 
L314<-ratesH[32142:33859]

as.matrix(ratesH)-> ratesHm
upperA_1<-c(1726,3670, 5341, 7002, 8661, 10362,12103, 13576, 15304, 16994, 18550, 20402, 21124, 21990, 23590, 24342, 26134, 27276, 28956, 30731, 32141, 33859)
lowerA_1<-c(0,1727, 3671, 5342, 7003, 8662,10363, 12104, 13577, 15305, 16995, 18551, 20403, 21125, 21991, 23591, 24343, 26135, 27277, 28957, 30732, 32142)
cbind(lowerA_1, upperA_1)->breaksH
defined.multi.profile(ratesHm,tree, breaksH)->ll9
write.table(ll9, file="profilebylocus", quote=FALSE, append=TRUE)



L315<-ratesI[1:1720] 
L316<-ratesI[1721:3516]
L317<-ratesI[3517:4945]
L318<-ratesI[4946:6375] 
L319<-ratesI[6376:7991] 
L320<-ratesI[7992:9624] 
L321<-ratesI[9625:11381] 
L322<-ratesI[11382:13087] 
L325<-ratesI[13088:14812] 
L327<-ratesI[14813:16353] 
L328<-ratesI[16354:18028] 
L329<-ratesI[18029:19779] 
L330<-ratesI[19780:21420] 
L331<-ratesI[21421:22620] 
L333<-ratesI[22621:24319] 
L334<-ratesI[24320:25487] 
L335<-ratesI[25488:27171] 
L336<-ratesI[27172:28622] 
L337<-ratesI[28623:30241] 
L338<-ratesI[30242:32164] 
L339<-ratesI[32165:33465] 
L340<-ratesI[33466:35125] 
L341<-ratesI[35126:36853] 
L342<-ratesI[36854:38590] 
L343<-ratesI[38591:39982] 
L344<-ratesI[39983:40942] 
L346<-ratesI[40943:42519] 
L347<-ratesI[42520:44188] 
L348<-ratesI[44189:45275] 
L349<-ratesI[45276:45636]

as.matrix(ratesI)-> ratesIm
upperA_1<-c(1720, 3516, 4945, 6375, 7991, 9624, 11381,13087, 14812, 16353, 18028, 19779, 21420, 22620, 24319, 25487, 27171, 28622, 30241, 32164, 33465, 35125, 36853, 38590, 39982, 40942, 42519, 44188, 45275, 45636)
lowerA_1<-c(0, 1721,3517, 4946, 6376, 7992, 9625,11382, 13088, 14813, 16354, 18029, 19780, 21421, 22621, 24320, 25488, 27172, 28623, 30242, 32165, 33466, 35126, 36854, 38591, 39983, 40943, 42520, 44189, 45276)
cbind(lowerA_1, upperA_1)->breaksI
defined.multi.profile(ratesIm,tree, breaksI)->ll10
write.table(ll10, file="profilebylocus", quote=FALSE, append=TRUE)






#######GET SIGNAL NOISE MONSTER HISTOGRAM!
##signal noise functions require matrix formatted site rates, so convert first
as.matrix(L1)->L1
as.matrix(L2)->L2
as.matrix(L3)->L3
as.matrix(L4)->L4
as.matrix(L5)->L5
as.matrix(L6)->L6
as.matrix(L7)->L7
as.matrix(L8)->L8
as.matrix(L9)->L9
as.matrix(L10)->L10
as.matrix(L11)->L11
as.matrix(L12)->L12
as.matrix(L21)->L21
as.matrix(L22)->L22
as.matrix(L24)->L24
as.matrix(L27)->L27
as.matrix(L28)->L28
as.matrix(L30)->L30
as.matrix(L31)->L31
as.matrix(L32)->L32
as.matrix(L34)->L34
as.matrix(L35)->L35
as.matrix(L36)->L36
as.matrix(L37)->L37
as.matrix(L38)->L38
as.matrix(L39)->L39
as.matrix(L41)->L41
as.matrix(L42)->L42
as.matrix(L43)->L43
as.matrix(L44)->L44
as.matrix(L48)->L48
as.matrix(L49)->L49
as.matrix(L50)->L50
as.matrix(L51)->L51
as.matrix(L53)->L53
as.matrix(L54)->L54
as.matrix(L57)->L57
as.matrix(L58)->L58
as.matrix(L59)->L59
as.matrix(L61)->L61
as.matrix(L62)->L62
as.matrix(L63)->L63
as.matrix(L64)->L64
as.matrix(L65)->L65
as.matrix(L67)->L67
as.matrix(L72)->L72
as.matrix(L73)->L73
as.matrix(L74)->L74
as.matrix(L75)->L75
as.matrix(L76)->L76
as.matrix(L77)->L77
as.matrix(L78)->L78
as.matrix(L79)->L79
as.matrix(L80)->L80
as.matrix(L81)->L81
as.matrix(L82)->L82
as.matrix(L85)->L85
as.matrix(L86)->L86
as.matrix(L87)->L87
as.matrix(L88)->L88
as.matrix(L89)->L89
as.matrix(L91)->L91
as.matrix(L94)->L94
as.matrix(L95)->L95
as.matrix(L98)->L98
as.matrix(L99)->L99
as.matrix(L100)->L100
as.matrix(L102)->L102
as.matrix(L103)->L103
as.matrix(L104)->L104
as.matrix(L105)->L105
as.matrix(L107)->L107
as.matrix(L108)->L108
as.matrix(L110)->L110
as.matrix(L111)->L111
as.matrix(L113)->L113
as.matrix(L115)->L115
as.matrix(L116)->L116
as.matrix(L118)->L118
as.matrix(L119)->L119
as.matrix(L120)->L120
as.matrix(L121)->L121
as.matrix(L122)->L122
as.matrix(L125)->L125
as.matrix(L126)->L126
as.matrix(L127)->L127
as.matrix(L128)->L128
as.matrix(L129)->L129
as.matrix(L130)->L130
as.matrix(L131)->L131
as.matrix(L133)->L133
as.matrix(L134)->L134
as.matrix(L135)->L135
as.matrix(L136)->L136
as.matrix(L137)->L137
as.matrix(L141)->L141
as.matrix(L142)->L142
as.matrix(L143)->L143
as.matrix(L145)->L145
as.matrix(L146)->L146
as.matrix(L147)->L147
as.matrix(L148)->L148
as.matrix(L149)->L149
as.matrix(L150)->L150
as.matrix(L151)->L151
as.matrix(L152)->L152
as.matrix(L153)->L153
as.matrix(L154)->L154
as.matrix(L155)->L155
as.matrix(L157)->L157
as.matrix(L158)->L158
as.matrix(L159)->L159
as.matrix(L161)->L161
as.matrix(L162)->L162
as.matrix(L163)->L163
as.matrix(L165)->L165
as.matrix(L169)->L169
as.matrix(L170)->L170
as.matrix(L171)->L171
as.matrix(L174)->L174
as.matrix(L175)->L175
as.matrix(L176)->L176
as.matrix(L177)->L177
as.matrix(L178)->L178
as.matrix(L180)->L180
as.matrix(L181)->L181
as.matrix(L182)->L182
as.matrix(L183)->L183
as.matrix(L184)->L184
as.matrix(L186)->L186
as.matrix(L187)->L187
as.matrix(L188)->L188
as.matrix(L189)->L189
as.matrix(L190)->L190
as.matrix(L192)->L192
as.matrix(L193)->L193
as.matrix(L194)->L194
as.matrix(L195)->L195
as.matrix(L196)->L196
as.matrix(L197)->L197
as.matrix(L201)->L201
as.matrix(L202)->L202
as.matrix(L203)->L203
as.matrix(L204)->L204
as.matrix(L205)->L205
as.matrix(L206)->L206
as.matrix(L207)->L207
as.matrix(L208)->L208
as.matrix(L209)->L209
as.matrix(L210)->L210
as.matrix(L211)->L211
as.matrix(L213)->L213
as.matrix(L214)->L214
as.matrix(L215)->L215
as.matrix(L216)->L216
as.matrix(L217)->L217
as.matrix(L219)->L219
as.matrix(L221)->L221
as.matrix(L222)->L222
as.matrix(L223)->L223
as.matrix(L224)->L224
as.matrix(L225)->L225
as.matrix(L226)->L226
as.matrix(L228)->L228
as.matrix(L229)->L229
as.matrix(L230)->L230
as.matrix(L231)->L231
as.matrix(L232)->L232
as.matrix(L233)->L233
as.matrix(L234)->L234
as.matrix(L235)->L235
as.matrix(L237)->L237
as.matrix(L239)->L239
as.matrix(L240)->L240
as.matrix(L241)->L241
as.matrix(L242)->L242
as.matrix(L243)->L243
as.matrix(L244)->L244
as.matrix(L245)->L245
as.matrix(L246)->L246
as.matrix(L247)->L247
as.matrix(L248)->L248
as.matrix(L249)->L249
as.matrix(L250)->L250
as.matrix(L251)->L251
as.matrix(L252)->L252
as.matrix(L253)->L253
as.matrix(L255)->L255
as.matrix(L256)->L256
as.matrix(L257)->L257
as.matrix(L258)->L258
as.matrix(L260)->L260
as.matrix(L261)->L261
as.matrix(L262)->L262
as.matrix(L264)->L264
as.matrix(L265)->L265
as.matrix(L266)->L266
as.matrix(L267)->L267
as.matrix(L269)->L269
as.matrix(L270)->L270
as.matrix(L271)->L271
as.matrix(L272)->L272
as.matrix(L275)->L275
as.matrix(L276)->L276
as.matrix(L277)->L277
as.matrix(L279)->L279
as.matrix(L280)->L280
as.matrix(L281)->L281
as.matrix(L282)->L282
as.matrix(L288)->L288
as.matrix(L289)->L289
as.matrix(L290)->L290
as.matrix(L292)->L292
as.matrix(L293)->L293
as.matrix(L294)->L294
as.matrix(L295)->L295
as.matrix(L296)->L296
as.matrix(L297)->L297
as.matrix(L298)->L298
as.matrix(L301)->L301
as.matrix(L302)->L302
as.matrix(L303)->L303
as.matrix(L305)->L305
as.matrix(L306)->L306
as.matrix(L307)->L307
as.matrix(L310)->L310
as.matrix(L311)->L311
as.matrix(L312)->L312
as.matrix(L314)->L314
as.matrix(L315)->L315
as.matrix(L316)->L316
as.matrix(L317)->L317
as.matrix(L318)->L318
as.matrix(L319)->L319
as.matrix(L320)->L320
as.matrix(L321)->L321
as.matrix(L322)->L322
as.matrix(L325)->L325
as.matrix(L327)->L327
as.matrix(L328)->L328
as.matrix(L329)->L329
as.matrix(L330)->L330
as.matrix(L331)->L331
as.matrix(L333)->L333
as.matrix(L334)->L334
as.matrix(L335)->L335
as.matrix(L336)->L336
as.matrix(L337)->L337
as.matrix(L338)->L338
as.matrix(L339)->L339
as.matrix(L340)->L340
as.matrix(L341)->L341
as.matrix(L342)->L342
as.matrix(L343)->L343
as.matrix(L344)->L344
as.matrix(L346)->L346
as.matrix(L347)->L347
as.matrix(L348)->L348
as.matrix(L349)->L349

###next calculate signal versus noise probabilities using the equations from Townsend et al. 2012

Approximator(.28,.032,L1,3)->L1pns
Approximator(.28,.032,L2,3)->L2pns
Approximator(.28,.032,L3,3)->L3pns
Approximator(.28,.032,L4,3)->L4pns
Approximator(.28,.032,L5,3)->L5pns
Approximator(.28,.032,L6,3)->L6pns
Approximator(.28,.032,L7,3)->L7pns
Approximator(.28,.032,L8,3)->L8pns
Approximator(.28,.032,L9,3)->L9pns
Approximator(.28,.032,L10,3)->L10pns
Approximator(.28,.032,L11,3)->L11pns
Approximator(.28,.032,L12,3)->L12pns
Approximator(.28,.032,L21,3)->L21pns
Approximator(.28,.032,L22,3)->L22pns
Approximator(.28,.032,L24,3)->L24pns
Approximator(.28,.032,L27,3)->L27pns
Approximator(.28,.032,L28,3)->L28pns
Approximator(.28,.032,L30,3)->L30pns
Approximator(.28,.032,L31,3)->L31pns
Approximator(.28,.032,L32,3)->L32pns
Approximator(.28,.032,L34,3)->L34pns
Approximator(.28,.032,L35,3)->L35pns
Approximator(.28,.032,L36,3)->L36pns
Approximator(.28,.032,L37,3)->L37pns
Approximator(.28,.032,L38,3)->L38pns
Approximator(.28,.032,L39,3)->L39pns
Approximator(.28,.032,L41,3)->L41pns
Approximator(.28,.032,L42,3)->L42pns
Approximator(.28,.032,L43,3)->L43pns
Approximator(.28,.032,L44,3)->L44pns
Approximator(.28,.032,L48,3)->L48pns
Approximator(.28,.032,L49,3)->L49pns
Approximator(.28,.032,L50,3)->L50pns
Approximator(.28,.032,L51,3)->L51pns
Approximator(.28,.032,L53,3)->L53pns
Approximator(.28,.032,L54,3)->L54pns
Approximator(.28,.032,L57,3)->L57pns
Approximator(.28,.032,L58,3)->L58pns
Approximator(.28,.032,L59,3)->L59pns
Approximator(.28,.032,L61,3)->L61pns
Approximator(.28,.032,L62,3)->L62pns
Approximator(.28,.032,L63,3)->L63pns
Approximator(.28,.032,L64,3)->L64pns
Approximator(.28,.032,L65,3)->L65pns
Approximator(.28,.032,L67,3)->L67pns
Approximator(.28,.032,L72,3)->L72pns
Approximator(.28,.032,L73,3)->L73pns
Approximator(.28,.032,L74,3)->L74pns
Approximator(.28,.032,L75,3)->L75pns
Approximator(.28,.032,L76,3)->L76pns
Approximator(.28,.032,L77,3)->L77pns
Approximator(.28,.032,L78,3)->L78pns
Approximator(.28,.032,L79,3)->L79pns
Approximator(.28,.032,L80,3)->L80pns
Approximator(.28,.032,L81,3)->L81pns
Approximator(.28,.032,L82,3)->L82pns
Approximator(.28,.032,L85,3)->L85pns
Approximator(.28,.032,L86,3)->L86pns
Approximator(.28,.032,L87,3)->L87pns
Approximator(.28,.032,L88,3)->L88pns
Approximator(.28,.032,L89,3)->L89pns
Approximator(.28,.032,L91,3)->L91pns
Approximator(.28,.032,L94,3)->L94pns
Approximator(.28,.032,L95,3)->L95pns
Approximator(.28,.032,L98,3)->L98pns
Approximator(.28,.032,L99,3)->L99pns
Approximator(.28,.032,L100,3)->L100pns
Approximator(.28,.032,L102,3)->L102pns
Approximator(.28,.032,L103,3)->L103pns
Approximator(.28,.032,L104,3)->L104pns
Approximator(.28,.032,L105,3)->L105pns
Approximator(.28,.032,L107,3)->L107pns
Approximator(.28,.032,L108,3)->L108pns
Approximator(.28,.032,L110,3)->L110pns
Approximator(.28,.032,L111,3)->L111pns
Approximator(.28,.032,L113,3)->L113pns
Approximator(.28,.032,L115,3)->L115pns
Approximator(.28,.032,L116,3)->L116pns
Approximator(.28,.032,L118,3)->L118pns
Approximator(.28,.032,L119,3)->L119pns
Approximator(.28,.032,L120,3)->L120pns
Approximator(.28,.032,L121,3)->L121pns
Approximator(.28,.032,L122,3)->L122pns
Approximator(.28,.032,L125,3)->L125pns
Approximator(.28,.032,L126,3)->L126pns
Approximator(.28,.032,L127,3)->L127pns
Approximator(.28,.032,L128,3)->L128pns
Approximator(.28,.032,L129,3)->L129pns
Approximator(.28,.032,L130,3)->L130pns
Approximator(.28,.032,L131,3)->L131pns
Approximator(.28,.032,L133,3)->L133pns
Approximator(.28,.032,L134,3)->L134pns
Approximator(.28,.032,L135,3)->L135pns
Approximator(.28,.032,L136,3)->L136pns
Approximator(.28,.032,L137,3)->L137pns
Approximator(.28,.032,L141,3)->L141pns
Approximator(.28,.032,L142,3)->L142pns
Approximator(.28,.032,L143,3)->L143pns
Approximator(.28,.032,L145,3)->L145pns
Approximator(.28,.032,L146,3)->L146pns
Approximator(.28,.032,L147,3)->L147pns
Approximator(.28,.032,L148,3)->L148pns
Approximator(.28,.032,L149,3)->L149pns
Approximator(.28,.032,L150,3)->L150pns
Approximator(.28,.032,L151,3)->L151pns
Approximator(.28,.032,L152,3)->L152pns
Approximator(.28,.032,L153,3)->L153pns
Approximator(.28,.032,L154,3)->L154pns
Approximator(.28,.032,L155,3)->L155pns
Approximator(.28,.032,L157,3)->L157pns
Approximator(.28,.032,L158,3)->L158pns
Approximator(.28,.032,L159,3)->L159pns
Approximator(.28,.032,L161,3)->L161pns
Approximator(.28,.032,L162,3)->L162pns
Approximator(.28,.032,L163,3)->L163pns
Approximator(.28,.032,L165,3)->L165pns
Approximator(.28,.032,L169,3)->L169pns
Approximator(.28,.032,L170,3)->L170pns
Approximator(.28,.032,L171,3)->L171pns
Approximator(.28,.032,L174,3)->L174pns
Approximator(.28,.032,L175,3)->L175pns
Approximator(.28,.032,L176,3)->L176pns
Approximator(.28,.032,L177,3)->L177pns
Approximator(.28,.032,L178,3)->L178pns
Approximator(.28,.032,L180,3)->L180pns
Approximator(.28,.032,L181,3)->L181pns
Approximator(.28,.032,L182,3)->L182pns
Approximator(.28,.032,L183,3)->L183pns
Approximator(.28,.032,L184,3)->L184pns
Approximator(.28,.032,L186,3)->L186pns
Approximator(.28,.032,L187,3)->L187pns
Approximator(.28,.032,L188,3)->L188pns
Approximator(.28,.032,L189,3)->L189pns
Approximator(.28,.032,L190,3)->L190pns
Approximator(.28,.032,L192,3)->L192pns
Approximator(.28,.032,L193,3)->L193pns
Approximator(.28,.032,L194,3)->L194pns
Approximator(.28,.032,L195,3)->L195pns
Approximator(.28,.032,L196,3)->L196pns
Approximator(.28,.032,L197,3)->L197pns
Approximator(.28,.032,L201,3)->L201pns
Approximator(.28,.032,L202,3)->L202pns
Approximator(.28,.032,L203,3)->L203pns
Approximator(.28,.032,L204,3)->L204pns
Approximator(.28,.032,L205,3)->L205pns
Approximator(.28,.032,L206,3)->L206pns
Approximator(.28,.032,L207,3)->L207pns
Approximator(.28,.032,L208,3)->L208pns
Approximator(.28,.032,L209,3)->L209pns
Approximator(.28,.032,L210,3)->L210pns
Approximator(.28,.032,L211,3)->L211pns
Approximator(.28,.032,L213,3)->L213pns
Approximator(.28,.032,L214,3)->L214pns
Approximator(.28,.032,L215,3)->L215pns
Approximator(.28,.032,L216,3)->L216pns
Approximator(.28,.032,L217,3)->L217pns
Approximator(.28,.032,L219,3)->L219pns
Approximator(.28,.032,L221,3)->L221pns
Approximator(.28,.032,L222,3)->L222pns
Approximator(.28,.032,L223,3)->L223pns
Approximator(.28,.032,L224,3)->L224pns
Approximator(.28,.032,L225,3)->L225pns
Approximator(.28,.032,L226,3)->L226pns
Approximator(.28,.032,L228,3)->L228pns
Approximator(.28,.032,L229,3)->L229pns
Approximator(.28,.032,L230,3)->L230pns
Approximator(.28,.032,L231,3)->L231pns
Approximator(.28,.032,L232,3)->L232pns
Approximator(.28,.032,L233,3)->L233pns
Approximator(.28,.032,L234,3)->L234pns
Approximator(.28,.032,L235,3)->L235pns
Approximator(.28,.032,L237,3)->L237pns
Approximator(.28,.032,L239,3)->L239pns
Approximator(.28,.032,L240,3)->L240pns
Approximator(.28,.032,L241,3)->L241pns
Approximator(.28,.032,L242,3)->L242pns
Approximator(.28,.032,L243,3)->L243pns
Approximator(.28,.032,L244,3)->L244pns
Approximator(.28,.032,L245,3)->L245pns
Approximator(.28,.032,L246,3)->L246pns
Approximator(.28,.032,L247,3)->L247pns
Approximator(.28,.032,L248,3)->L248pns
Approximator(.28,.032,L249,3)->L249pns
Approximator(.28,.032,L250,3)->L250pns
Approximator(.28,.032,L251,3)->L251pns
Approximator(.28,.032,L252,3)->L252pns
Approximator(.28,.032,L253,3)->L253pns
Approximator(.28,.032,L255,3)->L255pns
Approximator(.28,.032,L256,3)->L256pns
Approximator(.28,.032,L257,3)->L257pns
Approximator(.28,.032,L258,3)->L258pns
Approximator(.28,.032,L260,3)->L260pns
Approximator(.28,.032,L261,3)->L261pns
Approximator(.28,.032,L262,3)->L262pns
Approximator(.28,.032,L264,3)->L264pns
Approximator(.28,.032,L265,3)->L265pns
Approximator(.28,.032,L266,3)->L266pns
Approximator(.28,.032,L267,3)->L267pns
Approximator(.28,.032,L269,3)->L269pns
Approximator(.28,.032,L270,3)->L270pns
Approximator(.28,.032,L271,3)->L271pns
Approximator(.28,.032,L272,3)->L272pns
Approximator(.28,.032,L275,3)->L275pns
Approximator(.28,.032,L276,3)->L276pns
Approximator(.28,.032,L277,3)->L277pns
Approximator(.28,.032,L279,3)->L279pns
Approximator(.28,.032,L280,3)->L280pns
Approximator(.28,.032,L281,3)->L281pns
Approximator(.28,.032,L282,3)->L282pns
Approximator(.28,.032,L288,3)->L288pns
Approximator(.28,.032,L289,3)->L289pns
Approximator(.28,.032,L290,3)->L290pns
Approximator(.28,.032,L292,3)->L292pns
Approximator(.28,.032,L293,3)->L293pns
Approximator(.28,.032,L294,3)->L294pns
Approximator(.28,.032,L295,3)->L295pns
Approximator(.28,.032,L296,3)->L296pns
Approximator(.28,.032,L297,3)->L297pns
Approximator(.28,.032,L298,3)->L298pns
Approximator(.28,.032,L301,3)->L301pns
Approximator(.28,.032,L302,3)->L302pns
Approximator(.28,.032,L303,3)->L303pns
Approximator(.28,.032,L305,3)->L305pns
Approximator(.28,.032,L306,3)->L306pns
Approximator(.28,.032,L307,3)->L307pns
Approximator(.28,.032,L310,3)->L310pns
Approximator(.28,.032,L311,3)->L311pns
Approximator(.28,.032,L312,3)->L312pns
Approximator(.28,.032,L314,3)->L314pns
Approximator(.28,.032,L315,3)->L315pns
Approximator(.28,.032,L316,3)->L316pns
Approximator(.28,.032,L317,3)->L317pns
Approximator(.28,.032,L318,3)->L318pns
Approximator(.28,.032,L319,3)->L319pns
Approximator(.28,.032,L320,3)->L320pns
Approximator(.28,.032,L321,3)->L321pns
Approximator(.28,.032,L322,3)->L322pns
Approximator(.28,.032,L325,3)->L325pns
Approximator(.28,.032,L327,3)->L327pns
Approximator(.28,.032,L328,3)->L328pns
Approximator(.28,.032,L329,3)->L329pns
Approximator(.28,.032,L330,3)->L330pns
Approximator(.28,.032,L331,3)->L331pns
Approximator(.28,.032,L333,3)->L333pns
Approximator(.28,.032,L334,3)->L334pns
Approximator(.28,.032,L335,3)->L335pns
Approximator(.28,.032,L336,3)->L336pns
Approximator(.28,.032,L337,3)->L337pns
Approximator(.28,.032,L338,3)->L338pns
Approximator(.28,.032,L339,3)->L339pns
Approximator(.28,.032,L340,3)->L340pns
Approximator(.28,.032,L341,3)->L341pns
Approximator(.28,.032,L342,3)->L342pns
Approximator(.28,.032,L343,3)->L343pns
Approximator(.28,.032,L344,3)->L344pns
Approximator(.28,.032,L346,3)->L346pns
Approximator(.28,.032,L347,3)->L347pns
Approximator(.28,.032,L348,3)->L348pns
Approximator(.28,.032,L349,3)->L349pns
out2<-rbind(L1pns,L2pns,L3pns,L4pns,L5pns,L6pns,L7pns,L8pns,L9pns,L10pns,L11pns,L12pns,L21pns,L22pns,L24pns,L27pns,L28pns,L30pns,L31pns,L32pns,L34pns,L35pns,L36pns,L37pns,L38pns,L39pns,L41pns,L42pns,L43pns,L44pns,L48pns,L49pns,L50pns,L51pns,L53pns,L54pns,L57pns,L58pns,L59pns,L61pns,L62pns, L63pns, L64pns, L65pns, L67pns, L72pns, L73pns, L74pns, L75pns, L76pns, L77pns, L78pns, L79pns, L80pns, L81pns, L82pns, L85pns, L86pns, L87pns, L88pns, L89pns, L91pns, L94pns, L95pns, L98pns, L99pns, L100pns,L102pns, L103pns, L104pns, L105pns, L107pns, L108pns, L110pns, L111pns, L113pns, L115pns, L116pns, L118pns, L119pns, L120pns, L121pns, L122pns, L125pns, L126pns, L127pns, L128pns, L129pns, L130pns, L131pns, L133pns, L134pns, L135pns, L136pns, L137pns, L141pns, L142pns, L143pns, L145pns, L146pns, L147pns, L148pns, L149pns, L150pns, L151pns, L152pns, L153pns, L154pns, L155pns, L157pns, L158pns, L159pns, L161pns, L162pns, L163pns, L165pns, L169pns, L170pns, L171pns, L174pns, L175pns, L176pns, L177pns, L178pns, L180pns, L181pns, L182pns, L183pns, L184pns, L186pns, L187pns, L188pns, L189pns, L190pns, L192pns, L193pns, L194pns, L195pns, L196pns, L197pns, L201pns, L202pns, L203pns, L204pns, L205pns, L206pns, L207pns, L208pns, L209pns, L210pns, L211pns, L213pns, L214pns, L215pns, L216pns, L217pns, L219pns, L221pns, L222pns, L223pns, L224pns, L225pns, L226pns, L228pns, L229pns, L230pns, L231pns, L232pns, L233pns, L234pns, L235pns, L237pns, L239pns, L240pns, L241pns, L242pns, L243pns, L244pns, L245pns, L246pns, L247pns, L248pns, L249pns, L250pns, L251pns, L252pns, L253pns, L255pns, L256pns, L257pns, L258pns, L260pns, L261pns, L262pns, L264pns, L265pns, L266pns, L267pns, L269pns, L270pns, L271pns, L272pns, L275pns, L276pns, L277pns, L279pns, L280pns, L281pns, L282pns, L288pns, L289pns, L290pns, L292pns, L293pns, L294pns, L295pns, L296pns, L297pns, L298pns, L301pns, L302pns, L303pns, L305pns, L306pns, L307pns, L310pns, L311pns, L312pns, L314pns, L315pns, L316pns, L317pns, L318pns, L319pns, L320pns, L321pns, L322pns, L325pns, L327pns, L328pns, L329pns, L330pns, L331pns, L333pns, L334pns, L335pns, L336pns, L337pns, L338pns, L339pns, L340pns, L341pns, L342pns, L343pns, L344pns, L346pns, L347pns, L348pns, L349pns)

write.csv(out2,file="ANALYSIS_FINAL_rates_1_349_5MY_branch")
hist(out2[,1])


####now we are changing the internode to 6 million years

Approximator(.28,.038,L1,3)->L1pns6million
Approximator(.28,.038,L2,3)->L2pns6million
Approximator(.28,.038,L3,3)->L3pns6million
Approximator(.28,.038,L4,3)->L4pns6million
Approximator(.28,.038,L5,3)->L5pns6million
Approximator(.28,.038,L6,3)->L6pns6million
Approximator(.28,.038,L7,3)->L7pns6million
Approximator(.28,.038,L8,3)->L8pns6million
Approximator(.28,.038,L9,3)->L9pns6million
Approximator(.28,.038,L10,3)->L10pns6million
Approximator(.28,.038,L11,3)->L11pns6million
Approximator(.28,.038,L12,3)->L12pns6million
Approximator(.28,.038,L21,3)->L21pns6million
Approximator(.28,.038,L22,3)->L22pns6million
Approximator(.28,.038,L24,3)->L24pns6million
Approximator(.28,.038,L27,3)->L27pns6million
Approximator(.28,.038,L28,3)->L28pns6million
Approximator(.28,.038,L30,3)->L30pns6million
Approximator(.28,.038,L31,3)->L31pns6million
Approximator(.28,.038,L32,3)->L32pns6million
Approximator(.28,.038,L34,3)->L34pns6million
Approximator(.28,.038,L35,3)->L35pns6million
Approximator(.28,.038,L36,3)->L36pns6million
Approximator(.28,.038,L37,3)->L37pns6million
Approximator(.28,.038,L38,3)->L38pns6million
Approximator(.28,.038,L39,3)->L39pns6million
Approximator(.28,.038,L41,3)->L41pns6million
Approximator(.28,.038,L42,3)->L42pns6million
Approximator(.28,.038,L43,3)->L43pns6million
Approximator(.28,.038,L44,3)->L44pns6million
Approximator(.28,.038,L48,3)->L48pns6million
Approximator(.28,.038,L49,3)->L49pns6million
Approximator(.28,.038,L50,3)->L50pns6million
Approximator(.28,.038,L51,3)->L51pns6million
Approximator(.28,.038,L53,3)->L53pns6million
Approximator(.28,.038,L54,3)->L54pns6million
Approximator(.28,.038,L57,3)->L57pns6million
Approximator(.28,.038,L58,3)->L58pns6million
Approximator(.28,.038,L59,3)->L59pns6million
Approximator(.28,.038,L61,3)->L61pns6million
Approximator(.28,.038,L62,3)->L62pns6million
Approximator(.28,.038,L63,3)->L63pns6million
Approximator(.28,.038,L64,3)->L64pns6million
Approximator(.28,.038,L65,3)->L65pns6million
Approximator(.28,.038,L67,3)->L67pns6million
Approximator(.28,.038,L72,3)->L72pns6million
Approximator(.28,.038,L73,3)->L73pns6million
Approximator(.28,.038,L74,3)->L74pns6million
Approximator(.28,.038,L75,3)->L75pns6million
Approximator(.28,.038,L76,3)->L76pns6million
Approximator(.28,.038,L77,3)->L77pns6million
Approximator(.28,.038,L78,3)->L78pns6million
Approximator(.28,.038,L79,3)->L79pns6million
Approximator(.28,.038,L80,3)->L80pns6million
Approximator(.28,.038,L81,3)->L81pns6million
Approximator(.28,.038,L82,3)->L82pns6million
Approximator(.28,.038,L85,3)->L85pns6million
Approximator(.28,.038,L86,3)->L86pns6million
Approximator(.28,.038,L87,3)->L87pns6million
Approximator(.28,.038,L88,3)->L88pns6million
Approximator(.28,.038,L89,3)->L89pns6million
Approximator(.28,.038,L91,3)->L91pns6million
Approximator(.28,.038,L94,3)->L94pns6million
Approximator(.28,.038,L95,3)->L95pns6million
Approximator(.28,.038,L98,3)->L98pns6million
Approximator(.28,.038,L99,3)->L99pns6million
Approximator(.28,.038,L100,3)->L100pns6million
Approximator(.28,.038,L102,3)->L102pns6million
Approximator(.28,.038,L103,3)->L103pns6million
Approximator(.28,.038,L104,3)->L104pns6million
Approximator(.28,.038,L105,3)->L105pns6million
Approximator(.28,.038,L107,3)->L107pns6million
Approximator(.28,.038,L108,3)->L108pns6million
Approximator(.28,.038,L110,3)->L110pns6million
Approximator(.28,.038,L111,3)->L111pns6million
Approximator(.28,.038,L113,3)->L113pns6million
Approximator(.28,.038,L115,3)->L115pns6million
Approximator(.28,.038,L116,3)->L116pns6million
Approximator(.28,.038,L118,3)->L118pns6million
Approximator(.28,.038,L119,3)->L119pns6million
Approximator(.28,.038,L120,3)->L120pns6million
Approximator(.28,.038,L121,3)->L121pns6million
Approximator(.28,.038,L122,3)->L122pns6million
Approximator(.28,.038,L125,3)->L125pns6million
Approximator(.28,.038,L126,3)->L126pns6million
Approximator(.28,.038,L127,3)->L127pns6million
Approximator(.28,.038,L128,3)->L128pns6million
Approximator(.28,.038,L129,3)->L129pns6million
Approximator(.28,.038,L130,3)->L130pns6million
Approximator(.28,.038,L131,3)->L131pns6million
Approximator(.28,.038,L133,3)->L133pns6million
Approximator(.28,.038,L134,3)->L134pns6million
Approximator(.28,.038,L135,3)->L135pns6million
Approximator(.28,.038,L136,3)->L136pns6million
Approximator(.28,.038,L137,3)->L137pns6million
Approximator(.28,.038,L141,3)->L141pns6million
Approximator(.28,.038,L142,3)->L142pns6million
Approximator(.28,.038,L143,3)->L143pns6million
Approximator(.28,.038,L145,3)->L145pns6million
Approximator(.28,.038,L146,3)->L146pns6million
Approximator(.28,.038,L147,3)->L147pns6million
Approximator(.28,.038,L148,3)->L148pns6million
Approximator(.28,.038,L149,3)->L149pns6million
Approximator(.28,.038,L150,3)->L150pns6million
Approximator(.28,.038,L151,3)->L151pns6million
Approximator(.28,.038,L152,3)->L152pns6million
Approximator(.28,.038,L153,3)->L153pns6million
Approximator(.28,.038,L154,3)->L154pns6million
Approximator(.28,.038,L155,3)->L155pns6million
Approximator(.28,.038,L157,3)->L157pns6million
Approximator(.28,.038,L158,3)->L158pns6million
Approximator(.28,.038,L159,3)->L159pns6million
Approximator(.28,.038,L161,3)->L161pns6million
Approximator(.28,.038,L162,3)->L162pns6million
Approximator(.28,.038,L163,3)->L163pns6million
Approximator(.28,.038,L165,3)->L165pns6million
Approximator(.28,.038,L169,3)->L169pns6million
Approximator(.28,.038,L170,3)->L170pns6million
Approximator(.28,.038,L171,3)->L171pns6million
Approximator(.28,.038,L174,3)->L174pns6million
Approximator(.28,.038,L175,3)->L175pns6million
Approximator(.28,.038,L176,3)->L176pns6million
Approximator(.28,.038,L177,3)->L177pns6million
Approximator(.28,.038,L178,3)->L178pns6million
Approximator(.28,.038,L180,3)->L180pns6million
Approximator(.28,.038,L181,3)->L181pns6million
Approximator(.28,.038,L182,3)->L182pns6million
Approximator(.28,.038,L183,3)->L183pns6million
Approximator(.28,.038,L184,3)->L184pns6million
Approximator(.28,.038,L186,3)->L186pns6million
Approximator(.28,.038,L187,3)->L187pns6million
Approximator(.28,.038,L188,3)->L188pns6million
Approximator(.28,.038,L189,3)->L189pns6million
Approximator(.28,.038,L190,3)->L190pns6million
Approximator(.28,.038,L192,3)->L192pns6million
Approximator(.28,.038,L193,3)->L193pns6million
Approximator(.28,.038,L194,3)->L194pns6million
Approximator(.28,.038,L195,3)->L195pns6million
Approximator(.28,.038,L196,3)->L196pns6million
Approximator(.28,.038,L197,3)->L197pns6million
Approximator(.28,.038,L201,3)->L201pns6million
Approximator(.28,.038,L202,3)->L202pns6million
Approximator(.28,.038,L203,3)->L203pns6million
Approximator(.28,.038,L204,3)->L204pns6million
Approximator(.28,.038,L205,3)->L205pns6million
Approximator(.28,.038,L206,3)->L206pns6million
Approximator(.28,.038,L207,3)->L207pns6million
Approximator(.28,.038,L208,3)->L208pns6million
Approximator(.28,.038,L209,3)->L209pns6million
Approximator(.28,.038,L210,3)->L210pns6million
Approximator(.28,.038,L211,3)->L211pns6million
Approximator(.28,.038,L213,3)->L213pns6million
Approximator(.28,.038,L214,3)->L214pns6million
Approximator(.28,.038,L215,3)->L215pns6million
Approximator(.28,.038,L216,3)->L216pns6million
Approximator(.28,.038,L217,3)->L217pns6million
Approximator(.28,.038,L219,3)->L219pns6million
Approximator(.28,.038,L221,3)->L221pns6million
Approximator(.28,.038,L222,3)->L222pns6million
Approximator(.28,.038,L223,3)->L223pns6million
Approximator(.28,.038,L224,3)->L224pns6million
Approximator(.28,.038,L225,3)->L225pns6million
Approximator(.28,.038,L226,3)->L226pns6million
Approximator(.28,.038,L228,3)->L228pns6million
Approximator(.28,.038,L229,3)->L229pns6million
Approximator(.28,.038,L230,3)->L230pns6million
Approximator(.28,.038,L231,3)->L231pns6million
Approximator(.28,.038,L232,3)->L232pns6million
Approximator(.28,.038,L233,3)->L233pns6million
Approximator(.28,.038,L234,3)->L234pns6million
Approximator(.28,.038,L235,3)->L235pns6million
Approximator(.28,.038,L237,3)->L237pns6million
Approximator(.28,.038,L239,3)->L239pns6million
Approximator(.28,.038,L240,3)->L240pns6million
Approximator(.28,.038,L241,3)->L241pns6million
Approximator(.28,.038,L242,3)->L242pns6million
Approximator(.28,.038,L243,3)->L243pns6million
Approximator(.28,.038,L244,3)->L244pns6million
Approximator(.28,.038,L245,3)->L245pns6million
Approximator(.28,.038,L246,3)->L246pns6million
Approximator(.28,.038,L247,3)->L247pns6million
Approximator(.28,.038,L248,3)->L248pns6million
Approximator(.28,.038,L249,3)->L249pns6million
Approximator(.28,.038,L250,3)->L250pns6million
Approximator(.28,.038,L251,3)->L251pns6million
Approximator(.28,.038,L252,3)->L252pns6million
Approximator(.28,.038,L253,3)->L253pns6million
Approximator(.28,.038,L255,3)->L255pns6million
Approximator(.28,.038,L256,3)->L256pns6million
Approximator(.28,.038,L257,3)->L257pns6million
Approximator(.28,.038,L258,3)->L258pns6million
Approximator(.28,.038,L260,3)->L260pns6million
Approximator(.28,.038,L261,3)->L261pns6million
Approximator(.28,.038,L262,3)->L262pns6million
Approximator(.28,.038,L264,3)->L264pns6million
Approximator(.28,.038,L265,3)->L265pns6million
Approximator(.28,.038,L266,3)->L266pns6million
Approximator(.28,.038,L267,3)->L267pns6million
Approximator(.28,.038,L269,3)->L269pns6million
Approximator(.28,.038,L270,3)->L270pns6million
Approximator(.28,.038,L271,3)->L271pns6million
Approximator(.28,.038,L272,3)->L272pns6million
Approximator(.28,.038,L275,3)->L275pns6million
Approximator(.28,.038,L276,3)->L276pns6million
Approximator(.28,.038,L277,3)->L277pns6million
Approximator(.28,.038,L279,3)->L279pns6million
Approximator(.28,.038,L280,3)->L280pns6million
Approximator(.28,.038,L281,3)->L281pns6million
Approximator(.28,.038,L282,3)->L282pns6million
Approximator(.28,.038,L288,3)->L288pns6million
Approximator(.28,.038,L289,3)->L289pns6million
Approximator(.28,.038,L290,3)->L290pns6million
Approximator(.28,.038,L292,3)->L292pns6million
Approximator(.28,.038,L293,3)->L293pns6million
Approximator(.28,.038,L294,3)->L294pns6million
Approximator(.28,.038,L295,3)->L295pns6million
Approximator(.28,.038,L296,3)->L296pns6million
Approximator(.28,.038,L297,3)->L297pns6million
Approximator(.28,.038,L298,3)->L298pns6million
Approximator(.28,.038,L301,3)->L301pns6million
Approximator(.28,.038,L302,3)->L302pns6million
Approximator(.28,.038,L303,3)->L303pns6million
Approximator(.28,.038,L305,3)->L305pns6million
Approximator(.28,.038,L306,3)->L306pns6million
Approximator(.28,.038,L307,3)->L307pns6million
Approximator(.28,.038,L310,3)->L310pns6million
Approximator(.28,.038,L311,3)->L311pns6million
Approximator(.28,.038,L312,3)->L312pns6million
Approximator(.28,.038,L314,3)->L314pns6million
Approximator(.28,.038,L315,3)->L315pns6million
Approximator(.28,.038,L316,3)->L316pns6million
Approximator(.28,.038,L317,3)->L317pns6million
Approximator(.28,.038,L318,3)->L318pns6million
Approximator(.28,.038,L319,3)->L319pns6million
Approximator(.28,.038,L320,3)->L320pns6million
Approximator(.28,.038,L321,3)->L321pns6million
Approximator(.28,.038,L322,3)->L322pns6million
Approximator(.28,.038,L325,3)->L325pns6million
Approximator(.28,.038,L327,3)->L327pns6million
Approximator(.28,.038,L328,3)->L328pns6million
Approximator(.28,.038,L329,3)->L329pns6million
Approximator(.28,.038,L330,3)->L330pns6million
Approximator(.28,.038,L331,3)->L331pns6million
Approximator(.28,.038,L333,3)->L333pns6million
Approximator(.28,.038,L334,3)->L334pns6million
Approximator(.28,.038,L335,3)->L335pns6million
Approximator(.28,.038,L336,3)->L336pns6million
Approximator(.28,.038,L337,3)->L337pns6million
Approximator(.28,.038,L338,3)->L338pns6million
Approximator(.28,.038,L339,3)->L339pns6million
Approximator(.28,.038,L340,3)->L340pns6million
Approximator(.28,.038,L341,3)->L341pns6million
Approximator(.28,.038,L342,3)->L342pns6million
Approximator(.28,.038,L343,3)->L343pns6million
Approximator(.28,.038,L344,3)->L344pns6million
Approximator(.28,.038,L346,3)->L346pns6million
Approximator(.28,.038,L347,3)->L347pns6million
Approximator(.28,.038,L348,3)->L348pns6million
Approximator(.28,.038,L349,3)->L349pns6million
out6million<-rbind(L1pns6million,L2pns6million,L3pns6million,L4pns6million,L5pns6million,L6pns6million,L7pns6million,L8pns6million,L9pns6million,L10pns6million,L11pns6million,L12pns6million,L21pns6million,L22pns6million,L24pns6million,L27pns6million,L28pns6million,L30pns6million,L31pns6million,L32pns6million,L34pns6million,L35pns6million,L36pns6million,L37pns6million,L38pns6million,L39pns6million,L41pns6million,L42pns6million,L43pns6million,L44pns6million,L48pns6million,L49pns6million,L50pns6million,L51pns6million,L53pns6million,L54pns6million,L57pns6million,L58pns6million,L59pns6million,L61pns6million,L62pns6million, L63pns6million, L64pns6million, L65pns6million, L67pns6million, L72pns6million, L73pns6million, L74pns6million, L75pns6million, L76pns6million, L77pns6million, L78pns6million, L79pns6million, L80pns6million, L81pns6million, L82pns6million, L85pns6million, L86pns6million, L87pns6million, L88pns6million, L89pns6million, L91pns6million, L94pns6million, L95pns6million, L98pns6million, L99pns6million, L100pns6million,L102pns6million, L103pns6million, L104pns6million, L105pns6million, L107pns6million, L108pns6million, L110pns6million, L111pns6million, L113pns6million, L115pns6million, L116pns6million, L118pns6million, L119pns6million, L120pns6million, L121pns6million, L122pns6million, L125pns6million, L126pns6million, L127pns6million, L128pns6million, L129pns6million, L130pns6million, L131pns6million, L133pns6million, L134pns6million, L135pns6million, L136pns6million, L137pns6million, L141pns6million, L142pns6million, L143pns6million, L145pns6million, L146pns6million, L147pns6million, L148pns6million, L149pns6million, L150pns6million, L151pns6million, L152pns6million, L153pns6million, L154pns6million, L155pns6million, L157pns6million, L158pns6million, L159pns6million, L161pns6million, L162pns6million, L163pns6million, L165pns6million, L169pns6million, L170pns6million, L171pns6million, L174pns6million, L175pns6million, L176pns6million, L177pns6million, L178pns6million, L180pns6million, L181pns6million, L182pns6million, L183pns6million, L184pns6million, L186pns6million, L187pns6million, L188pns6million, L189pns6million, L190pns6million, L192pns6million, L193pns6million, L194pns6million, L195pns6million, L196pns6million, L197pns6million, L201pns6million, L202pns6million, L203pns6million, L204pns6million, L205pns6million, L206pns6million, L207pns6million, L208pns6million, L209pns6million, L210pns6million, L211pns6million, L213pns6million, L214pns6million, L215pns6million, L216pns6million, L217pns6million, L219pns6million, L221pns6million, L222pns6million, L223pns6million, L224pns6million, L225pns6million, L226pns6million, L228pns6million, L229pns6million, L230pns6million, L231pns6million, L232pns6million, L233pns6million, L234pns6million, L235pns6million, L237pns6million, L239pns6million, L240pns6million, L241pns6million, L242pns6million, L243pns6million, L244pns6million, L245pns6million, L246pns6million, L247pns6million, L248pns6million, L249pns6million, L250pns6million, L251pns6million, L252pns6million, L253pns6million, L255pns6million, L256pns6million, L257pns6million, L258pns6million, L260pns6million, L261pns6million, L262pns6million, L264pns6million, L265pns6million, L266pns6million, L267pns6million, L269pns6million, L270pns6million, L271pns6million, L272pns6million, L275pns6million, L276pns6million, L277pns6million, L279pns6million, L280pns6million, L281pns6million, L282pns6million, L288pns6million, L289pns6million, L290pns6million, L292pns6million, L293pns6million, L294pns6million, L295pns6million, L296pns6million, L297pns6million, L298pns6million, L301pns6million, L302pns6million, L303pns6million, L305pns6million, L306pns6million, L307pns6million, L310pns6million, L311pns6million, L312pns6million, L314pns6million, L315pns6million, L316pns6million, L317pns6million, L318pns6million, L319pns6million, L320pns6million, L321pns6million, L322pns6million, L325pns6million, L327pns6million, L328pns6million, L329pns6million, L330pns6million, L331pns6million, L333pns6million, L334pns6million, L335pns6million, L336pns6million, L337pns6million, L338pns6million, L339pns6million, L340pns6million, L341pns6million, L342pns6million, L343pns6million, L344pns6million, L346pns6million, L347pns6million, L348pns6million, L349pns6million)
write.csv(out6million,file="ANALYSIS_FINAL_rates_1_349_6MY_branch")
hist(out6million[,1])




##### Now we can look at the partitions, first assemble them and convert them to a matrix format
part1<-c(L4,L1,L348,L305)
part2<-c(L10,L171,L269,L197,L79)
part3<-c(L216,L111,L100)
part4<-c(L102,L280,L75,L73,L303,L205,L234)
part5<-c(L50,L181,L178,L103,L86)
part6<-c(L232,L104,L39,L59,L242)
part7<-c(L105,L193)
part8<-c(L241,L41,L107,L257)
part9<-c(L2,L180,L310,L120,L49,L78,L108,L22)
part10<-c(L189,L219,L88,L9,L11)
part11<-c(L110,L331,L37)
part12<-c(L131,L113,L36,L64,L188,L290)
part13<-c(L217,L115,L210)
part14<-c(L297,L116,L82,L255,L134,L195,L233)
part15<-c(L277,L209,L141,L118,L121)
part16<-c(L266,L248,L119,L252,L250,L165,L147)
part17<-c(L12,L244,L175)
part18<-c(L122,L61)
part19<-c(L125,L319,L157)
part20<-c(L21,L126)
part21<-c(L228,L339,L230,L349,L130,L127)
part22<-c(L314,L272,L264,L128,L261,L30)
part23<-c(L129,L28,L162)
part24<-c(L196,L176,L229,L235,L133,L320,L44,L43,L341)
part25<-c(L32,L243,L135)
part26<-c(L246,L251,L38,L136,L311,L146)
part27<-c(L282,L137,L148,L245,L192,L294)
part28<-c(L317,L333,L204,L142,L337)
part29<-c(L174,L279,L143)
part30<-c(L145)
part31<-c(L312,L57,L256,L91,L149,L186,L270,L67)
part32<-c(L31,L150)
part33<-c(L315,L267,L151)
part34<-c(L58,L98,L152)
part35<-c(L153,L155,L258,L213,L95,L177)
part36<-c(L35,L316,L154,L169)
part37<-c(L158,L62)
part38<-c(L159,L343,L99,L336)
part39<-c(L211,L338,L292,L247,L161)
part40<-c(L163)
part41<-c(L53,L295,L340,L276,L281,L170,L302)
part42<-c(L231,L190,L182,L187)
part43<-c(L183)
part44<-c(L184)
part45<-c(L80,L194,L325,L87)
part46<-c(L201,L34,L74)
part47<-c(L48,L202)
part48<-c(L293,L203,L322)
part49<-c(L6,L327,L206,L249)
part50<-c(L24,L207,L224,L27)
part51<-c(L208,L239,L307)
part52<-c(L214)
part53<-c(L215,L237)
part54<-c(L89,L221)
part55<-c(L42,L222,L7)
part56<-c(L328,L76,L223)
part57<-c(L271,L225)
part58<-c(L226,L329)
part59<-c(L296,L240,L346,L347)
part60<-c(L321,L253,L330)
part61<-c(L260)
part62<-c(L334,L298,L72,L262)
part63<-c(L265)
part64<-c(L275,L51)
part65<-c(L288)
part66<-c(L342,L289)
part67<-c(L8,L54,L3,L94)
part68<-c(L301)
part69<-c(L306)
part70<-c(L5,L318)
part71<-c(L335)
part72<-c(L344)
part73<-c(L63,L65)
part74<-c(L85,L77)
part75<-c(L81)

as.matrix(part1)->part1
as.matrix(part2)->part2
as.matrix(part3)->part3
as.matrix(part4)->part4
as.matrix(part5)->part5
as.matrix(part6)->part6
as.matrix(part7)->part7
as.matrix(part8)->part8
as.matrix(part9)->part9
as.matrix(part10)->part10
as.matrix(part11)->part11
as.matrix(part12)->part12
as.matrix(part13)->part13
as.matrix(part14)->part14
as.matrix(part15)->part15
as.matrix(part16)->part16
as.matrix(part17)->part17
as.matrix(part18)->part18
as.matrix(part19)->part19
as.matrix(part20)->part20
as.matrix(part21)->part21
as.matrix(part22)->part22
as.matrix(part23)->part23
as.matrix(part24)->part24
as.matrix(part25)->part25
as.matrix(part26)->part26
as.matrix(part27)->part27
as.matrix(part28)->part28
as.matrix(part29)->part29
as.matrix(part30)->part30
as.matrix(part31)->part31
as.matrix(part32)->part32
as.matrix(part33)->part33
as.matrix(part34)->part34
as.matrix(part35)->part35
as.matrix(part36)->part36
as.matrix(part37)->part37
as.matrix(part38)->part38
as.matrix(part39)->part39
as.matrix(part40)->part40
as.matrix(part41)->part41
as.matrix(part42)->part42
as.matrix(part43)->part43
as.matrix(part44)->part44
as.matrix(part45)->part45
as.matrix(part46)->part46
as.matrix(part47)->part47
as.matrix(part48)->part48
as.matrix(part49)->part49
as.matrix(part50)->part50
as.matrix(part51)->part51
as.matrix(part52)->part52
as.matrix(part53)->part53
as.matrix(part54)->part54
as.matrix(part55)->part55
as.matrix(part56)->part56
as.matrix(part57)->part57
as.matrix(part58)->part58
as.matrix(part59)->part59
as.matrix(part60)->part60
as.matrix(part61)->part61
as.matrix(part62)->part62
as.matrix(part63)->part63
as.matrix(part64)->part64
as.matrix(part65)->part65
as.matrix(part66)->part66
as.matrix(part67)->part67
as.matrix(part68)->part68
as.matrix(part69)->part69
as.matrix(part70)->part70
as.matrix(part71)->part71
as.matrix(part72)->part72
as.matrix(part73)->part73
as.matrix(part74)->part74
as.matrix(part75)->part75

Approximator(.28,.038,part1,3)->part1pns6million
Approximator(.28,.038,part2,3)->part2pns6million
Approximator(.28,.038,part3,3)->part3pns6million
Approximator(.28,.038,part4,3)->part4pns6million
Approximator(.28,.038,part5,3)->part5pns6million
Approximator(.28,.038,part6,3)->part6pns6million
Approximator(.28,.038,part7,3)->part7pns6million
Approximator(.28,.038,part8,3)->part8pns6million
Approximator(.28,.038,part9,3)->part9pns6million
Approximator(.28,.038,part10,3)->part10pns6million
Approximator(.28,.038,part11,3)->part11pns6million
Approximator(.28,.038,part12,3)->part12pns6million
Approximator(.28,.038,part13,3)->part13pns6million
Approximator(.28,.038,part14,3)->part14pns6million
Approximator(.28,.038,part15,3)->part15pns6million
Approximator(.28,.038,part16,3)->part16pns6million
Approximator(.28,.038,part17,3)->part17pns6million
Approximator(.28,.038,part18,3)->part18pns6million
Approximator(.28,.038,part19,3)->part19pns6million
Approximator(.28,.038,part20,3)->part20pns6million
Approximator(.28,.038,part21,3)->part21pns6million
Approximator(.28,.038,part22,3)->part22pns6million
Approximator(.28,.038,part23,3)->part23pns6million
Approximator(.28,.038,part24,3)->part24pns6million
Approximator(.28,.038,part25,3)->part25pns6million
Approximator(.28,.038,part26,3)->part26pns6million
Approximator(.28,.038,part27,3)->part27pns6million
Approximator(.28,.038,part28,3)->part28pns6million
Approximator(.28,.038,part29,3)->part29pns6million
Approximator(.28,.038,part30,3)->part30pns6million
Approximator(.28,.038,part31,3)->part31pns6million
Approximator(.28,.038,part32,3)->part32pns6million
Approximator(.28,.038,part33,3)->part33pns6million
Approximator(.28,.038,part34,3)->part34pns6million
Approximator(.28,.038,part35,3)->part35pns6million
Approximator(.28,.038,part36,3)->part36pns6million
Approximator(.28,.038,part37,3)->part37pns6million
Approximator(.28,.038,part38,3)->part38pns6million
Approximator(.28,.038,part39,3)->part39pns6million
Approximator(.28,.038,part40,3)->part40pns6million
Approximator(.28,.038,part41,3)->part41pns6million
Approximator(.28,.038,part42,3)->part42pns6million
Approximator(.28,.038,part43,3)->part43pns6million
Approximator(.28,.038,part44,3)->part44pns6million
Approximator(.28,.038,part45,3)->part45pns6million
Approximator(.28,.038,part46,3)->part46pns6million
Approximator(.28,.038,part47,3)->part47pns6million
Approximator(.28,.038,part48,3)->part48pns6million
Approximator(.28,.038,part49,3)->part49pns6million
Approximator(.28,.038,part50,3)->part50pns6million
Approximator(.28,.038,part51,3)->part51pns6million
Approximator(.28,.038,part52,3)->part52pns6million
Approximator(.28,.038,part53,3)->part53pns6million
Approximator(.28,.038,part54,3)->part54pns6million
Approximator(.28,.038,part55,3)->part55pns6million
Approximator(.28,.038,part56,3)->part56pns6million
Approximator(.28,.038,part57,3)->part57pns6million
Approximator(.28,.038,part58,3)->part58pns6million
Approximator(.28,.038,part59,3)->part59pns6million
Approximator(.28,.038,part60,3)->part60pns6million
Approximator(.28,.038,part61,3)->part61pns6million
Approximator(.28,.038,part62,3)->part62pns6million
Approximator(.28,.038,part63,3)->part63pns6million
Approximator(.28,.038,part64,3)->part64pns6million
Approximator(.28,.038,part65,3)->part65pns6million
Approximator(.28,.038,part66,3)->part66pns6million
Approximator(.28,.038,part67,3)->part67pns6million
Approximator(.28,.038,part68,3)->part68pns6million
Approximator(.28,.038,part69,3)->part69pns6million
Approximator(.28,.038,part70,3)->part70pns6million
Approximator(.28,.038,part71,3)->part71pns6million
Approximator(.28,.038,part72,3)->part72pns6million
Approximator(.28,.038,part73,3)->part73pns6million
Approximator(.28,.038,part74,3)->part74pns6million
Approximator(.28,.038,part75,3)->part75pns6million




out_Parts_6million<-rbind(part1pns6million,part2pns6million,part3pns6million,part4pns6million,part5pns6million,part6pns6million,part7pns6million,part8pns6million,part9pns6million,part10pns6million,part11pns6million,part12pns6million,part13pns6million,part14pns6million,part15pns6million,part16pns6million,part17pns6million,part18pns6million,part19pns6million,part20pns6million,part21pns6million,part22pns6million,part23pns6million,part24pns6million,part25pns6million,part26pns6million,part27pns6million,part28pns6million,part29pns6million,part30pns6million,part31pns6million,part32pns6million,part33pns6million,part34pns6million,part35pns6million,part36pns6million,part37pns6million,part38pns6million,part39pns6million,part40pns6million,part41pns6million,part42pns6million,part43pns6million,part44pns6million,part45pns6million,part46pns6million,part47pns6million,part48pns6million,part49pns6million,part50pns6million,part51pns6million,part52pns6million,part53pns6million,part54pns6million,part55pns6million, part56pns6million,part57pns6million,part58pns6million,part59pns6million,part60pns6million,part61pns6million,part62pns6million,part63pns6million,part64pns6million,part65pns6million,part66pns6million,part67pns6million,part68pns6million,part69pns6million,part70pns6million,part71pns6million,part72pns6million,part73pns6million,part74pns6million,part75pns6million)
write.csv(out_Parts_6million,file="ANALYSIS_FINAL_rates_parts_1_349_6MY_branch")
hist(out_Parts_6million[,1])

Approximator(.28,.032,part1,3)->part1pns555
Approximator(.28,.032,part2,3)->part2pns555
Approximator(.28,.032,part3,3)->part3pns555
Approximator(.28,.032,part4,3)->part4pns555
Approximator(.28,.032,part5,3)->part5pns555
Approximator(.28,.032,part6,3)->part6pns555
Approximator(.28,.032,part7,3)->part7pns555
Approximator(.28,.032,part8,3)->part8pns555
Approximator(.28,.032,part9,3)->part9pns555
Approximator(.28,.032,part10,3)->part10pns555
Approximator(.28,.032,part11,3)->part11pns555
Approximator(.28,.032,part12,3)->part12pns555
Approximator(.28,.032,part13,3)->part13pns555
Approximator(.28,.032,part14,3)->part14pns555
Approximator(.28,.032,part15,3)->part15pns555
Approximator(.28,.032,part16,3)->part16pns555
Approximator(.28,.032,part17,3)->part17pns555
Approximator(.28,.032,part18,3)->part18pns555
Approximator(.28,.032,part19,3)->part19pns555
Approximator(.28,.032,part20,3)->part20pns555
Approximator(.28,.032,part21,3)->part21pns555
Approximator(.28,.032,part22,3)->part22pns555
Approximator(.28,.032,part23,3)->part23pns555
Approximator(.28,.032,part24,3)->part24pns555
Approximator(.28,.032,part25,3)->part25pns555
Approximator(.28,.032,part26,3)->part26pns555
Approximator(.28,.032,part27,3)->part27pns555
Approximator(.28,.032,part28,3)->part28pns555
Approximator(.28,.032,part29,3)->part29pns555
Approximator(.28,.032,part30,3)->part30pns555
Approximator(.28,.032,part31,3)->part31pns555
Approximator(.28,.032,part32,3)->part32pns555
Approximator(.28,.032,part33,3)->part33pns555
Approximator(.28,.032,part34,3)->part34pns555
Approximator(.28,.032,part35,3)->part35pns555
Approximator(.28,.032,part36,3)->part36pns555
Approximator(.28,.032,part37,3)->part37pns555
Approximator(.28,.032,part38,3)->part38pns555
Approximator(.28,.032,part39,3)->part39pns555
Approximator(.28,.032,part40,3)->part40pns555
Approximator(.28,.032,part41,3)->part41pns555
Approximator(.28,.032,part42,3)->part42pns555
Approximator(.28,.032,part43,3)->part43pns555
Approximator(.28,.032,part44,3)->part44pns555
Approximator(.28,.032,part45,3)->part45pns555
Approximator(.28,.032,part46,3)->part46pns555
Approximator(.28,.032,part47,3)->part47pns555
Approximator(.28,.032,part48,3)->part48pns555
Approximator(.28,.032,part49,3)->part49pns555
Approximator(.28,.032,part50,3)->part50pns555
Approximator(.28,.032,part51,3)->part51pns555
Approximator(.28,.032,part52,3)->part52pns555
Approximator(.28,.032,part53,3)->part53pns555
Approximator(.28,.032,part54,3)->part54pns555
Approximator(.28,.032,part55,3)->part55pns555
Approximator(.28,.032,part56,3)->part56pns555
Approximator(.28,.032,part57,3)->part57pns555
Approximator(.28,.032,part58,3)->part58pns555
Approximator(.28,.032,part59,3)->part59pns555
Approximator(.28,.032,part60,3)->part60pns555
Approximator(.28,.032,part61,3)->part61pns555
Approximator(.28,.032,part62,3)->part62pns555
Approximator(.28,.032,part63,3)->part63pns555
Approximator(.28,.032,part64,3)->part64pns555
Approximator(.28,.032,part65,3)->part65pns555
Approximator(.28,.032,part66,3)->part66pns555
Approximator(.28,.032,part67,3)->part67pns555
Approximator(.28,.032,part68,3)->part68pns555
Approximator(.28,.032,part69,3)->part69pns555
Approximator(.28,.032,part70,3)->part70pns555
Approximator(.28,.032,part71,3)->part71pns555
Approximator(.28,.032,part72,3)->part72pns555
Approximator(.28,.032,part73,3)->part73pns555
Approximator(.28,.032,part74,3)->part74pns555
Approximator(.28,.032,part75,3)->part75pns555


out_Parts_555<-rbind(part1pns555,part2pns555,part3pns555,part4pns555,part5pns555,part6pns555,part7pns555,part8pns555,part9pns555,part10pns555,part11pns555,part12pns555,part13pns555,part14pns555,part15pns555,part16pns555,part17pns555,part18pns555,part19pns555,part20pns555,part21pns555,part22pns555,part23pns555,part24pns555,part25pns555,part26pns555,part27pns555,part28pns555,part29pns555,part30pns555,part31pns555,part32pns555,part33pns555,part34pns555,part35pns555,part36pns555,part37pns555,part38pns555,part39pns555,part40pns555,part41pns555,part42pns555,part43pns555,part44pns555,part45pns555,part46pns555,part47pns555,part48pns555,part49pns555,part50pns555,part51pns555,part52pns555,part53pns555,part54pns555,part55pns555, part56pns555,part57pns555,part58pns555,part59pns555,part60pns555,part61pns555,part62pns555,part63pns555,part64pns555,part65pns555,part66pns555,part67pns555,part68pns555,part69pns555,part70pns555,part71pns555,part72pns555,part73pns555,part74pns555,part75pns555)
write.csv(out_Parts_555,file="ANALYSIS_FINAL_rates_parts_1_349_5MY_branch")
hist(out_Parts_555[,1])


##Now let's generate PI profiles
#Note these files are just cleaned up versions of the output generated above.
read.csv("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/ALL_PI_Data.csv", header=TRUE)->PI
read.csv("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/ALL_PI_Names.csv")->PI_names
row.names(PI)<-PI_names[,1]

PI[1,]->branching_times

###also look at neoAves
read.csv("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/ALL_PI_Data_to_plot_NeoAves.csv", header=FALSE)->Navplot
dim(Navplot)

##first let's see what the worst loci in terms of Signal look like
##i did this externally really quickly using the output tables

#We have these loci as "worst": 349,272,193,208,260,214,223,152,305,344
##these are the "best": 184,81,266,213,150,121,338,216,182,149

##Here we are just doing a matrix transpose and naming the columns to get this formated correctly for PI functions
t(Navplot)->col_navplot
as.data.frame(col_navplot)->df_colnavplot
names(df_colnavplot)<-PI_names[,1]


##getting ready to plot, just some more housekeeping
round(max(df_colnavplot[,1]))->upper
upper/5->by.this

cbind(df_colnavplot$L349, df_colnavplot$L272,df_colnavplot$L193,df_colnavplot$L208,df_colnavplot$L260,df_colnavplot$L214,df_colnavplot$L223,df_colnavplot$L152,df_colnavplot$L305,df_colnavplot$L344)->worst

round(max(worst))->uppery
uppery/10->by.y


##ok now plot, note that I pull figures into illustrator anyway so I did not move the 0 of the x axis over to the left. If you are using graphics directly out of R you will need to do this.

plot(df_colnavplot[,1], worst[,1],pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))

axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
c("black","blue","gray","green","purple","brown","azure","red","yellow","pink")->colors
c("L349","L272","L193","L208","L260","L214","L223","L152","L305","L344")->leglab
c(1,2,3,1,2,3,1,2,3,1)->style
c(2,2,3,2,2,3,2,2,3,2)->thickness
legend("topright",y=NULL,leglab,lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
for (i in 1:10){
worst[,i]->inform.at.time	
yy <-predict(interpSpline(df_colnavplot[,1], inform.at.time))

lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)

}

##repeat with the highest prob loci
cbind(df_colnavplot$L184, df_colnavplot$L81,df_colnavplot$L266,df_colnavplot$L213,df_colnavplot$L150,df_colnavplot$L121,df_colnavplot$L338,df_colnavplot$L216,df_colnavplot$L182,df_colnavplot$L149)->best

round(max(best))->uppery
uppery/10->by.y

plot(df_colnavplot[,1], worst[,1],pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))

axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
c("black","blue","gray","green","purple","brown","azure","red","yellow","pink")->colors
c("L184","L81","L266","L213","L150","L121","L338","L216","L182","L149")->leglab
c(1,2,3,1,2,3,1,2,3,1)->style
c(2,2,3,2,2,3,2,2,3,2)->thickness
legend("topright",y=NULL,leglab,lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
for (i in 1:10){
best[,i]->inform.at.time	
yy <-predict(interpSpline(df_colnavplot[,1], inform.at.time))

lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)

}

###ok lets see which profiles have the lowest ratio of peak to root informativeness, again this is just cleaned up output

read.csv("~/Dropbox/avian-phylogenomics/Zenodo/Phylogenetic_Informativeness/ALL_PI_Data_to_plot.csv", header=FALSE)->ALLplot
dim(ALLplot)
t(ALLplot)->col_ALLplot
as.data.frame(col_ALLplot)->df_ALLnavplot
names(df_ALLnavplot)<-PI_names[,1]

dim(df_ALLnavplot)


###okie dokie lets output the profiles of the partitions, i'm just going to do this quick here using greps in textwrangler rather than writing an efficient function. I'm using the profile generator function in the sourced functions, inferring each profile and then writing them into a text file so we can use the same workflow as above

##first part might take a couple of minutes

inform.profile.generator2(part1,tree)->part1pnsprofile_out
inform.profile.generator2(part2,tree)->part2pnsprofile_out
inform.profile.generator2(part3,tree)->part3pnsprofile_out
inform.profile.generator2(part4,tree)->part4pnsprofile_out
inform.profile.generator2(part5,tree)->part5pnsprofile_out
inform.profile.generator2(part6,tree)->part6pnsprofile_out
inform.profile.generator2(part7,tree)->part7pnsprofile_out
inform.profile.generator2(part8,tree)->part8pnsprofile_out
inform.profile.generator2(part9,tree)->part9pnsprofile_out
inform.profile.generator2(part10,tree)->part10pnsprofile_out
inform.profile.generator2(part11,tree)->part11pnsprofile_out
inform.profile.generator2(part12,tree)->part12pnsprofile_out
inform.profile.generator2(part13,tree)->part13pnsprofile_out
inform.profile.generator2(part14,tree)->part14pnsprofile_out
inform.profile.generator2(part15,tree)->part15pnsprofile_out
inform.profile.generator2(part16,tree)->part16pnsprofile_out
inform.profile.generator2(part17,tree)->part17pnsprofile_out
inform.profile.generator2(part18,tree)->part18pnsprofile_out
inform.profile.generator2(part19,tree)->part19pnsprofile_out
inform.profile.generator2(part20,tree)->part20pnsprofile_out
inform.profile.generator2(part21,tree)->part21pnsprofile_out
inform.profile.generator2(part22,tree)->part22pnsprofile_out
inform.profile.generator2(part23,tree)->part23pnsprofile_out
inform.profile.generator2(part24,tree)->part24pnsprofile_out
inform.profile.generator2(part25,tree)->part25pnsprofile_out
inform.profile.generator2(part26,tree)->part26pnsprofile_out
inform.profile.generator2(part27,tree)->part27pnsprofile_out
inform.profile.generator2(part28,tree)->part28pnsprofile_out
inform.profile.generator2(part29,tree)->part29pnsprofile_out
inform.profile.generator2(part30,tree)->part30pnsprofile_out
inform.profile.generator2(part31,tree)->part31pnsprofile_out
inform.profile.generator2(part32,tree)->part32pnsprofile_out
inform.profile.generator2(part33,tree)->part33pnsprofile_out
inform.profile.generator2(part34,tree)->part34pnsprofile_out
inform.profile.generator2(part35,tree)->part35pnsprofile_out
inform.profile.generator2(part36,tree)->part36pnsprofile_out
inform.profile.generator2(part37,tree)->part37pnsprofile_out
inform.profile.generator2(part38,tree)->part38pnsprofile_out
inform.profile.generator2(part39,tree)->part39pnsprofile_out
inform.profile.generator2(part40,tree)->part40pnsprofile_out
inform.profile.generator2(part41,tree)->part41pnsprofile_out
inform.profile.generator2(part42,tree)->part42pnsprofile_out
inform.profile.generator2(part43,tree)->part43pnsprofile_out
inform.profile.generator2(part44,tree)->part44pnsprofile_out
inform.profile.generator2(part45,tree)->part45pnsprofile_out
inform.profile.generator2(part46,tree)->part46pnsprofile_out
inform.profile.generator2(part47,tree)->part47pnsprofile_out
inform.profile.generator2(part48,tree)->part48pnsprofile_out
inform.profile.generator2(part49,tree)->part49pnsprofile_out
inform.profile.generator2(part50,tree)->part50pnsprofile_out
inform.profile.generator2(part51,tree)->part51pnsprofile_out
inform.profile.generator2(part52,tree)->part52pnsprofile_out
inform.profile.generator2(part53,tree)->part53pnsprofile_out
inform.profile.generator2(part54,tree)->part54pnsprofile_out
inform.profile.generator2(part55,tree)->part55pnsprofile_out
inform.profile.generator2(part56,tree)->part56pnsprofile_out
inform.profile.generator2(part57,tree)->part57pnsprofile_out
inform.profile.generator2(part58,tree)->part58pnsprofile_out
inform.profile.generator2(part59,tree)->part59pnsprofile_out
inform.profile.generator2(part60,tree)->part60pnsprofile_out
inform.profile.generator2(part61,tree)->part61pnsprofile_out
inform.profile.generator2(part62,tree)->part62pnsprofile_out
inform.profile.generator2(part63,tree)->part63pnsprofile_out
inform.profile.generator2(part64,tree)->part64pnsprofile_out
inform.profile.generator2(part65,tree)->part65pnsprofile_out
inform.profile.generator2(part66,tree)->part66pnsprofile_out
inform.profile.generator2(part67,tree)->part67pnsprofile_out
inform.profile.generator2(part68,tree)->part68pnsprofile_out
inform.profile.generator2(part69,tree)->part69pnsprofile_out
inform.profile.generator2(part70,tree)->part70pnsprofile_out
inform.profile.generator2(part71,tree)->part71pnsprofile_out
inform.profile.generator2(part72,tree)->part72pnsprofile_out
inform.profile.generator2(part73,tree)->part73pnsprofile_out
inform.profile.generator2(part74,tree)->part74pnsprofile_out
inform.profile.generator2(part75,tree)->part75pnsprofile_out


##assemble!
rbind(part1pnsprofile_out, part2pnsprofile_out, part3pnsprofile_out, part4pnsprofile_out, part5pnsprofile_out, part6pnsprofile_out, part7pnsprofile_out, part8pnsprofile_out, part9pnsprofile_out, part10pnsprofile_out, part11pnsprofile_out, part12pnsprofile_out, part13pnsprofile_out, part14pnsprofile_out, part15pnsprofile_out, part16pnsprofile_out, part17pnsprofile_out, part18pnsprofile_out, part19pnsprofile_out, part20pnsprofile_out, part21pnsprofile_out, part22pnsprofile_out, part23pnsprofile_out, part24pnsprofile_out, part25pnsprofile_out, part26pnsprofile_out, part27pnsprofile_out, part28pnsprofile_out, part29pnsprofile_out, part30pnsprofile_out, part31pnsprofile_out, part32pnsprofile_out, part33pnsprofile_out, part34pnsprofile_out, part35pnsprofile_out, part36pnsprofile_out, part37pnsprofile_out, part38pnsprofile_out, part39pnsprofile_out, part40pnsprofile_out, part41pnsprofile_out, part42pnsprofile_out, part43pnsprofile_out, part44pnsprofile_out, part45pnsprofile_out, part46pnsprofile_out, part47pnsprofile_out, part48pnsprofile_out, part49pnsprofile_out, part50pnsprofile_out, part51pnsprofile_out, part52pnsprofile_out, part53pnsprofile_out, part54pnsprofile_out, part55pnsprofile_out, part56pnsprofile_out, part57pnsprofile_out, part58pnsprofile_out, part59pnsprofile_out, part60pnsprofile_out, part61pnsprofile_out, part62pnsprofile_out, part63pnsprofile_out, part64pnsprofile_out, part65pnsprofile_out, part66pnsprofile_out, part67pnsprofile_out, part68pnsprofile_out, part69pnsprofile_out, part70pnsprofile_out, part71pnsprofile_out, part72pnsprofile_out, part73pnsprofile_out, part74pnsprofile_out, part75pnsprofile_out)->profiles_of_partitions

write.table(profiles_of_partitions, file="profilebyPARTITION", quote=FALSE, append=TRUE)




### lets get a sense of the variance!

profiles_of_partitions[,204]->root_inform
profiles_of_partitions[,193]->neoaves_inform
##lets look at the ratio of informativeness at the crown neoaves versus the root of aves: near 1 is constant, greater than 1 is a decline, less than 1 is a rise
neoaves_inform/root_inform->ratio
rationames<-c("part1","part2","part3","part4","part5","part6","part7","part8","part9","part10","part11","part12","part13","part14","part15","part16","part17","part18","part19","part20","part21","part22","part23","part24","part25","part26","part27","part28","part29","part30","part31","part32","part33","part34","part35","part36","part37","part38","part39","part40","part41","part42","part43","part44","part45","part46","part47","part48","part49","part50","part51","part52","part53","part54","part55","part56","part57","part58","part59","part60","part61","part62","part63","part64","part65","part66","part67","part68","part69","part70","part71","part72","part73","part74","part75")
names(ratio)<-rationames
sort(ratio)
##top partitions
##part61    part37    part13    part10    part72    part17 
##0.8943941 1.0043514 1.0312579 1.0419695 1.0650173 1.1024419 
##   part60    part22    part33    part30    part44    part12 
##1.1064528 1.1348308 1.1349526 1.1582414 1.1584780 1.1675807 
##   part16    part75    part38    part39    part74    part48 
##1.1692707 1.1815259 1.1846315 1.2081084 1.2093177 1.2130554 
##   part26    part34    part32    part29    part47    part31 
##1.2234491 1.2263921 1.2311682 1.2326462 1.2537351 1.2639625 
##    part4    part42    part24    part71    part11    part27 
##1.2721031 1.2832110 1.2936762 1.2994874 1.3078838 1.3111915 
##   part73    part45    part28    part15    part35    part36 
##1.3112310 1.3129048 1.3162551 1.3251077 1.3363793 1.3461475 
##   part46    part23    part53    part18     part3    part40 
##1.3465043 1.3588619 1.3728418 1.3741448 1.3779201 1.3942025 
##   part21    part69     part8     part6     part5    part41 
##1.4063893 1.4212117 1.4288790 1.4336324 1.4350433 1.4369905 
##   part43    part19    part64    part66    part57    part49 
##1.4425234 1.4599935 1.4662153 1.4777498 1.4877277 1.4911740 
##   part65    part54    part14    part50    part67    part59 
##1.4984127 1.5046317 1.5229979 1.5309421 1.5414402 1.5497066 
##   part70     part9     part2    part52    part25    part55 
##1.5655939 1.5713477 1.5860405 1.5965080 1.6065983 1.6383587 
##   part56    part63    part20    part68    part62    part58 
##1.6757564 1.7178926 1.7239050 1.7565023 1.7677280 1.8230559 
##    part7    part51     part1 
##1.8821471 1.9929689 2.1133824 








rbind(part1pnsprofile_out, part51pnsprofile_out, part7pnsprofile_out, part58pnsprofile_out, part62pnsprofile_out, part68pnsprofile_out, part20pnsprofile_out, part63pnsprofile_out, part56pnsprofile_out, part55pnsprofile_out)->highest
rbind(part34pnsprofile_out,part29pnsprofile_out,part4pnsprofile_out,part71pnsprofile_out, part28pnsprofile_out, part36pnsprofile_out, part23pnsprofile_out, part21pnsprofile_out,part69pnsprofile_out,part8pnsprofile_out)->midrange
rbind(part61pnsprofile_out, part37pnsprofile_out, part13pnsprofile_out, part10pnsprofile_out, part72pnsprofile_out, part17pnsprofile_out, part60pnsprofile_out, part22pnsprofile_out, part33pnsprofile_out, part30pnsprofile_out)->lowest


##plot!

#first transpose the above
t(highest)->tprof
t(midrange)->tprof2
t(lowest)->tprof3



round(max(highest))->uppery
uppery/10->by.y
sort(branching.times(tree))->x
x[1:203]->x
c(0,x)->x

max(x)->upper
upper/10-> by.this
tprof[,1]->temp
temp[1:204]->temp2
plot(x, temp2,pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))

axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
c("black","blue","gray","green","purple","brown","azure","red","yellow","pink")->colors
c("part1","part51","part7","part58","part62","part68","part20","part63","part56","part55")->leglab
c(1,2,3,1,2,3,1,2,3,1)->style
c(2,2,3,2,2,3,2,2,3,2)->thickness
legend("topright",y=NULL,leglab,lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
for (i in 1:10){
tprof[,i]->inform.at.time
inform.at.time[1:204]->crown.aves	

yy <-predict(interpSpline(x, crown.aves))

lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)

}

##midrange
round(max(midrange))->uppery
uppery/10->by.y
tprof2[,1]->temp
temp[1:204]->temp2
plot(x, temp2,pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))
axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
c("black","blue","gray","green","purple","brown","azure","red","yellow","pink")->colors
c("part34","part29","part4","part71","part28","part36","part23","part21","part69","part8")->leglab
c(1,2,3,1,2,3,1,2,3,1)->style
c(2,2,3,2,2,3,2,2,3,2)->thickness
legend("topright",y=NULL,leglab,lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
for (i in 1:10){
tprof2[,i]->inform.at.time
inform.at.time[1:204]->crown.aves	

yy <-predict(interpSpline(x, crown.aves))

lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)

}

##best
round(max(lowest))->uppery
uppery/10->by.y
tprof3[,1]->temp
temp[1:204]->temp2
plot(x, temp2,pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))
axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
c("black","blue","gray","green","purple","brown","azure","red","yellow","pink")->colors
c("part61","part37","part13","part10","part72","part17","part60","part22","part33","part30")->leglab
c(1,2,3,1,2,3,1,2,3,1)->style
c(2,2,3,2,2,3,2,2,3,2)->thickness
legend("topright",y=NULL,leglab,lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
for (i in 1:10){
tprof3[,i]->inform.at.time
inform.at.time[1:204]->crown.aves	

yy <-predict(interpSpline(x, crown.aves))

lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)

}


##lets try the whole thing!! Crazy!
all<-c(part1,part2,part3,part4,part5,part6,part7,part8,part9,part10,part11,part12,part13,part14,part15,part16,part17,part18,part19,part20,part21,part22,part23,part24,part25,part26,part27,part28,part29,part30,part31,part32,part33,part34,part35,part36,part37,part38,part39,part40,part41,part42,part43,part44,part45,part46,part47,part48,part49,part50,part51,part52,part53,part54,part55,part56,part57,part58,part59,part60,part61,part62,part63,part64,part65,part66,part67,part68,part69,part70,part71,part72,part73,part74,part75)
top65<-c(part2,part3,part4,part5,part6,part8,part9,part10,part11,part12,part13,part14,part15,part16,part17,part18,part19,part21,part22,part23,part24,part25,part26,part27,part28,part29,part30,part31,part32,part33,part34,part35,part36,part37,part38,part39,part40,part41,part42,part43,part44,part45,part46,part47,part48,part49,part50,part52,part53,part54,part57,part59,part60,part61,part64,part65,part66,part67,part69,part70,part71,part72,part73,part74,part75)
inform.profile.generator2(top65,tree)->mega_pruned_profile
write.table(mega_pruned_profile, file="profile_of_pruneddataset", quote=FALSE)
t(megaprofile)->tprof4
t(mega_pruned_profile)->tprof5

round(max(tprof4))->uppery
uppery/10->by.y
tprof4[,1]->temp
temp[1:204]->temp2
plot(x, temp2,pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))
axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
c("black","blue","gray","green","purple","brown","azure","red","yellow","pink")->colors
c("all_loci","part37","part13","part10","part72","part17","part60","part22","part33","part30")->leglab
c(1,2)->style
c(2)->thickness
legend("topright",y=NULL,leglab,lty=style,col=colors,lwd=thickness,title="Partition PI Profile")


yy <-predict(interpSpline(x, temp2))

lines(yy, pch=NA_integer_, col=colors[1],lty=style[1],)

tprof5[,1]->temp3
temp3[1:204]->temp4
yy <-predict(interpSpline(x, temp4))

lines(yy, pch=NA_integer_, col=colors[2],lty=style[2],)


##those plots werent that useful at all, lets see what the dive time profiles look like as a whole



rbind(part3pnsprofile_out,part11pnsprofile_out,part13pnsprofile_out, part17pnsprofile_out, part18pnsprofile_out, part19pnsprofile_out, part29pnsprofile_out,part30pnsprofile_out,part32pnsprofile_out, part33pnsprofile_out, part34pnsprofile_out, part40pnsprofile_out, part43pnsprofile_out, part44pnsprofile_out, part69pnsprofile_out, part70pnsprofile_out, part71pnsprofile_out)->panelA
rbind(part46pnsprofile_out,part47pnsprofile_out,part48pnsprofile_out,part13pnsprofile_out, part52pnsprofile_out, part53pnsprofile_out, part54pnsprofile_out, part57pnsprofile_out,part58pnsprofile_out,part60pnsprofile_out, part61pnsprofile_out, part64pnsprofile_out, part65pnsprofile_out, part66pnsprofile_out, part68pnsprofile_out, part72pnsprofile_out, part73pnsprofile_out, part74pnsprofile_out, part75pnsprofile_out)->panelB


##plot!






#first transpose the above
t(panelA)->tdiv
t(panelB)->tdiv2

round(max(tdiv))->uppery
uppery/10->by.y
tdiv[,1]->temp
temp[1:204]->temp2
plot(x, temp2,pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))
axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
c("black","blue","gray","green","purple","brown","azure","red","yellow","pink","black","blue","gray","green","purple","brown","azure")->colors
c("part3","part11","part13", "part17", "part18", "part19", "part29","part30","part32", "part33", "part34", "part40"," part43", "part44", "part69", "part70", "part71")->leglab
c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1)->style
c(2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2,2,2,3)->thickness
legend("topright",y=NULL,leglab,lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
for (i in 1:length(leglab)){
tdiv[,i]->inform.at.time
inform.at.time[1:204]->crown.aves	

yy <-predict(interpSpline(x, crown.aves))

lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)

}




round(max(tdiv2))->uppery
uppery/10->by.y
tdiv2[,1]->temp
temp[1:204]->temp2
plot(x, temp2,pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))
axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
c("black","blue","gray","green","purple","brown","azure","red","yellow","pink","black","blue","gray","green","purple","brown","azure","red","yellow","pink")->colors
c("part46","part47","part48","part13","part52","part53","part54","part57","part58","part60","part61","part64","part65","part66","part68","part72","part73","part74","part75")->leglab
c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1)->style
c(2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2)->thickness
legend("topright",y=NULL,leglab,lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
for (i in 1:length(leglab)){
tdiv2[,i]->inform.at.time
inform.at.time[1:204]->crown.aves	

yy <-predict(interpSpline(x, crown.aves))

lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)

}



###OK, now lets make some heat maps! Just repeat the below for different versions and use space.maker.narrow for a zoomed in perspective

allrates<-c(L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L21,L22,L24,L27,L28,L30,L31,L32,L34,L35,L36,L37,L38,L39,L41,L42,L43,L44,L48,L49,L50,L51,L53,L54,L57,L58,L59,L61,L62,L63,L64,L65,L67,L72,L73,L74,L75,L76,L77,L78,L79,L80,L81,L82,L85,L86,L87,L88,L89,L91,L94,L95,L98,L99,L100,L102,L103,L104,L105,L107,L108,L110,L111,L113,L115,L116,L118,L119,L120,L121,L122,L125,L126,L127,L128,L129,L130,L131,L133,L134,L135,L136,L137,L141,L142,L143,L145,L146,L147,L148,L149,L150,L151,L152,L153,L154,L155,L157,L158,L159,L161,L162,L163,L165,L169,L170,L171,L174,L175,L176,L177,L178,L180,L181,L182,L183,L184,L186,L187,L188,L189,L190,L192,L193,L194,L195,L196,L197,L201,L202,L203,L204,L205,L206,L207,L208,L209,L210,L211,L213,L214,L215,L216,L217,L219,L221,L222,L223,L224,L225,L226,L228,L229,L230,L231,L232,L233,L234,L235,L237,L239,L240,L241,L242,L243,L244,L245,L246,L247,L248,L249,L250,L251,L252,L253,L255,L256,L257,L258,L260,L261,L262,L264,L265,L266,L267,L269,L270,L271,L272,L275,L276,L277,L279,L280,L281,L282,L288,L289,L290,L292,L293,L294,L295,L296,L297,L298,L301,L302,L303,L305,L306,L307,L310,L311,L312,L314,L315,L316,L317,L318,L319,L320,L321,L322,L325,L327,L328,L329,L330,L331,L333,L334,L335,L336,L337,L338,L339,L340,L341,L342,L343,L344,L346,L347,L348,L349)
as.matrix(allrates)->allrates
###sortby length
length(L216)->leL216
length(L182)->leL182
length(L149)->leL149
length(L213)->leL213
length(L184)->leL184
length(L223)->leL223
length(L305)->leL305
length(L8)->leL8
length(L203)->leL203
length(L95)->leL95
length(L41)->leL41
length(L295)->leL295
length(L107)->leL107
length(L163)->leL163
length(L21)->leL21
length(L2)->leL2
length(L325)->leL325
length(L80)->leL80
length(L28)->leL28

c(leL216,leL182,leL149,leL213,leL184,leL223,leL305,leL8,leL203,leL95,leL41,leL295,leL107,leL163,leL21,leL2,leL325,leL80,leL28)->ll
names(ll)<-c("leL216","leL182","leL149","leL213","leL184","leL223","leL305","leL8","leL203","leL95","leL41","leL295","leL107","leL163","leL21","leL2","leL325","leL80","leL28")
sort(ll)
###lets do the crown of neoaves @ 0.30
space.maker(allrates,.30,3)->a1
space.maker(L213,.30,3)->a2
space.maker(L182,.30,3)->a3
space.maker(L203,.30,3)->a4
space.maker(L149,.30,3)->a5
space.maker(L2,.30,3)->a6
space.maker(L295,.30,3)->a7
space.maker(L325,.30,3)->a8
space.maker(L107,.30,3)->a9
space.maker(L41,.30,3)->a10
space.maker(L216,.30,3)->a11
space.maker(L28,.30,3)->a12
space.maker(L184,.30,3)->a13
space.maker(L8,.30,3)->a14
space.maker(L95,.30,3)->a15
space.maker(L80,.30,3)->a16
space.maker(L163,.30,3)->a17
space.maker(L21,.30,3)->a18
space.maker(L305,.30,3)->a19
space.maker(L223,.30,3)->a20

rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)->demo
as.matrix(demo)->demo2
row.names(demo2)<-c("allrates","L213","leL182","leL203","leL149","leL2","leL295","leL325","leL107","leL41","leL216","leL28",
"leL184","leL8","leL95","leL80","leL163","leL21","leL305","leL223")

.30/20->by.this
seq(by.this,0.30-0.0001,by=by.this)->lilts

colnames(demo2)<-lilts
heatmap.2(demo2, Colv=F,Rowv=F, scale='none')