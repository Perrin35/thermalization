OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.78491640090942) q[0];
sx q[0];
rz(1.69503215153749) q[0];
sx q[0];
rz(8.52300814389392) q[0];
rz(4.47373676300049) q[1];
sx q[1];
rz(2.69822150667245) q[1];
sx q[1];
rz(8.65952095984622) q[1];
cx q[1],q[0];
rz(-0.0888521000742912) q[0];
sx q[0];
rz(5.60616269906098) q[0];
sx q[0];
rz(11.9637343645017) q[0];
rz(-1.32773625850677) q[2];
sx q[2];
rz(4.38366607029969) q[2];
sx q[2];
rz(13.0619899988095) q[2];
cx q[2],q[1];
rz(2.26970767974854) q[1];
sx q[1];
rz(2.80083400209481) q[1];
sx q[1];
rz(11.1158413648526) q[1];
rz(0.428496897220612) q[3];
sx q[3];
rz(3.44660729368264) q[3];
sx q[3];
rz(11.1712916850965) q[3];
cx q[3],q[2];
rz(0.296099334955215) q[2];
sx q[2];
rz(5.19503608544404) q[2];
sx q[2];
rz(10.8968305349271) q[2];
rz(1.71464014053345) q[3];
sx q[3];
rz(4.78574827511842) q[3];
sx q[3];
rz(8.89105401038333) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.18328368663788) q[0];
sx q[0];
rz(4.7754135449701) q[0];
sx q[0];
rz(11.7175931692044) q[0];
rz(-0.0350041836500168) q[1];
sx q[1];
rz(5.13634720643098) q[1];
sx q[1];
rz(8.47120598553821) q[1];
cx q[1],q[0];
rz(2.78244209289551) q[0];
sx q[0];
rz(1.11385646660859) q[0];
sx q[0];
rz(10.3600345611493) q[0];
rz(-2.4875819683075) q[2];
sx q[2];
rz(3.61870446999604) q[2];
sx q[2];
rz(9.5185110181491) q[2];
cx q[2],q[1];
rz(-1.19989049434662) q[1];
sx q[1];
rz(-0.395597068471364) q[1];
sx q[1];
rz(12.0381920099179) q[1];
rz(0.816166579723358) q[3];
sx q[3];
rz(1.07936945756013) q[3];
sx q[3];
rz(8.0243129491727) q[3];
cx q[3],q[2];
rz(1.16350519657135) q[2];
sx q[2];
rz(4.92299512227113) q[2];
sx q[2];
rz(11.7168221235196) q[2];
rz(0.165440112352371) q[3];
sx q[3];
rz(1.44355073769624) q[3];
sx q[3];
rz(11.8166022062223) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.77168321609497) q[0];
sx q[0];
rz(2.21943798859651) q[0];
sx q[0];
rz(12.9682936429898) q[0];
rz(-0.153321281075478) q[1];
sx q[1];
rz(6.46005025704438) q[1];
sx q[1];
rz(10.5609707593839) q[1];
cx q[1],q[0];
rz(0.763218522071838) q[0];
sx q[0];
rz(-0.399706689519338) q[0];
sx q[0];
rz(9.45565907693609) q[0];
rz(2.86747455596924) q[2];
sx q[2];
rz(5.18042865593965) q[2];
sx q[2];
rz(10.0016815423886) q[2];
cx q[2],q[1];
rz(0.786277770996094) q[1];
sx q[1];
rz(4.48002913792665) q[1];
sx q[1];
rz(6.03654882907077) q[1];
rz(1.10688126087189) q[3];
sx q[3];
rz(1.21642103989656) q[3];
sx q[3];
rz(9.84062710999652) q[3];
cx q[3],q[2];
rz(-0.0841269642114639) q[2];
sx q[2];
rz(4.02015403111512) q[2];
sx q[2];
rz(9.4798346221368) q[2];
rz(1.34752404689789) q[3];
sx q[3];
rz(4.95305064518983) q[3];
sx q[3];
rz(8.42048881053134) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.364430010318756) q[0];
sx q[0];
rz(2.93685201008851) q[0];
sx q[0];
rz(7.85077641009494) q[0];
rz(-0.691529989242554) q[1];
sx q[1];
rz(5.61389222939546) q[1];
sx q[1];
rz(9.23935209809943) q[1];
cx q[1],q[0];
rz(1.79509925842285) q[0];
sx q[0];
rz(-0.43432625929778) q[0];
sx q[0];
rz(7.44745252131625) q[0];
rz(0.859472095966339) q[2];
sx q[2];
rz(2.53441158135469) q[2];
sx q[2];
rz(11.39274761676) q[2];
cx q[2],q[1];
rz(2.59258270263672) q[1];
sx q[1];
rz(1.82357779343659) q[1];
sx q[1];
rz(9.74624899624988) q[1];
rz(-1.02493989467621) q[3];
sx q[3];
rz(4.53510835965211) q[3];
sx q[3];
rz(10.7703457832257) q[3];
cx q[3],q[2];
rz(1.73762476444244) q[2];
sx q[2];
rz(4.57058850129182) q[2];
sx q[2];
rz(9.41954755959242) q[2];
rz(0.387492030858994) q[3];
sx q[3];
rz(3.64627775748307) q[3];
sx q[3];
rz(10.7277828216474) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.72652792930603) q[0];
sx q[0];
rz(4.12512758572633) q[0];
sx q[0];
rz(10.0407521486203) q[0];
rz(5.09432506561279) q[1];
sx q[1];
rz(5.66483870347077) q[1];
sx q[1];
rz(9.9380856513898) q[1];
cx q[1],q[0];
rz(0.0307260416448116) q[0];
sx q[0];
rz(6.56488290627534) q[0];
sx q[0];
rz(6.25659630297824) q[0];
rz(5.43816518783569) q[2];
sx q[2];
rz(1.35755244095857) q[2];
sx q[2];
rz(14.8010320424955) q[2];
cx q[2],q[1];
rz(-1.75558805465698) q[1];
sx q[1];
rz(4.64034047921235) q[1];
sx q[1];
rz(10.3603853940885) q[1];
rz(2.67424178123474) q[3];
sx q[3];
rz(7.9597374518686) q[3];
sx q[3];
rz(9.08504167794391) q[3];
cx q[3],q[2];
rz(1.00802278518677) q[2];
sx q[2];
rz(1.27613893349702) q[2];
sx q[2];
rz(8.68230066298648) q[2];
rz(1.69784688949585) q[3];
sx q[3];
rz(4.87390402157838) q[3];
sx q[3];
rz(7.38452169894382) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.831085741519928) q[0];
sx q[0];
rz(2.09489122231538) q[0];
sx q[0];
rz(9.06496933697864) q[0];
rz(0.85044914484024) q[1];
sx q[1];
rz(4.32281139691407) q[1];
sx q[1];
rz(6.74901602267429) q[1];
cx q[1],q[0];
rz(-0.820792615413666) q[0];
sx q[0];
rz(2.38938364584977) q[0];
sx q[0];
rz(9.74397478102847) q[0];
rz(-1.4288387298584) q[2];
sx q[2];
rz(3.90989420016343) q[2];
sx q[2];
rz(12.463709807388) q[2];
cx q[2],q[1];
rz(1.17594313621521) q[1];
sx q[1];
rz(4.66138699849183) q[1];
sx q[1];
rz(7.61292538642093) q[1];
rz(-2.39218878746033) q[3];
sx q[3];
rz(2.38300493557984) q[3];
sx q[3];
rz(11.8489770650785) q[3];
cx q[3],q[2];
rz(0.742385268211365) q[2];
sx q[2];
rz(1.96116045315797) q[2];
sx q[2];
rz(12.9307243585508) q[2];
rz(-0.243891671299934) q[3];
sx q[3];
rz(6.23359313805635) q[3];
sx q[3];
rz(11.3490918636243) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.57049822807312) q[0];
sx q[0];
rz(2.16223928530748) q[0];
sx q[0];
rz(11.4003479242246) q[0];
rz(2.12655377388) q[1];
sx q[1];
rz(5.74615112145478) q[1];
sx q[1];
rz(12.1799194574277) q[1];
cx q[1],q[0];
rz(3.72341918945313) q[0];
sx q[0];
rz(6.4950579722696) q[0];
sx q[0];
rz(9.27853407560989) q[0];
rz(-0.256623029708862) q[2];
sx q[2];
rz(4.41176739533479) q[2];
sx q[2];
rz(13.8929323911588) q[2];
cx q[2],q[1];
rz(3.84289002418518) q[1];
sx q[1];
rz(6.36313548882539) q[1];
sx q[1];
rz(12.8909549474637) q[1];
rz(3.68389320373535) q[3];
sx q[3];
rz(1.46315482457215) q[3];
sx q[3];
rz(9.44309572166904) q[3];
cx q[3],q[2];
rz(-0.405479073524475) q[2];
sx q[2];
rz(1.49612859089906) q[2];
sx q[2];
rz(12.3386716604154) q[2];
rz(-2.06857109069824) q[3];
sx q[3];
rz(5.53438988526399) q[3];
sx q[3];
rz(9.13120681642696) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.881840825080872) q[0];
sx q[0];
rz(4.04143193562562) q[0];
sx q[0];
rz(9.85279068946048) q[0];
rz(4.15987873077393) q[1];
sx q[1];
rz(3.57794004877145) q[1];
sx q[1];
rz(10.2958087682645) q[1];
cx q[1],q[0];
rz(-0.508650958538055) q[0];
sx q[0];
rz(2.06744602520997) q[0];
sx q[0];
rz(13.214098429672) q[0];
rz(3.49892354011536) q[2];
sx q[2];
rz(7.66010490258271) q[2];
sx q[2];
rz(9.27844903468295) q[2];
cx q[2],q[1];
rz(0.582951366901398) q[1];
sx q[1];
rz(4.59145227273042) q[1];
sx q[1];
rz(12.7084572076718) q[1];
rz(-0.915240406990051) q[3];
sx q[3];
rz(3.9375537951761) q[3];
sx q[3];
rz(11.3812784910123) q[3];
cx q[3],q[2];
rz(-3.27387928962708) q[2];
sx q[2];
rz(5.48761049111421) q[2];
sx q[2];
rz(13.8450383901517) q[2];
rz(3.72006678581238) q[3];
sx q[3];
rz(5.65067282517488) q[3];
sx q[3];
rz(11.4123676776807) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.46716749668121) q[0];
sx q[0];
rz(3.35867965419824) q[0];
sx q[0];
rz(11.9665632009427) q[0];
rz(1.8425874710083) q[1];
sx q[1];
rz(4.75195113022859) q[1];
sx q[1];
rz(5.64131734370395) q[1];
cx q[1],q[0];
rz(3.22109246253967) q[0];
sx q[0];
rz(2.53895172675187) q[0];
sx q[0];
rz(11.6261904001157) q[0];
rz(3.90688681602478) q[2];
sx q[2];
rz(1.46647158463533) q[2];
sx q[2];
rz(12.8217627763669) q[2];
cx q[2],q[1];
rz(-0.363433450460434) q[1];
sx q[1];
rz(1.55234077771241) q[1];
sx q[1];
rz(8.56177721022769) q[1];
rz(-0.526858866214752) q[3];
sx q[3];
rz(5.39682117302949) q[3];
sx q[3];
rz(8.14524970053836) q[3];
cx q[3],q[2];
rz(0.820152640342712) q[2];
sx q[2];
rz(1.71653107007081) q[2];
sx q[2];
rz(13.1766326188962) q[2];
rz(-0.520147025585175) q[3];
sx q[3];
rz(4.19616201718385) q[3];
sx q[3];
rz(8.93125826715633) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.80101203918457) q[0];
sx q[0];
rz(4.38511792023713) q[0];
sx q[0];
rz(10.3627946734349) q[0];
rz(0.341713786125183) q[1];
sx q[1];
rz(4.52370968659455) q[1];
sx q[1];
rz(9.58206327854797) q[1];
cx q[1],q[0];
rz(-0.926259160041809) q[0];
sx q[0];
rz(3.79660263855989) q[0];
sx q[0];
rz(10.6732671022336) q[0];
rz(-2.32882308959961) q[2];
sx q[2];
rz(4.77384391625459) q[2];
sx q[2];
rz(9.36628453656241) q[2];
cx q[2],q[1];
rz(-2.59636354446411) q[1];
sx q[1];
rz(5.66139808495576) q[1];
sx q[1];
rz(10.4607180118482) q[1];
rz(0.408779889345169) q[3];
sx q[3];
rz(4.47803130944306) q[3];
sx q[3];
rz(11.9294447660367) q[3];
cx q[3],q[2];
rz(-3.36050581932068) q[2];
sx q[2];
rz(5.93813124497468) q[2];
sx q[2];
rz(10.6986576080243) q[2];
rz(0.00565654691308737) q[3];
sx q[3];
rz(4.8417171557718) q[3];
sx q[3];
rz(7.26421473025485) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.4411678314209) q[0];
sx q[0];
rz(1.52529743512208) q[0];
sx q[0];
rz(6.68367216586276) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.87144160270691) q[1];
sx q[1];
rz(8.47952190239961) q[1];
sx q[1];
rz(9.70323169826671) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(2.01329016685486) q[2];
sx q[2];
rz(0.873947056131907) q[2];
sx q[2];
rz(5.44734117983981) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-2.70399832725525) q[3];
sx q[3];
rz(4.63481930096681) q[3];
sx q[3];
rz(12.491354918472) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
