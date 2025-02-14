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
rz(-2.01055359840393) q[0];
sx q[0];
rz(3.717488976317) q[0];
sx q[0];
rz(5.15626380442783) q[0];
rz(2.19088864326477) q[1];
sx q[1];
rz(3.74284407694871) q[1];
sx q[1];
rz(9.75847015380069) q[1];
cx q[1],q[0];
rz(5.343337059021) q[0];
sx q[0];
rz(2.69320953090722) q[0];
sx q[0];
rz(9.77634356021091) q[0];
rz(2.68621802330017) q[2];
sx q[2];
rz(4.47135201294953) q[2];
sx q[2];
rz(9.49318194984599) q[2];
cx q[2],q[1];
rz(-0.0444943904876709) q[1];
sx q[1];
rz(1.35840288003022) q[1];
sx q[1];
rz(13.7863502263944) q[1];
rz(3.41435790061951) q[3];
sx q[3];
rz(2.35415098269517) q[3];
sx q[3];
rz(10.7637104749601) q[3];
cx q[3],q[2];
rz(-7.54328155517578) q[2];
sx q[2];
rz(8.34918084939057) q[2];
sx q[2];
rz(2.09137008189365) q[2];
rz(2.84607338905334) q[3];
sx q[3];
rz(7.05204239686067) q[3];
sx q[3];
rz(13.9356627225797) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.16939008235931) q[0];
sx q[0];
rz(1.01935950120027) q[0];
sx q[0];
rz(10.0920347332875) q[0];
rz(3.29229474067688) q[1];
sx q[1];
rz(4.19639423687989) q[1];
sx q[1];
rz(9.73642513751193) q[1];
cx q[1],q[0];
rz(1.22094762325287) q[0];
sx q[0];
rz(4.19157460530336) q[0];
sx q[0];
rz(14.0558852910916) q[0];
rz(3.73685240745544) q[2];
sx q[2];
rz(0.515805872278758) q[2];
sx q[2];
rz(11.7467262506406) q[2];
cx q[2],q[1];
rz(-4.43944120407104) q[1];
sx q[1];
rz(9.19357076485688) q[1];
sx q[1];
rz(8.50592432021304) q[1];
rz(3.06875419616699) q[3];
sx q[3];
rz(1.17411163647706) q[3];
sx q[3];
rz(8.2829479932706) q[3];
cx q[3],q[2];
rz(0.475058913230896) q[2];
sx q[2];
rz(4.03502401907975) q[2];
sx q[2];
rz(7.05082080363437) q[2];
rz(-3.62248635292053) q[3];
sx q[3];
rz(2.37002626259858) q[3];
sx q[3];
rz(7.87007806300327) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.84592205286026) q[0];
sx q[0];
rz(4.74557379086549) q[0];
sx q[0];
rz(10.2043604016225) q[0];
rz(-5.91145372390747) q[1];
sx q[1];
rz(4.14916101296479) q[1];
sx q[1];
rz(10.1857595801274) q[1];
cx q[1],q[0];
rz(-0.261393815279007) q[0];
sx q[0];
rz(6.72191801865632) q[0];
sx q[0];
rz(10.1608087181966) q[0];
rz(0.518634676933289) q[2];
sx q[2];
rz(5.5274716933542) q[2];
sx q[2];
rz(10.5664906263272) q[2];
cx q[2],q[1];
rz(1.07989883422852) q[1];
sx q[1];
rz(1.73622802098329) q[1];
sx q[1];
rz(2.41676995753452) q[1];
rz(-3.14725875854492) q[3];
sx q[3];
rz(3.85969636042649) q[3];
sx q[3];
rz(15.2606668233792) q[3];
cx q[3],q[2];
rz(-0.804282128810883) q[2];
sx q[2];
rz(0.90039840539033) q[2];
sx q[2];
rz(6.73835990428134) q[2];
rz(3.83353662490845) q[3];
sx q[3];
rz(2.42718360026414) q[3];
sx q[3];
rz(10.2335711479108) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.968750178813934) q[0];
sx q[0];
rz(4.48854878743226) q[0];
sx q[0];
rz(12.8226585149686) q[0];
rz(4.8484582901001) q[1];
sx q[1];
rz(8.0043422301584) q[1];
sx q[1];
rz(12.6851253271024) q[1];
cx q[1],q[0];
rz(-1.25601208209991) q[0];
sx q[0];
rz(2.66544541914994) q[0];
sx q[0];
rz(7.51008567809268) q[0];
rz(-1.44271492958069) q[2];
sx q[2];
rz(4.42150739033753) q[2];
sx q[2];
rz(5.16748616694614) q[2];
cx q[2],q[1];
rz(-3.2657642364502) q[1];
sx q[1];
rz(1.33914366562898) q[1];
sx q[1];
rz(10.2005277037542) q[1];
rz(-0.815779209136963) q[3];
sx q[3];
rz(4.99259975750978) q[3];
sx q[3];
rz(11.2356723308484) q[3];
cx q[3],q[2];
rz(0.287483006715775) q[2];
sx q[2];
rz(3.41526034672792) q[2];
sx q[2];
rz(7.33886859416171) q[2];
rz(1.69151973724365) q[3];
sx q[3];
rz(4.58886423905427) q[3];
sx q[3];
rz(13.4153892755429) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.57978224754333) q[0];
sx q[0];
rz(6.10523691971833) q[0];
sx q[0];
rz(7.81127414702579) q[0];
rz(-1.96445071697235) q[1];
sx q[1];
rz(4.41168585618074) q[1];
sx q[1];
rz(7.86156163214847) q[1];
cx q[1],q[0];
rz(1.62306833267212) q[0];
sx q[0];
rz(4.08777889807756) q[0];
sx q[0];
rz(9.51018854825898) q[0];
rz(-1.98149478435516) q[2];
sx q[2];
rz(5.01514533360536) q[2];
sx q[2];
rz(8.13561377524539) q[2];
cx q[2],q[1];
rz(-7.56124019622803) q[1];
sx q[1];
rz(-1.71249326865142) q[1];
sx q[1];
rz(8.48911116122409) q[1];
rz(4.70621919631958) q[3];
sx q[3];
rz(3.55100915034349) q[3];
sx q[3];
rz(5.31678102015659) q[3];
cx q[3],q[2];
rz(-0.135811805725098) q[2];
sx q[2];
rz(1.15671888192231) q[2];
sx q[2];
rz(10.7269453763883) q[2];
rz(-5.94323396682739) q[3];
sx q[3];
rz(5.476466806727) q[3];
sx q[3];
rz(10.0871665239255) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.00789167732000351) q[0];
sx q[0];
rz(4.73871293862397) q[0];
sx q[0];
rz(10.6245240926664) q[0];
rz(-1.78027355670929) q[1];
sx q[1];
rz(2.16845503647859) q[1];
sx q[1];
rz(11.9884865045468) q[1];
cx q[1],q[0];
rz(1.78521299362183) q[0];
sx q[0];
rz(4.72616532643373) q[0];
sx q[0];
rz(14.1307391881864) q[0];
rz(5.14522743225098) q[2];
sx q[2];
rz(1.71464005311067) q[2];
sx q[2];
rz(7.55401632785007) q[2];
cx q[2],q[1];
rz(-5.43402719497681) q[1];
sx q[1];
rz(-5.80739149252837) q[1];
sx q[1];
rz(11.295837378494) q[1];
rz(1.54505610466003) q[3];
sx q[3];
rz(-1.86960682074492) q[3];
sx q[3];
rz(10.0810429215352) q[3];
cx q[3],q[2];
rz(-0.631523847579956) q[2];
sx q[2];
rz(-0.574545947713307) q[2];
sx q[2];
rz(19.3230104207914) q[2];
rz(1.27246415615082) q[3];
sx q[3];
rz(4.89055660565431) q[3];
sx q[3];
rz(9.20232438146278) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.94721394777298) q[0];
sx q[0];
rz(6.94094768364961) q[0];
sx q[0];
rz(9.75791037677928) q[0];
rz(-3.37014007568359) q[1];
sx q[1];
rz(4.48177042801911) q[1];
sx q[1];
rz(5.70752618312045) q[1];
cx q[1],q[0];
rz(0.270283371210098) q[0];
sx q[0];
rz(5.51881256897981) q[0];
sx q[0];
rz(15.8612885236661) q[0];
rz(0.982515871524811) q[2];
sx q[2];
rz(0.555983932810374) q[2];
sx q[2];
rz(8.88291249274417) q[2];
cx q[2],q[1];
rz(3.68811631202698) q[1];
sx q[1];
rz(5.03822651703889) q[1];
sx q[1];
rz(4.92343566416904) q[1];
rz(-1.08955013751984) q[3];
sx q[3];
rz(4.93988076050813) q[3];
sx q[3];
rz(12.7216810941617) q[3];
cx q[3],q[2];
rz(2.87547421455383) q[2];
sx q[2];
rz(5.73175779183442) q[2];
sx q[2];
rz(7.9486441373746) q[2];
rz(-1.10097813606262) q[3];
sx q[3];
rz(8.36153999169404) q[3];
sx q[3];
rz(9.94277886151477) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.54408931732178) q[0];
sx q[0];
rz(4.02792981465394) q[0];
sx q[0];
rz(6.5903062581937) q[0];
rz(0.330445259809494) q[1];
sx q[1];
rz(4.68529871304566) q[1];
sx q[1];
rz(14.3428272962491) q[1];
cx q[1],q[0];
rz(0.489560514688492) q[0];
sx q[0];
rz(0.083261640863963) q[0];
sx q[0];
rz(7.63750717639133) q[0];
rz(1.22273254394531) q[2];
sx q[2];
rz(4.53530731995637) q[2];
sx q[2];
rz(12.9754316568296) q[2];
cx q[2],q[1];
rz(2.1077892780304) q[1];
sx q[1];
rz(6.79188433487947) q[1];
sx q[1];
rz(6.83539292811557) q[1];
rz(1.71242618560791) q[3];
sx q[3];
rz(-0.838115779561452) q[3];
sx q[3];
rz(10.8576406001966) q[3];
cx q[3],q[2];
rz(9.70030784606934) q[2];
sx q[2];
rz(-0.587746469182424) q[2];
sx q[2];
rz(9.67347667216464) q[2];
rz(-1.21348774433136) q[3];
sx q[3];
rz(4.58601191838319) q[3];
sx q[3];
rz(9.4864117115657) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.989517748355865) q[0];
sx q[0];
rz(2.23534205754335) q[0];
sx q[0];
rz(10.11137596368) q[0];
rz(7.39612674713135) q[1];
sx q[1];
rz(2.33908918698365) q[1];
sx q[1];
rz(5.96755311488315) q[1];
cx q[1],q[0];
rz(3.5154709815979) q[0];
sx q[0];
rz(10.9495603164011) q[0];
sx q[0];
rz(9.10164628028079) q[0];
rz(-0.248387932777405) q[2];
sx q[2];
rz(5.45666423638398) q[2];
sx q[2];
rz(12.9704816102903) q[2];
cx q[2],q[1];
rz(-1.31798470020294) q[1];
sx q[1];
rz(5.61416688759858) q[1];
sx q[1];
rz(3.29375550746127) q[1];
rz(1.64884519577026) q[3];
sx q[3];
rz(2.40253177483613) q[3];
sx q[3];
rz(8.54261026381656) q[3];
cx q[3],q[2];
rz(0.453720659017563) q[2];
sx q[2];
rz(6.84998241265351) q[2];
sx q[2];
rz(8.74020675419971) q[2];
rz(1.94139301776886) q[3];
sx q[3];
rz(2.01224163373048) q[3];
sx q[3];
rz(6.462758755676) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(4.2028489112854) q[0];
sx q[0];
rz(3.62273287971551) q[0];
sx q[0];
rz(6.28801867961093) q[0];
rz(-1.62394106388092) q[1];
sx q[1];
rz(7.02712622483308) q[1];
sx q[1];
rz(6.02947423457309) q[1];
cx q[1],q[0];
rz(4.08769083023071) q[0];
sx q[0];
rz(4.79774061043794) q[0];
sx q[0];
rz(12.7775633096616) q[0];
rz(-2.08480548858643) q[2];
sx q[2];
rz(0.804826410608836) q[2];
sx q[2];
rz(9.88064000605747) q[2];
cx q[2],q[1];
rz(6.2898645401001) q[1];
sx q[1];
rz(-0.906139222783498) q[1];
sx q[1];
rz(13.8877482175748) q[1];
rz(1.02030718326569) q[3];
sx q[3];
rz(1.74096778233583) q[3];
sx q[3];
rz(8.97733867763683) q[3];
cx q[3],q[2];
rz(0.978105127811432) q[2];
sx q[2];
rz(5.04818907578523) q[2];
sx q[2];
rz(4.32394168376132) q[2];
rz(0.676650822162628) q[3];
sx q[3];
rz(3.52815321286256) q[3];
sx q[3];
rz(10.438579773895) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-4.30792570114136) q[0];
sx q[0];
rz(4.50689152081544) q[0];
sx q[0];
rz(13.3730580568235) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.649163901805878) q[1];
sx q[1];
rz(-0.876666394872121) q[1];
sx q[1];
rz(9.62226240932151) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.5388223528862) q[2];
sx q[2];
rz(4.50139072735841) q[2];
sx q[2];
rz(9.73559380172893) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(2.02228260040283) q[3];
sx q[3];
rz(2.71970856388146) q[3];
sx q[3];
rz(7.34359142779514) q[3];
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
