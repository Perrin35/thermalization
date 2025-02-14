OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.651694595813751) q[0];
sx q[0];
rz(5.40138474305207) q[0];
sx q[0];
rz(9.28009792267486) q[0];
rz(3.55968499183655) q[1];
sx q[1];
rz(3.24536446680362) q[1];
sx q[1];
rz(11.0544626474301) q[1];
cx q[1],q[0];
rz(-0.666833162307739) q[0];
sx q[0];
rz(1.07682433922822) q[0];
sx q[0];
rz(7.98846802710696) q[0];
rz(-1.51060771942139) q[2];
sx q[2];
rz(7.01758304436738) q[2];
sx q[2];
rz(13.3028268575589) q[2];
cx q[2],q[1];
rz(-0.50070333480835) q[1];
sx q[1];
rz(1.90240600903565) q[1];
sx q[1];
rz(9.01476267575427) q[1];
rz(0.898999631404877) q[3];
sx q[3];
rz(6.90590921242768) q[3];
sx q[3];
rz(10.5485343694608) q[3];
cx q[3],q[2];
rz(3.45802998542786) q[2];
sx q[2];
rz(1.87081590493257) q[2];
sx q[2];
rz(10.0801197647969) q[2];
rz(1.75730276107788) q[3];
sx q[3];
rz(4.10610053141648) q[3];
sx q[3];
rz(8.6038172006528) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.692193865776062) q[0];
sx q[0];
rz(1.37575021584565) q[0];
sx q[0];
rz(9.93291912078067) q[0];
rz(-1.38057684898376) q[1];
sx q[1];
rz(5.70180550416047) q[1];
sx q[1];
rz(8.21520588397189) q[1];
cx q[1],q[0];
rz(-0.660994529724121) q[0];
sx q[0];
rz(3.75995508034761) q[0];
sx q[0];
rz(9.10504672526523) q[0];
rz(0.0176969282329082) q[2];
sx q[2];
rz(4.76597669919068) q[2];
sx q[2];
rz(5.82770869731113) q[2];
cx q[2],q[1];
rz(0.758250713348389) q[1];
sx q[1];
rz(2.03133252461488) q[1];
sx q[1];
rz(9.56762616931602) q[1];
rz(0.534044921398163) q[3];
sx q[3];
rz(8.74885669549043) q[3];
sx q[3];
rz(9.01830939053699) q[3];
cx q[3],q[2];
rz(-1.93805825710297) q[2];
sx q[2];
rz(4.53761568863923) q[2];
sx q[2];
rz(6.90008375643894) q[2];
rz(-1.53617024421692) q[3];
sx q[3];
rz(4.18574574788148) q[3];
sx q[3];
rz(9.91561198829814) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.15983748435974) q[0];
sx q[0];
rz(1.4977389891916) q[0];
sx q[0];
rz(8.69971141814395) q[0];
rz(-1.26953506469727) q[1];
sx q[1];
rz(3.81924912531907) q[1];
sx q[1];
rz(10.3839533686559) q[1];
cx q[1],q[0];
rz(5.27801847457886) q[0];
sx q[0];
rz(1.10148969491059) q[0];
sx q[0];
rz(10.5584176540296) q[0];
rz(4.2964825630188) q[2];
sx q[2];
rz(4.35403695900971) q[2];
sx q[2];
rz(12.791115975372) q[2];
cx q[2],q[1];
rz(0.196344643831253) q[1];
sx q[1];
rz(4.52913192112977) q[1];
sx q[1];
rz(10.879039502136) q[1];
rz(-1.40149164199829) q[3];
sx q[3];
rz(4.89219978650148) q[3];
sx q[3];
rz(10.81952891349) q[3];
cx q[3],q[2];
rz(1.8150669336319) q[2];
sx q[2];
rz(4.40837159951264) q[2];
sx q[2];
rz(9.16777775286838) q[2];
rz(1.08109307289124) q[3];
sx q[3];
rz(0.52101460297639) q[3];
sx q[3];
rz(7.84602258204623) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.08398878574371) q[0];
sx q[0];
rz(2.86315176089341) q[0];
sx q[0];
rz(8.24688682555362) q[0];
rz(-0.150116935372353) q[1];
sx q[1];
rz(3.67373672326142) q[1];
sx q[1];
rz(10.9410905599515) q[1];
cx q[1],q[0];
rz(2.02228426933289) q[0];
sx q[0];
rz(4.43065026600892) q[0];
sx q[0];
rz(9.08019137977763) q[0];
rz(-0.728677213191986) q[2];
sx q[2];
rz(4.42999485333497) q[2];
sx q[2];
rz(8.90803173779651) q[2];
cx q[2],q[1];
rz(3.14322972297668) q[1];
sx q[1];
rz(1.47125843365724) q[1];
sx q[1];
rz(8.18083891867801) q[1];
rz(0.189900934696198) q[3];
sx q[3];
rz(4.07737294037873) q[3];
sx q[3];
rz(9.03249779938861) q[3];
cx q[3],q[2];
rz(-0.842450439929962) q[2];
sx q[2];
rz(3.71148112614686) q[2];
sx q[2];
rz(11.6432635545652) q[2];
rz(-0.959215700626373) q[3];
sx q[3];
rz(4.15425637562806) q[3];
sx q[3];
rz(9.20853633283778) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.320903092622757) q[0];
sx q[0];
rz(4.0454262812906) q[0];
sx q[0];
rz(13.4492816686551) q[0];
rz(1.29549050331116) q[1];
sx q[1];
rz(0.774064453440257) q[1];
sx q[1];
rz(8.92980707286998) q[1];
cx q[1],q[0];
rz(0.519091010093689) q[0];
sx q[0];
rz(1.25357309182221) q[0];
sx q[0];
rz(11.3003659009854) q[0];
rz(-0.963168084621429) q[2];
sx q[2];
rz(3.94286677439744) q[2];
sx q[2];
rz(10.0281117916028) q[2];
cx q[2],q[1];
rz(2.04308295249939) q[1];
sx q[1];
rz(7.94175496895845) q[1];
sx q[1];
rz(10.3261590361516) q[1];
rz(-0.576822459697723) q[3];
sx q[3];
rz(-0.601895419759206) q[3];
sx q[3];
rz(9.68921092747852) q[3];
cx q[3],q[2];
rz(4.95145845413208) q[2];
sx q[2];
rz(3.26114689012105) q[2];
sx q[2];
rz(8.1565960407178) q[2];
rz(-0.0825033336877823) q[3];
sx q[3];
rz(1.63763299782807) q[3];
sx q[3];
rz(10.0818099141042) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.65620112419128) q[0];
sx q[0];
rz(3.34842923481996) q[0];
sx q[0];
rz(7.92326793669864) q[0];
rz(-5.22071027755737) q[1];
sx q[1];
rz(1.1524008830362) q[1];
sx q[1];
rz(10.6001267194669) q[1];
cx q[1],q[0];
rz(2.99060153961182) q[0];
sx q[0];
rz(6.66272464592988) q[0];
sx q[0];
rz(7.22373459338352) q[0];
rz(-2.84789991378784) q[2];
sx q[2];
rz(-0.898646203679494) q[2];
sx q[2];
rz(11.3791022062223) q[2];
cx q[2],q[1];
rz(-2.29860019683838) q[1];
sx q[1];
rz(4.82958105404908) q[1];
sx q[1];
rz(11.2489305496137) q[1];
rz(-0.810445249080658) q[3];
sx q[3];
rz(4.93107250531251) q[3];
sx q[3];
rz(6.04547975062534) q[3];
cx q[3],q[2];
rz(-0.330927342176437) q[2];
sx q[2];
rz(1.28454795678193) q[2];
sx q[2];
rz(11.5786878824155) q[2];
rz(-2.2641613483429) q[3];
sx q[3];
rz(2.87434470851953) q[3];
sx q[3];
rz(13.2964389085691) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.18386089801788) q[0];
sx q[0];
rz(1.53936246235902) q[0];
sx q[0];
rz(12.4275111913602) q[0];
rz(-1.21443700790405) q[1];
sx q[1];
rz(4.64934590657289) q[1];
sx q[1];
rz(8.60818514823123) q[1];
cx q[1],q[0];
rz(4.02745580673218) q[0];
sx q[0];
rz(9.65039172967012) q[0];
sx q[0];
rz(9.33950510471269) q[0];
rz(3.94414520263672) q[2];
sx q[2];
rz(4.90159133275086) q[2];
sx q[2];
rz(9.60652121006652) q[2];
cx q[2],q[1];
rz(1.54677557945251) q[1];
sx q[1];
rz(5.56344405015046) q[1];
sx q[1];
rz(9.48163967057272) q[1];
rz(0.904666364192963) q[3];
sx q[3];
rz(2.17628476222093) q[3];
sx q[3];
rz(8.98358387350246) q[3];
cx q[3],q[2];
rz(0.821780860424042) q[2];
sx q[2];
rz(0.374151619272777) q[2];
sx q[2];
rz(8.94011936186954) q[2];
rz(-0.365684926509857) q[3];
sx q[3];
rz(4.50607100327546) q[3];
sx q[3];
rz(11.4075617551725) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.13002383708954) q[0];
sx q[0];
rz(9.11365619500215) q[0];
sx q[0];
rz(9.52200960218116) q[0];
rz(0.104879960417747) q[1];
sx q[1];
rz(3.92129132350022) q[1];
sx q[1];
rz(10.7394555568616) q[1];
cx q[1],q[0];
rz(2.07856941223145) q[0];
sx q[0];
rz(3.13483506267006) q[0];
sx q[0];
rz(9.78335700034305) q[0];
rz(-2.90005230903625) q[2];
sx q[2];
rz(4.92976227601106) q[2];
sx q[2];
rz(7.83067963122531) q[2];
cx q[2],q[1];
rz(4.56316041946411) q[1];
sx q[1];
rz(2.08892575104768) q[1];
sx q[1];
rz(6.22253391741916) q[1];
rz(1.48363626003265) q[3];
sx q[3];
rz(1.03565517266328) q[3];
sx q[3];
rz(11.9700886964719) q[3];
cx q[3],q[2];
rz(0.176819920539856) q[2];
sx q[2];
rz(4.53918078740174) q[2];
sx q[2];
rz(12.4131045102994) q[2];
rz(1.2166473865509) q[3];
sx q[3];
rz(6.06511941750581) q[3];
sx q[3];
rz(10.0409638643186) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.0230872631073) q[0];
sx q[0];
rz(3.14707195119048) q[0];
sx q[0];
rz(10.9294785022657) q[0];
rz(-0.798323750495911) q[1];
sx q[1];
rz(1.61838045914704) q[1];
sx q[1];
rz(9.76009274124309) q[1];
cx q[1],q[0];
rz(0.178517952561378) q[0];
sx q[0];
rz(4.97346106370027) q[0];
sx q[0];
rz(10.9302229642789) q[0];
rz(1.54456090927124) q[2];
sx q[2];
rz(3.89215693076188) q[2];
sx q[2];
rz(7.26613495349094) q[2];
cx q[2],q[1];
rz(-0.112869523465633) q[1];
sx q[1];
rz(4.99929598172242) q[1];
sx q[1];
rz(11.3508388757627) q[1];
rz(1.7280900478363) q[3];
sx q[3];
rz(5.28406635125215) q[3];
sx q[3];
rz(7.52729580401584) q[3];
cx q[3],q[2];
rz(3.01407837867737) q[2];
sx q[2];
rz(4.31002691586549) q[2];
sx q[2];
rz(12.1287774801175) q[2];
rz(-2.09949994087219) q[3];
sx q[3];
rz(4.54691747029359) q[3];
sx q[3];
rz(13.1088287591855) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.15227782726288) q[0];
sx q[0];
rz(2.25700256426866) q[0];
sx q[0];
rz(13.0159837961118) q[0];
rz(2.16470289230347) q[1];
sx q[1];
rz(4.77927437623078) q[1];
sx q[1];
rz(9.92599908112689) q[1];
cx q[1],q[0];
rz(1.65274429321289) q[0];
sx q[0];
rz(2.79786589940126) q[0];
sx q[0];
rz(10.0196058511655) q[0];
rz(-0.451480001211166) q[2];
sx q[2];
rz(6.74806466897065) q[2];
sx q[2];
rz(7.63820407389804) q[2];
cx q[2],q[1];
rz(-3.52952980995178) q[1];
sx q[1];
rz(6.81051650841767) q[1];
sx q[1];
rz(13.0503375291745) q[1];
rz(3.7595386505127) q[3];
sx q[3];
rz(6.83496204217012) q[3];
sx q[3];
rz(8.48078939913913) q[3];
cx q[3],q[2];
rz(8.40126514434814) q[2];
sx q[2];
rz(4.03395781119401) q[2];
sx q[2];
rz(6.06546399592563) q[2];
rz(1.90675413608551) q[3];
sx q[3];
rz(3.80559149582917) q[3];
sx q[3];
rz(10.3454301714818) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.12130260467529) q[0];
sx q[0];
rz(1.19422844250733) q[0];
sx q[0];
rz(6.1303384065549) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(3.33788776397705) q[1];
sx q[1];
rz(2.79604822595651) q[1];
sx q[1];
rz(8.34358928202792) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.196227759122849) q[2];
sx q[2];
rz(1.42304304440553) q[2];
sx q[2];
rz(10.5460710287015) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(2.61182856559753) q[3];
sx q[3];
rz(3.40986693103845) q[3];
sx q[3];
rz(19.8024549245755) q[3];
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
