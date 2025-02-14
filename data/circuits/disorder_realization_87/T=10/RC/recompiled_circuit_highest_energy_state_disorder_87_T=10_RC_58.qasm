OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6753766) q[0];
sx q[0];
rz(-2.197062) q[0];
sx q[0];
rz(-0.20139995) q[0];
rz(-2.0693076) q[1];
sx q[1];
rz(-2.3787002) q[1];
sx q[1];
rz(1.720517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.759255) q[0];
sx q[0];
rz(-1.5159642) q[0];
sx q[0];
rz(-2.1817529) q[0];
rz(-0.51789639) q[2];
sx q[2];
rz(-2.1764048) q[2];
sx q[2];
rz(-2.6014858) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0396467) q[1];
sx q[1];
rz(-1.1803487) q[1];
sx q[1];
rz(0.70266188) q[1];
rz(-pi) q[2];
rz(-2.8086923) q[3];
sx q[3];
rz(-1.9431291) q[3];
sx q[3];
rz(-1.2941192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1066771) q[2];
sx q[2];
rz(-2.3372529) q[2];
sx q[2];
rz(-2.8625281) q[2];
rz(-1.8883102) q[3];
sx q[3];
rz(-2.8918355) q[3];
sx q[3];
rz(2.9899924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6107553) q[0];
sx q[0];
rz(-2.294367) q[0];
sx q[0];
rz(0.24638677) q[0];
rz(-1.1950182) q[1];
sx q[1];
rz(-0.89391005) q[1];
sx q[1];
rz(-1.6349207) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5824597) q[0];
sx q[0];
rz(-1.3252859) q[0];
sx q[0];
rz(-1.1051635) q[0];
x q[1];
rz(1.7072524) q[2];
sx q[2];
rz(-2.8359957) q[2];
sx q[2];
rz(1.157925) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40827258) q[1];
sx q[1];
rz(-2.773914) q[1];
sx q[1];
rz(1.7321083) q[1];
rz(-0.43466062) q[3];
sx q[3];
rz(-2.0512329) q[3];
sx q[3];
rz(-0.91033376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40111497) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(1.9096036) q[2];
rz(1.1022107) q[3];
sx q[3];
rz(-3.1372742) q[3];
sx q[3];
rz(-1.8050885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6983637) q[0];
sx q[0];
rz(-0.6969499) q[0];
sx q[0];
rz(0.90782905) q[0];
rz(0.56697956) q[1];
sx q[1];
rz(-2.1588529) q[1];
sx q[1];
rz(-0.10279113) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1452652) q[0];
sx q[0];
rz(-1.780637) q[0];
sx q[0];
rz(2.5681433) q[0];
rz(-pi) q[1];
rz(3.0605761) q[2];
sx q[2];
rz(-1.4668494) q[2];
sx q[2];
rz(0.31652094) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72702209) q[1];
sx q[1];
rz(-2.1466781) q[1];
sx q[1];
rz(-2.058279) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0166753) q[3];
sx q[3];
rz(-1.301487) q[3];
sx q[3];
rz(0.96162063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28321442) q[2];
sx q[2];
rz(-2.3329222) q[2];
sx q[2];
rz(1.4374479) q[2];
rz(2.7490859) q[3];
sx q[3];
rz(-1.1350574) q[3];
sx q[3];
rz(0.2984305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0482386) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(1.1887953) q[0];
rz(-2.8554754) q[1];
sx q[1];
rz(-0.61758271) q[1];
sx q[1];
rz(1.6292705) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28933576) q[0];
sx q[0];
rz(-1.2931543) q[0];
sx q[0];
rz(-2.1136978) q[0];
rz(-0.34660201) q[2];
sx q[2];
rz(-0.79814974) q[2];
sx q[2];
rz(2.5426898) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0551853) q[1];
sx q[1];
rz(-2.1595104) q[1];
sx q[1];
rz(2.9179395) q[1];
rz(2.9305601) q[3];
sx q[3];
rz(-1.0730626) q[3];
sx q[3];
rz(-0.56609234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7901223) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(-2.4646087) q[2];
rz(-1.2466768) q[3];
sx q[3];
rz(-1.8691984) q[3];
sx q[3];
rz(-1.5164794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077654) q[0];
sx q[0];
rz(-0.25660577) q[0];
sx q[0];
rz(-3.1166792) q[0];
rz(-1.0554396) q[1];
sx q[1];
rz(-2.1565304) q[1];
sx q[1];
rz(-0.28344646) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6993857) q[0];
sx q[0];
rz(-1.3106579) q[0];
sx q[0];
rz(-0.26630638) q[0];
rz(-pi) q[1];
rz(1.3749429) q[2];
sx q[2];
rz(-1.699866) q[2];
sx q[2];
rz(1.376525) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.746096) q[1];
sx q[1];
rz(-1.0701792) q[1];
sx q[1];
rz(-3.1067645) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.153454) q[3];
sx q[3];
rz(-0.71708369) q[3];
sx q[3];
rz(-2.5337608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8289566) q[2];
sx q[2];
rz(-0.41746155) q[2];
sx q[2];
rz(0.32583315) q[2];
rz(-1.8587941) q[3];
sx q[3];
rz(-1.157434) q[3];
sx q[3];
rz(1.5475984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40389898) q[0];
sx q[0];
rz(-1.4827381) q[0];
sx q[0];
rz(2.7666336) q[0];
rz(-2.6629579) q[1];
sx q[1];
rz(-0.67239434) q[1];
sx q[1];
rz(-0.37014827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8106812) q[0];
sx q[0];
rz(-3.1157012) q[0];
sx q[0];
rz(2.6736892) q[0];
rz(-pi) q[1];
rz(3.0022718) q[2];
sx q[2];
rz(-1.3874386) q[2];
sx q[2];
rz(0.56688165) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83997516) q[1];
sx q[1];
rz(-0.8704005) q[1];
sx q[1];
rz(-0.61995929) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5657827) q[3];
sx q[3];
rz(-1.1429938) q[3];
sx q[3];
rz(1.4354777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7438573) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.7605304) q[2];
rz(-2.0928275) q[3];
sx q[3];
rz(-1.7966813) q[3];
sx q[3];
rz(-2.5019808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25310707) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(-2.2174368) q[0];
rz(-2.2249075) q[1];
sx q[1];
rz(-2.075383) q[1];
sx q[1];
rz(-1.7402657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9144672) q[0];
sx q[0];
rz(-1.0773398) q[0];
sx q[0];
rz(-1.6363793) q[0];
x q[1];
rz(-1.7426874) q[2];
sx q[2];
rz(-0.92280932) q[2];
sx q[2];
rz(-1.6164219) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7622213) q[1];
sx q[1];
rz(-1.4109542) q[1];
sx q[1];
rz(0.74121468) q[1];
rz(-pi) q[2];
rz(-1.0823233) q[3];
sx q[3];
rz(-1.3783749) q[3];
sx q[3];
rz(0.50979641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2283198) q[2];
sx q[2];
rz(-1.4618123) q[2];
sx q[2];
rz(-1.224996) q[2];
rz(-0.0082958881) q[3];
sx q[3];
rz(-0.010992916) q[3];
sx q[3];
rz(2.2744961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55815721) q[0];
sx q[0];
rz(-0.83034101) q[0];
sx q[0];
rz(0.48026568) q[0];
rz(-0.14446124) q[1];
sx q[1];
rz(-2.2339349) q[1];
sx q[1];
rz(-2.2772148) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1594161) q[0];
sx q[0];
rz(-1.5242531) q[0];
sx q[0];
rz(-1.7345364) q[0];
rz(-pi) q[1];
rz(-1.6046674) q[2];
sx q[2];
rz(-1.665417) q[2];
sx q[2];
rz(-2.692344) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0675812) q[1];
sx q[1];
rz(-2.1076822) q[1];
sx q[1];
rz(-0.28089653) q[1];
x q[2];
rz(-2.8671667) q[3];
sx q[3];
rz(-1.84019) q[3];
sx q[3];
rz(-1.5116949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20624837) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(0.63507357) q[2];
rz(-2.5687929) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(0.87535453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11005814) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(3.0173259) q[0];
rz(2.7763413) q[1];
sx q[1];
rz(-1.4271913) q[1];
sx q[1];
rz(1.431538) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99125049) q[0];
sx q[0];
rz(-1.7479154) q[0];
sx q[0];
rz(2.9464673) q[0];
rz(-pi) q[1];
rz(-0.65406873) q[2];
sx q[2];
rz(-2.0383012) q[2];
sx q[2];
rz(-1.2142177) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46043432) q[1];
sx q[1];
rz(-1.8755113) q[1];
sx q[1];
rz(-2.4439993) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67691524) q[3];
sx q[3];
rz(-0.64161086) q[3];
sx q[3];
rz(2.7144002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8753836) q[2];
sx q[2];
rz(-2.2201846) q[2];
sx q[2];
rz(1.9688152) q[2];
rz(-1.4069936) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(1.2146568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59239546) q[0];
sx q[0];
rz(-1.1347436) q[0];
sx q[0];
rz(2.5469575) q[0];
rz(1.0058962) q[1];
sx q[1];
rz(-1.9391831) q[1];
sx q[1];
rz(-0.88919052) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2834463) q[0];
sx q[0];
rz(-2.915417) q[0];
sx q[0];
rz(1.7986138) q[0];
rz(1.4290733) q[2];
sx q[2];
rz(-2.6772039) q[2];
sx q[2];
rz(-2.1187256) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8264849) q[1];
sx q[1];
rz(-1.6651622) q[1];
sx q[1];
rz(-3.0815275) q[1];
x q[2];
rz(0.066915705) q[3];
sx q[3];
rz(-0.51735462) q[3];
sx q[3];
rz(-2.2676615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95120007) q[2];
sx q[2];
rz(-1.9659646) q[2];
sx q[2];
rz(0.54538837) q[2];
rz(2.178318) q[3];
sx q[3];
rz(-1.1759718) q[3];
sx q[3];
rz(1.4639328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49190285) q[0];
sx q[0];
rz(-2.2408673) q[0];
sx q[0];
rz(1.3229205) q[0];
rz(-0.84359618) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(2.0100935) q[2];
sx q[2];
rz(-1.7514624) q[2];
sx q[2];
rz(-3.024586) q[2];
rz(2.5582377) q[3];
sx q[3];
rz(-1.1370549) q[3];
sx q[3];
rz(2.1255253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
