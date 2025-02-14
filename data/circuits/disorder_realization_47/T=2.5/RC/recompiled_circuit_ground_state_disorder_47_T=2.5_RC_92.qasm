OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.2831777) q[0];
sx q[0];
rz(-1.5812961) q[0];
sx q[0];
rz(2.4726782) q[0];
rz(0.34612292) q[1];
sx q[1];
rz(-2.3200413) q[1];
sx q[1];
rz(-0.93924826) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72059435) q[0];
sx q[0];
rz(-1.8587374) q[0];
sx q[0];
rz(1.1714515) q[0];
rz(-1.461636) q[2];
sx q[2];
rz(-2.7353035) q[2];
sx q[2];
rz(-3.0325505) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2841442) q[1];
sx q[1];
rz(-1.836187) q[1];
sx q[1];
rz(1.1737203) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9272618) q[3];
sx q[3];
rz(-0.89205304) q[3];
sx q[3];
rz(-1.9744622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5657438) q[2];
sx q[2];
rz(-2.1243024) q[2];
sx q[2];
rz(3.1080833) q[2];
rz(1.9401898) q[3];
sx q[3];
rz(-1.7824495) q[3];
sx q[3];
rz(2.0638154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031438436) q[0];
sx q[0];
rz(-2.8903676) q[0];
sx q[0];
rz(2.187619) q[0];
rz(-1.6961478) q[1];
sx q[1];
rz(-2.0868389) q[1];
sx q[1];
rz(-1.1711858) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32098371) q[0];
sx q[0];
rz(-2.1760328) q[0];
sx q[0];
rz(0.73802276) q[0];
rz(-1.8656125) q[2];
sx q[2];
rz(-2.1871242) q[2];
sx q[2];
rz(0.64586879) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3747159) q[1];
sx q[1];
rz(-1.7991814) q[1];
sx q[1];
rz(1.5414052) q[1];
x q[2];
rz(2.0061139) q[3];
sx q[3];
rz(-1.2655228) q[3];
sx q[3];
rz(2.0729617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90616068) q[2];
sx q[2];
rz(-1.7110598) q[2];
sx q[2];
rz(1.1492427) q[2];
rz(-3.1084642) q[3];
sx q[3];
rz(-1.5002316) q[3];
sx q[3];
rz(-0.59605014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0843435) q[0];
sx q[0];
rz(-0.79279041) q[0];
sx q[0];
rz(2.1113915) q[0];
rz(1.9016117) q[1];
sx q[1];
rz(-1.7844618) q[1];
sx q[1];
rz(-2.1485567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1112615) q[0];
sx q[0];
rz(-1.5227888) q[0];
sx q[0];
rz(2.6859849) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6898174) q[2];
sx q[2];
rz(-2.2348711) q[2];
sx q[2];
rz(1.9805976) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5807016) q[1];
sx q[1];
rz(-2.1296394) q[1];
sx q[1];
rz(2.2528095) q[1];
x q[2];
rz(2.7168324) q[3];
sx q[3];
rz(-0.83213193) q[3];
sx q[3];
rz(-2.5855541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21732907) q[2];
sx q[2];
rz(-2.0254878) q[2];
sx q[2];
rz(2.7868311) q[2];
rz(2.1445856) q[3];
sx q[3];
rz(-1.487178) q[3];
sx q[3];
rz(-0.94064373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724801) q[0];
sx q[0];
rz(-1.4262154) q[0];
sx q[0];
rz(-1.1068363) q[0];
rz(-2.2726982) q[1];
sx q[1];
rz(-0.51742253) q[1];
sx q[1];
rz(2.4443764) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91928673) q[0];
sx q[0];
rz(-0.59658748) q[0];
sx q[0];
rz(-0.35937341) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12733404) q[2];
sx q[2];
rz(-0.72605726) q[2];
sx q[2];
rz(-0.51148326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6706411) q[1];
sx q[1];
rz(-1.1837675) q[1];
sx q[1];
rz(2.7426038) q[1];
x q[2];
rz(1.5460062) q[3];
sx q[3];
rz(-1.5953334) q[3];
sx q[3];
rz(0.33543521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8640459) q[2];
sx q[2];
rz(-0.20865455) q[2];
sx q[2];
rz(-0.29402688) q[2];
rz(-2.0781519) q[3];
sx q[3];
rz(-1.682351) q[3];
sx q[3];
rz(2.3991876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71642891) q[0];
sx q[0];
rz(-0.18121457) q[0];
sx q[0];
rz(0.46318769) q[0];
rz(-0.40395346) q[1];
sx q[1];
rz(-1.5084167) q[1];
sx q[1];
rz(-2.892866) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8489784) q[0];
sx q[0];
rz(-0.19751829) q[0];
sx q[0];
rz(-1.5794483) q[0];
rz(1.8291446) q[2];
sx q[2];
rz(-1.0444006) q[2];
sx q[2];
rz(-1.0963944) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32837069) q[1];
sx q[1];
rz(-1.2036474) q[1];
sx q[1];
rz(0.35069938) q[1];
rz(-pi) q[2];
rz(0.71813138) q[3];
sx q[3];
rz(-2.2635824) q[3];
sx q[3];
rz(0.75961514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2948461) q[2];
sx q[2];
rz(-1.3174026) q[2];
sx q[2];
rz(-1.741629) q[2];
rz(1.8585662) q[3];
sx q[3];
rz(-1.3834407) q[3];
sx q[3];
rz(-2.9156901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95752174) q[0];
sx q[0];
rz(-2.3029843) q[0];
sx q[0];
rz(2.7958909) q[0];
rz(3.0905981) q[1];
sx q[1];
rz(-1.2268927) q[1];
sx q[1];
rz(-2.3695703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2564023) q[0];
sx q[0];
rz(-1.2202383) q[0];
sx q[0];
rz(1.8392483) q[0];
rz(-pi) q[1];
rz(3.0446485) q[2];
sx q[2];
rz(-0.92099316) q[2];
sx q[2];
rz(-2.52663) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.009338) q[1];
sx q[1];
rz(-1.8937832) q[1];
sx q[1];
rz(-2.1732974) q[1];
rz(-pi) q[2];
rz(2.8528522) q[3];
sx q[3];
rz(-1.0613228) q[3];
sx q[3];
rz(0.71392194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2766075) q[2];
sx q[2];
rz(-1.7165246) q[2];
sx q[2];
rz(2.9841606) q[2];
rz(0.43846798) q[3];
sx q[3];
rz(-0.74123588) q[3];
sx q[3];
rz(-2.1854775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(0.081414374) q[0];
sx q[0];
rz(-1.0670476) q[0];
sx q[0];
rz(9/(11*pi)) q[0];
rz(-0.57834894) q[1];
sx q[1];
rz(-1.3713505) q[1];
sx q[1];
rz(-2.9885805) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6575847) q[0];
sx q[0];
rz(-1.0562684) q[0];
sx q[0];
rz(0.57782509) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0566759) q[2];
sx q[2];
rz(-1.6667637) q[2];
sx q[2];
rz(1.3451479) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0205903) q[1];
sx q[1];
rz(-0.7210702) q[1];
sx q[1];
rz(0.64863689) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97665321) q[3];
sx q[3];
rz(-0.53721957) q[3];
sx q[3];
rz(1.2402759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1795307) q[2];
sx q[2];
rz(-0.099055812) q[2];
sx q[2];
rz(2.8774101) q[2];
rz(-1.8386486) q[3];
sx q[3];
rz(-1.7641726) q[3];
sx q[3];
rz(-1.1072268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95689479) q[0];
sx q[0];
rz(-0.95320025) q[0];
sx q[0];
rz(-1.4554998) q[0];
rz(-2.1799344) q[1];
sx q[1];
rz(-1.7627629) q[1];
sx q[1];
rz(-0.89967322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1048516) q[0];
sx q[0];
rz(-2.0213958) q[0];
sx q[0];
rz(1.3356871) q[0];
rz(-pi) q[1];
rz(1.3191965) q[2];
sx q[2];
rz(-0.74079047) q[2];
sx q[2];
rz(-2.0513926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4845061) q[1];
sx q[1];
rz(-2.5528347) q[1];
sx q[1];
rz(-1.6480873) q[1];
rz(1.7877696) q[3];
sx q[3];
rz(-1.5343175) q[3];
sx q[3];
rz(-2.77751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6962894) q[2];
sx q[2];
rz(-0.60089198) q[2];
sx q[2];
rz(-0.75330934) q[2];
rz(0.2937915) q[3];
sx q[3];
rz(-1.1283504) q[3];
sx q[3];
rz(-1.1118579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488572) q[0];
sx q[0];
rz(-2.1519372) q[0];
sx q[0];
rz(-2.2897172) q[0];
rz(1.4082255) q[1];
sx q[1];
rz(-1.6586761) q[1];
sx q[1];
rz(-2.7588989) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9155884) q[0];
sx q[0];
rz(-1.3289641) q[0];
sx q[0];
rz(2.7126724) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23613249) q[2];
sx q[2];
rz(-2.3575767) q[2];
sx q[2];
rz(1.2003984) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5671317) q[1];
sx q[1];
rz(-1.8715011) q[1];
sx q[1];
rz(1.6733132) q[1];
x q[2];
rz(-2.6555041) q[3];
sx q[3];
rz(-1.6440142) q[3];
sx q[3];
rz(1.40772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7644299) q[2];
sx q[2];
rz(-2.1775776) q[2];
sx q[2];
rz(-2.5386179) q[2];
rz(-2.073334) q[3];
sx q[3];
rz(-1.4385185) q[3];
sx q[3];
rz(-2.300613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.1339486) q[0];
sx q[0];
rz(-1.4643865) q[0];
sx q[0];
rz(0.35368791) q[0];
rz(1.1324646) q[1];
sx q[1];
rz(-2.1734889) q[1];
sx q[1];
rz(-0.38945928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48943993) q[0];
sx q[0];
rz(-1.8021694) q[0];
sx q[0];
rz(-2.6315297) q[0];
rz(-pi) q[1];
rz(3.0385397) q[2];
sx q[2];
rz(-1.3259058) q[2];
sx q[2];
rz(0.51765294) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8824401) q[1];
sx q[1];
rz(-2.5969567) q[1];
sx q[1];
rz(2.5914521) q[1];
rz(0.52822379) q[3];
sx q[3];
rz(-1.6162795) q[3];
sx q[3];
rz(0.32873617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0103717) q[2];
sx q[2];
rz(-1.6348569) q[2];
sx q[2];
rz(1.9949251) q[2];
rz(1.6049339) q[3];
sx q[3];
rz(-2.6585572) q[3];
sx q[3];
rz(0.35227942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1491886) q[0];
sx q[0];
rz(-2.438899) q[0];
sx q[0];
rz(-2.6371523) q[0];
rz(3.0814677) q[1];
sx q[1];
rz(-0.45777121) q[1];
sx q[1];
rz(-0.4578185) q[1];
rz(-1.6788775) q[2];
sx q[2];
rz(-1.6440132) q[2];
sx q[2];
rz(-3.0655412) q[2];
rz(0.89636421) q[3];
sx q[3];
rz(-1.5980362) q[3];
sx q[3];
rz(-2.7841795) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
