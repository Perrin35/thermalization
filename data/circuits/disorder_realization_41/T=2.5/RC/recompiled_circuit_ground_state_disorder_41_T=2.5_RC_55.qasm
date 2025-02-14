OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1908258) q[0];
sx q[0];
rz(-1.289225) q[0];
sx q[0];
rz(3.0486795) q[0];
rz(3.0463123) q[1];
sx q[1];
rz(-2.4089101) q[1];
sx q[1];
rz(-1.9021775) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5372972) q[0];
sx q[0];
rz(-1.3193466) q[0];
sx q[0];
rz(2.9049113) q[0];
rz(-1.3746512) q[2];
sx q[2];
rz(-1.8298379) q[2];
sx q[2];
rz(2.1262622) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2762294) q[1];
sx q[1];
rz(-1.6003834) q[1];
sx q[1];
rz(-0.33025708) q[1];
rz(0.57115023) q[3];
sx q[3];
rz(-1.1451654) q[3];
sx q[3];
rz(-1.7895123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6713082) q[2];
sx q[2];
rz(-1.6613864) q[2];
sx q[2];
rz(-1.3489464) q[2];
rz(0.74364439) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(2.9644137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.141356) q[0];
sx q[0];
rz(-1.5970705) q[0];
sx q[0];
rz(2.8479688) q[0];
rz(-1.0307505) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(-2.7986599) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9259277) q[0];
sx q[0];
rz(-2.5747888) q[0];
sx q[0];
rz(-1.180992) q[0];
x q[1];
rz(0.12983506) q[2];
sx q[2];
rz(-1.9399683) q[2];
sx q[2];
rz(-0.8952199) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82417497) q[1];
sx q[1];
rz(-1.3567827) q[1];
sx q[1];
rz(0.70862464) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2003056) q[3];
sx q[3];
rz(-1.3664748) q[3];
sx q[3];
rz(-0.91703892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12884101) q[2];
sx q[2];
rz(-1.7495456) q[2];
sx q[2];
rz(-0.7592321) q[2];
rz(3.1208842) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52221209) q[0];
sx q[0];
rz(-0.601957) q[0];
sx q[0];
rz(3.1371327) q[0];
rz(0.64747539) q[1];
sx q[1];
rz(-0.59440333) q[1];
sx q[1];
rz(-2.3562145) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42703907) q[0];
sx q[0];
rz(-1.4903965) q[0];
sx q[0];
rz(0.15958448) q[0];
x q[1];
rz(2.9663229) q[2];
sx q[2];
rz(-1.2614377) q[2];
sx q[2];
rz(-0.64495211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8632311) q[1];
sx q[1];
rz(-1.0309217) q[1];
sx q[1];
rz(3.0616772) q[1];
rz(-pi) q[2];
rz(0.66297279) q[3];
sx q[3];
rz(-1.0043274) q[3];
sx q[3];
rz(2.278355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7711827) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(1.3578337) q[2];
rz(2.8010098) q[3];
sx q[3];
rz(-0.95652306) q[3];
sx q[3];
rz(-0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1444645) q[0];
sx q[0];
rz(-1.0839394) q[0];
sx q[0];
rz(-2.2564364) q[0];
rz(0.74686933) q[1];
sx q[1];
rz(-0.62364686) q[1];
sx q[1];
rz(-2.0451827) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88056662) q[0];
sx q[0];
rz(-1.9184904) q[0];
sx q[0];
rz(0.33190042) q[0];
x q[1];
rz(0.030641067) q[2];
sx q[2];
rz(-1.5113792) q[2];
sx q[2];
rz(-2.5031646) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1433574) q[1];
sx q[1];
rz(-1.1117742) q[1];
sx q[1];
rz(1.8602536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.40989) q[3];
sx q[3];
rz(-1.1771132) q[3];
sx q[3];
rz(1.7997774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61802822) q[2];
sx q[2];
rz(-0.21054331) q[2];
sx q[2];
rz(0.024624126) q[2];
rz(2.5060182) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(-2.2583101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6898952) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(-0.46689335) q[0];
rz(3.0310071) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(2.0707524) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12884451) q[0];
sx q[0];
rz(-1.8086047) q[0];
sx q[0];
rz(-1.1046011) q[0];
rz(0.51009615) q[2];
sx q[2];
rz(-0.3613216) q[2];
sx q[2];
rz(-2.9053743) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3339798) q[1];
sx q[1];
rz(-2.4845504) q[1];
sx q[1];
rz(0.64589898) q[1];
x q[2];
rz(-2.4038845) q[3];
sx q[3];
rz(-2.2900351) q[3];
sx q[3];
rz(-1.4054474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.11833) q[2];
sx q[2];
rz(-1.7795965) q[2];
sx q[2];
rz(0.85232097) q[2];
rz(0.037633745) q[3];
sx q[3];
rz(-1.9084385) q[3];
sx q[3];
rz(0.5948624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0090050176) q[0];
sx q[0];
rz(-1.9629033) q[0];
sx q[0];
rz(1.8527385) q[0];
rz(1.5124849) q[1];
sx q[1];
rz(-2.0178724) q[1];
sx q[1];
rz(1.1423133) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2347533) q[0];
sx q[0];
rz(-1.2182353) q[0];
sx q[0];
rz(-1.4666739) q[0];
rz(3.0794607) q[2];
sx q[2];
rz(-1.3634063) q[2];
sx q[2];
rz(-2.1211989) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21396046) q[1];
sx q[1];
rz(-2.7134656) q[1];
sx q[1];
rz(2.5090748) q[1];
rz(1.3942424) q[3];
sx q[3];
rz(-1.8774967) q[3];
sx q[3];
rz(-0.97489417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(-0.12040559) q[2];
rz(-1.1558007) q[3];
sx q[3];
rz(-1.5294231) q[3];
sx q[3];
rz(2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33893809) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(-1.4539723) q[0];
rz(-1.5441719) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(2.6905751) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9748443) q[0];
sx q[0];
rz(-1.5112425) q[0];
sx q[0];
rz(-0.3015428) q[0];
x q[1];
rz(-2.9092009) q[2];
sx q[2];
rz(-1.2133779) q[2];
sx q[2];
rz(-3.0136303) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4828795) q[1];
sx q[1];
rz(-0.45156839) q[1];
sx q[1];
rz(0.2803448) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3312206) q[3];
sx q[3];
rz(-1.5640508) q[3];
sx q[3];
rz(0.41616671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1958127) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(1.3607402) q[2];
rz(2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(1.3895234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1239531) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(-0.70607591) q[0];
rz(-1.5977244) q[1];
sx q[1];
rz(-1.4223301) q[1];
sx q[1];
rz(-2.2854038) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5917011) q[0];
sx q[0];
rz(-0.020792637) q[0];
sx q[0];
rz(-1.9338305) q[0];
rz(-0.82243408) q[2];
sx q[2];
rz(-1.8959799) q[2];
sx q[2];
rz(2.6564244) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4256712) q[1];
sx q[1];
rz(-1.3678663) q[1];
sx q[1];
rz(-0.58260609) q[1];
rz(-pi) q[2];
rz(2.059465) q[3];
sx q[3];
rz(-1.4525982) q[3];
sx q[3];
rz(-1.4335872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84264821) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(2.974158) q[2];
rz(-1.5873448) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(-2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3779959) q[0];
sx q[0];
rz(-1.6908228) q[0];
sx q[0];
rz(-2.7698621) q[0];
rz(1.4338214) q[1];
sx q[1];
rz(-1.6038409) q[1];
sx q[1];
rz(0.30002123) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8512631) q[0];
sx q[0];
rz(-1.4545049) q[0];
sx q[0];
rz(-0.28684464) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6113564) q[2];
sx q[2];
rz(-0.99405655) q[2];
sx q[2];
rz(2.4941913) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3206818) q[1];
sx q[1];
rz(-1.6947985) q[1];
sx q[1];
rz(0.65202608) q[1];
rz(3.0172273) q[3];
sx q[3];
rz(-2.2699286) q[3];
sx q[3];
rz(1.3209526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6649449) q[2];
sx q[2];
rz(-0.46755329) q[2];
sx q[2];
rz(0.42759582) q[2];
rz(-1.556373) q[3];
sx q[3];
rz(-1.1170324) q[3];
sx q[3];
rz(0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72117358) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(-2.6721201) q[0];
rz(0.92533127) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(-0.023177711) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4172702) q[0];
sx q[0];
rz(-0.5529772) q[0];
sx q[0];
rz(-2.824845) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2397348) q[2];
sx q[2];
rz(-2.0624071) q[2];
sx q[2];
rz(0.41504809) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7851023) q[1];
sx q[1];
rz(-1.7447326) q[1];
sx q[1];
rz(0.16696232) q[1];
x q[2];
rz(2.9523207) q[3];
sx q[3];
rz(-2.2918275) q[3];
sx q[3];
rz(1.4355575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8699441) q[2];
sx q[2];
rz(-1.9359438) q[2];
sx q[2];
rz(-0.49087697) q[2];
rz(-1.0902181) q[3];
sx q[3];
rz(-0.96929437) q[3];
sx q[3];
rz(2.8276665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8563817) q[0];
sx q[0];
rz(-2.2638392) q[0];
sx q[0];
rz(1.0765156) q[0];
rz(-2.8648227) q[1];
sx q[1];
rz(-2.3634187) q[1];
sx q[1];
rz(1.707911) q[1];
rz(-2.4527373) q[2];
sx q[2];
rz(-1.8102472) q[2];
sx q[2];
rz(2.3285338) q[2];
rz(3.1076486) q[3];
sx q[3];
rz(-0.53474075) q[3];
sx q[3];
rz(-0.026215601) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
