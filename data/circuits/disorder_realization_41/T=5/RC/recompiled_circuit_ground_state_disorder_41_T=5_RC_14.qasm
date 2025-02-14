OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5858894) q[0];
sx q[0];
rz(-1.1893505) q[0];
sx q[0];
rz(1.9777634) q[0];
rz(-0.7467421) q[1];
sx q[1];
rz(-0.26489869) q[1];
sx q[1];
rz(0.17118153) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8876182) q[0];
sx q[0];
rz(-1.7897535) q[0];
sx q[0];
rz(-0.15677126) q[0];
rz(-2.3642478) q[2];
sx q[2];
rz(-1.0413578) q[2];
sx q[2];
rz(-2.6482761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92058006) q[1];
sx q[1];
rz(-1.5713931) q[1];
sx q[1];
rz(1.8662054) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5653651) q[3];
sx q[3];
rz(-1.5539546) q[3];
sx q[3];
rz(-1.668949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.091398679) q[2];
sx q[2];
rz(-1.5768496) q[2];
sx q[2];
rz(-1.0696627) q[2];
rz(1.7012677) q[3];
sx q[3];
rz(-1.0245208) q[3];
sx q[3];
rz(1.9694156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1759724) q[0];
sx q[0];
rz(-1.6930641) q[0];
sx q[0];
rz(0.58797055) q[0];
rz(-1.8486456) q[1];
sx q[1];
rz(-2.1473532) q[1];
sx q[1];
rz(0.15486823) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42255369) q[0];
sx q[0];
rz(-2.6309359) q[0];
sx q[0];
rz(0.22061781) q[0];
rz(1.7172408) q[2];
sx q[2];
rz(-1.3322209) q[2];
sx q[2];
rz(-1.5244791) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1701291) q[1];
sx q[1];
rz(-1.9057353) q[1];
sx q[1];
rz(1.0748802) q[1];
rz(-pi) q[2];
rz(-2.5741379) q[3];
sx q[3];
rz(-2.8914872) q[3];
sx q[3];
rz(-0.21931822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4984442) q[2];
sx q[2];
rz(-2.4109106) q[2];
sx q[2];
rz(-3.0291962) q[2];
rz(-0.82320881) q[3];
sx q[3];
rz(-1.9197437) q[3];
sx q[3];
rz(-2.6196151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4560029) q[0];
sx q[0];
rz(-1.0887479) q[0];
sx q[0];
rz(-0.96813694) q[0];
rz(1.118842) q[1];
sx q[1];
rz(-1.6066931) q[1];
sx q[1];
rz(1.8082089) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3385425) q[0];
sx q[0];
rz(-0.97369872) q[0];
sx q[0];
rz(-0.19495585) q[0];
rz(-1.1022262) q[2];
sx q[2];
rz(-1.0194687) q[2];
sx q[2];
rz(-2.9561549) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1044429) q[1];
sx q[1];
rz(-2.5698476) q[1];
sx q[1];
rz(1.7556095) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3465856) q[3];
sx q[3];
rz(-2.7751659) q[3];
sx q[3];
rz(2.623705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94990388) q[2];
sx q[2];
rz(-2.1046941) q[2];
sx q[2];
rz(2.8766768) q[2];
rz(2.8252937) q[3];
sx q[3];
rz(-2.8079872) q[3];
sx q[3];
rz(1.3914289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3259657) q[0];
sx q[0];
rz(-2.9209324) q[0];
sx q[0];
rz(1.6478446) q[0];
rz(1.7296467) q[1];
sx q[1];
rz(-1.0271881) q[1];
sx q[1];
rz(2.484201) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18497224) q[0];
sx q[0];
rz(-1.4669384) q[0];
sx q[0];
rz(1.4968027) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78423402) q[2];
sx q[2];
rz(-1.1335229) q[2];
sx q[2];
rz(-1.0904877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3616391) q[1];
sx q[1];
rz(-1.3758389) q[1];
sx q[1];
rz(-3.1038398) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3300276) q[3];
sx q[3];
rz(-2.2144284) q[3];
sx q[3];
rz(-2.0763458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.21939453) q[2];
sx q[2];
rz(-1.2224835) q[2];
sx q[2];
rz(-0.40327367) q[2];
rz(1.203631) q[3];
sx q[3];
rz(-1.5091242) q[3];
sx q[3];
rz(-2.6327314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4598684) q[0];
sx q[0];
rz(-2.7175856) q[0];
sx q[0];
rz(-2.5176609) q[0];
rz(-2.4195747) q[1];
sx q[1];
rz(-1.3474418) q[1];
sx q[1];
rz(3.0375979) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926804) q[0];
sx q[0];
rz(-0.70250547) q[0];
sx q[0];
rz(-2.9215432) q[0];
x q[1];
rz(1.9172759) q[2];
sx q[2];
rz(-1.9110381) q[2];
sx q[2];
rz(0.97803309) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.69332921) q[1];
sx q[1];
rz(-0.86719705) q[1];
sx q[1];
rz(-2.2660822) q[1];
x q[2];
rz(0.73443074) q[3];
sx q[3];
rz(-1.0356257) q[3];
sx q[3];
rz(1.2169692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7907052) q[2];
sx q[2];
rz(-0.71983379) q[2];
sx q[2];
rz(-0.78901115) q[2];
rz(-1.2645432) q[3];
sx q[3];
rz(-1.0435373) q[3];
sx q[3];
rz(2.1544382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3037381) q[0];
sx q[0];
rz(-0.721295) q[0];
sx q[0];
rz(-0.20027941) q[0];
rz(-1.4664949) q[1];
sx q[1];
rz(-1.3267582) q[1];
sx q[1];
rz(1.5178348) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90676722) q[0];
sx q[0];
rz(-1.2041766) q[0];
sx q[0];
rz(0.54327048) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0671704) q[2];
sx q[2];
rz(-1.4563515) q[2];
sx q[2];
rz(-1.9527854) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1863757) q[1];
sx q[1];
rz(-2.5506335) q[1];
sx q[1];
rz(0.0062828961) q[1];
x q[2];
rz(-0.79849859) q[3];
sx q[3];
rz(-1.3991465) q[3];
sx q[3];
rz(0.6397748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.62310654) q[2];
sx q[2];
rz(-0.97272626) q[2];
sx q[2];
rz(-0.34131193) q[2];
rz(0.85706472) q[3];
sx q[3];
rz(-1.9144446) q[3];
sx q[3];
rz(-0.46522337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.051006) q[0];
sx q[0];
rz(-1.0336646) q[0];
sx q[0];
rz(1.8940014) q[0];
rz(0.11820758) q[1];
sx q[1];
rz(-1.8659614) q[1];
sx q[1];
rz(3.1059713) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8310332) q[0];
sx q[0];
rz(-0.11971029) q[0];
sx q[0];
rz(1.381078) q[0];
rz(-2.9526677) q[2];
sx q[2];
rz(-1.9686724) q[2];
sx q[2];
rz(2.973345) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5069711) q[1];
sx q[1];
rz(-1.8160607) q[1];
sx q[1];
rz(2.5140106) q[1];
rz(-pi) q[2];
rz(2.3603975) q[3];
sx q[3];
rz(-1.6577814) q[3];
sx q[3];
rz(-1.9613822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9901765) q[2];
sx q[2];
rz(-1.4750865) q[2];
sx q[2];
rz(-2.0659633) q[2];
rz(0.52078024) q[3];
sx q[3];
rz(-0.62028378) q[3];
sx q[3];
rz(-1.5889997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2369279) q[0];
sx q[0];
rz(-1.0228782) q[0];
sx q[0];
rz(0.9032332) q[0];
rz(1.5029933) q[1];
sx q[1];
rz(-1.6939949) q[1];
sx q[1];
rz(-1.2618056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.107511) q[0];
sx q[0];
rz(-0.26824646) q[0];
sx q[0];
rz(-2.111042) q[0];
x q[1];
rz(2.4440358) q[2];
sx q[2];
rz(-0.64854014) q[2];
sx q[2];
rz(-1.4319789) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8324229) q[1];
sx q[1];
rz(-1.5201525) q[1];
sx q[1];
rz(3.0502351) q[1];
rz(1.0699959) q[3];
sx q[3];
rz(-1.0251932) q[3];
sx q[3];
rz(0.16014447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.5791851) q[2];
sx q[2];
rz(-0.6520485) q[2];
sx q[2];
rz(-0.11422608) q[2];
rz(-1.7646029) q[3];
sx q[3];
rz(-1.1403133) q[3];
sx q[3];
rz(-2.637114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1219516) q[0];
sx q[0];
rz(-0.99948245) q[0];
sx q[0];
rz(1.1720538) q[0];
rz(1.3837586) q[1];
sx q[1];
rz(-1.9969321) q[1];
sx q[1];
rz(-0.11464548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4534476) q[0];
sx q[0];
rz(-2.3575961) q[0];
sx q[0];
rz(2.408124) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9681712) q[2];
sx q[2];
rz(-2.6422814) q[2];
sx q[2];
rz(-0.87022802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9108565) q[1];
sx q[1];
rz(-1.641526) q[1];
sx q[1];
rz(1.8190228) q[1];
x q[2];
rz(-1.678336) q[3];
sx q[3];
rz(-2.3643374) q[3];
sx q[3];
rz(-1.4449256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5961479) q[2];
sx q[2];
rz(-0.20716509) q[2];
sx q[2];
rz(0.88085112) q[2];
rz(2.7018069) q[3];
sx q[3];
rz(-1.1652911) q[3];
sx q[3];
rz(-2.0161276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3038444) q[0];
sx q[0];
rz(-1.4864018) q[0];
sx q[0];
rz(0.098175511) q[0];
rz(-2.5671666) q[1];
sx q[1];
rz(-1.4593294) q[1];
sx q[1];
rz(-2.548545) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2132872) q[0];
sx q[0];
rz(-0.79003626) q[0];
sx q[0];
rz(1.9844878) q[0];
rz(-2.7213989) q[2];
sx q[2];
rz(-2.0717142) q[2];
sx q[2];
rz(-0.80244697) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8349466) q[1];
sx q[1];
rz(-3.1241841) q[1];
sx q[1];
rz(-2.2724292) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2393059) q[3];
sx q[3];
rz(-0.96486366) q[3];
sx q[3];
rz(0.87406033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1526996) q[2];
sx q[2];
rz(-1.1836735) q[2];
sx q[2];
rz(-0.93842554) q[2];
rz(1.3931795) q[3];
sx q[3];
rz(-1.3104855) q[3];
sx q[3];
rz(2.9032629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68168454) q[0];
sx q[0];
rz(-0.60612283) q[0];
sx q[0];
rz(-1.7932307) q[0];
rz(-2.0943191) q[1];
sx q[1];
rz(-0.92887639) q[1];
sx q[1];
rz(2.1705719) q[1];
rz(2.1793096) q[2];
sx q[2];
rz(-0.078799993) q[2];
sx q[2];
rz(-2.7344631) q[2];
rz(-1.3503475) q[3];
sx q[3];
rz(-1.7488283) q[3];
sx q[3];
rz(-1.8976952) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
