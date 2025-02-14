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
rz(1.1531416) q[0];
sx q[0];
rz(-0.81557953) q[0];
sx q[0];
rz(2.3834035) q[0];
rz(-0.012501333) q[1];
sx q[1];
rz(-1.3036417) q[1];
sx q[1];
rz(1.5703896) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2822097) q[0];
sx q[0];
rz(-0.73732418) q[0];
sx q[0];
rz(-0.54873772) q[0];
x q[1];
rz(1.1711898) q[2];
sx q[2];
rz(-1.5087224) q[2];
sx q[2];
rz(-2.2860179) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6627854) q[1];
sx q[1];
rz(-2.7764455) q[1];
sx q[1];
rz(-1.410233) q[1];
rz(1.2615573) q[3];
sx q[3];
rz(-1.6282363) q[3];
sx q[3];
rz(-0.32306898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6115173) q[2];
sx q[2];
rz(-3.1334183) q[2];
sx q[2];
rz(-0.44069904) q[2];
rz(3.0618482) q[3];
sx q[3];
rz(-0.00010448797) q[3];
sx q[3];
rz(1.124148) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1752862) q[0];
sx q[0];
rz(-0.056862406) q[0];
sx q[0];
rz(-2.9735907) q[0];
rz(-3.1213144) q[1];
sx q[1];
rz(-2.8333277) q[1];
sx q[1];
rz(-1.6049989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8533596) q[0];
sx q[0];
rz(-2.8194966) q[0];
sx q[0];
rz(-2.3951985) q[0];
rz(-pi) q[1];
rz(3.1390983) q[2];
sx q[2];
rz(-1.5606784) q[2];
sx q[2];
rz(-1.5960787) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55024338) q[1];
sx q[1];
rz(-0.004799407) q[1];
sx q[1];
rz(-1.8018434) q[1];
rz(-pi) q[2];
x q[2];
rz(0.042293799) q[3];
sx q[3];
rz(-1.5618069) q[3];
sx q[3];
rz(-0.90153722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7961879) q[2];
sx q[2];
rz(-0.91190839) q[2];
sx q[2];
rz(-1.7339285) q[2];
rz(1.0450854) q[3];
sx q[3];
rz(-0.049523517) q[3];
sx q[3];
rz(-2.869587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8318091) q[0];
sx q[0];
rz(-2.164916) q[0];
sx q[0];
rz(-2.580544) q[0];
rz(0.27684119) q[1];
sx q[1];
rz(-3.1287153) q[1];
sx q[1];
rz(1.3078088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5878752) q[0];
sx q[0];
rz(-0.29969117) q[0];
sx q[0];
rz(1.7100348) q[0];
rz(-pi) q[1];
rz(3.1409227) q[2];
sx q[2];
rz(-3.1342034) q[2];
sx q[2];
rz(2.0787897) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7737116) q[1];
sx q[1];
rz(-2.5653172) q[1];
sx q[1];
rz(-1.4632844) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6827319) q[3];
sx q[3];
rz(-0.94597497) q[3];
sx q[3];
rz(0.030908728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7967367) q[2];
sx q[2];
rz(-3.1414746) q[2];
sx q[2];
rz(-2.5522088) q[2];
rz(1.899259) q[3];
sx q[3];
rz(-3.1292249) q[3];
sx q[3];
rz(1.3602268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0068483343) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(-1.7929329) q[0];
rz(-3.1349365) q[1];
sx q[1];
rz(-1.8204047) q[1];
sx q[1];
rz(-0.033500813) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071246192) q[0];
sx q[0];
rz(-1.7435562) q[0];
sx q[0];
rz(2.2281649) q[0];
rz(-pi) q[1];
rz(-3.1042023) q[2];
sx q[2];
rz(-3.0218736) q[2];
sx q[2];
rz(0.38116383) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0211612) q[1];
sx q[1];
rz(-0.26909262) q[1];
sx q[1];
rz(-1.6090367) q[1];
rz(2.2073645) q[3];
sx q[3];
rz(-0.2033955) q[3];
sx q[3];
rz(2.0980841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.032430705) q[2];
sx q[2];
rz(-0.0065294821) q[2];
sx q[2];
rz(-2.8429441) q[2];
rz(-1.2700891) q[3];
sx q[3];
rz(-0.01615571) q[3];
sx q[3];
rz(3.0884009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.771516) q[0];
sx q[0];
rz(-1.6271485) q[0];
sx q[0];
rz(-2.694743) q[0];
rz(2.9554548) q[1];
sx q[1];
rz(-0.061807241) q[1];
sx q[1];
rz(1.4131379) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7639358) q[0];
sx q[0];
rz(-1.3400338) q[0];
sx q[0];
rz(-3.1070437) q[0];
rz(-2.9252626) q[2];
sx q[2];
rz(-1.9842752) q[2];
sx q[2];
rz(2.095393) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6525938) q[1];
sx q[1];
rz(-1.6203383) q[1];
sx q[1];
rz(0.045157305) q[1];
rz(-pi) q[2];
rz(1.8522315) q[3];
sx q[3];
rz(-1.7858943) q[3];
sx q[3];
rz(-0.18751442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3073005) q[2];
sx q[2];
rz(-1.5520381) q[2];
sx q[2];
rz(-2.6429122) q[2];
rz(-2.5636766) q[3];
sx q[3];
rz(-0.48360616) q[3];
sx q[3];
rz(-0.59670603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96108288) q[0];
sx q[0];
rz(-1.120765) q[0];
sx q[0];
rz(-2.7951796) q[0];
rz(-0.60180426) q[1];
sx q[1];
rz(-1.5806942) q[1];
sx q[1];
rz(-2.3880889) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6974073) q[0];
sx q[0];
rz(-1.7427398) q[0];
sx q[0];
rz(2.9471022) q[0];
rz(-3.0664938) q[2];
sx q[2];
rz(-1.4595928) q[2];
sx q[2];
rz(0.62125833) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9648887) q[1];
sx q[1];
rz(-0.93550013) q[1];
sx q[1];
rz(-2.6179362) q[1];
rz(1.7338792) q[3];
sx q[3];
rz(-1.661851) q[3];
sx q[3];
rz(-1.0331105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.57009131) q[2];
sx q[2];
rz(-0.0034905958) q[2];
sx q[2];
rz(-1.6174779) q[2];
rz(0.12024719) q[3];
sx q[3];
rz(-0.0032987981) q[3];
sx q[3];
rz(2.6038468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26871249) q[0];
sx q[0];
rz(-2.217642) q[0];
sx q[0];
rz(0.22126108) q[0];
rz(-1.4606754) q[1];
sx q[1];
rz(-0.9333846) q[1];
sx q[1];
rz(-3.0642919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6344096) q[0];
sx q[0];
rz(-1.612525) q[0];
sx q[0];
rz(3.1401278) q[0];
rz(-pi) q[1];
rz(-0.004782025) q[2];
sx q[2];
rz(-1.5796697) q[2];
sx q[2];
rz(0.21172571) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0127924) q[1];
sx q[1];
rz(-2.9523053) q[1];
sx q[1];
rz(-0.38893338) q[1];
rz(-pi) q[2];
rz(3.0139913) q[3];
sx q[3];
rz(-1.7831047) q[3];
sx q[3];
rz(-1.9636167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7920502) q[2];
sx q[2];
rz(-0.011186102) q[2];
sx q[2];
rz(-0.95996094) q[2];
rz(-2.8121484) q[3];
sx q[3];
rz(-3.1335148) q[3];
sx q[3];
rz(0.85028696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9462117) q[0];
sx q[0];
rz(-2.52849) q[0];
sx q[0];
rz(-0.10928133) q[0];
rz(-2.7682313) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(1.9111309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41540158) q[0];
sx q[0];
rz(-2.2080732) q[0];
sx q[0];
rz(0.7755875) q[0];
rz(-2.3745499) q[2];
sx q[2];
rz(-0.27077507) q[2];
sx q[2];
rz(-0.8052288) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8671678) q[1];
sx q[1];
rz(-1.503957) q[1];
sx q[1];
rz(-0.011456077) q[1];
rz(1.968574) q[3];
sx q[3];
rz(-2.4511485) q[3];
sx q[3];
rz(-3.0544314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5665148) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(1.8129978) q[2];
rz(1.7447507) q[3];
sx q[3];
rz(-3.1378523) q[3];
sx q[3];
rz(-1.0218792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063865572) q[0];
sx q[0];
rz(-1.6985748) q[0];
sx q[0];
rz(2.5685837) q[0];
rz(0.30814463) q[1];
sx q[1];
rz(-0.40987086) q[1];
sx q[1];
rz(2.134197) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0125748) q[0];
sx q[0];
rz(-1.5652302) q[0];
sx q[0];
rz(1.6773423) q[0];
rz(3.0976866) q[2];
sx q[2];
rz(-1.708235) q[2];
sx q[2];
rz(-0.034417987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7511661) q[1];
sx q[1];
rz(-1.5301643) q[1];
sx q[1];
rz(-1.4875814) q[1];
x q[2];
rz(3.1281578) q[3];
sx q[3];
rz(-1.5406431) q[3];
sx q[3];
rz(2.9363901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8241626) q[2];
sx q[2];
rz(-2.514826) q[2];
sx q[2];
rz(2.7516348) q[2];
rz(3.0725078) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(-2.8100584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1214509) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(0.48594117) q[0];
rz(-2.2700229) q[1];
sx q[1];
rz(-1.3078682) q[1];
sx q[1];
rz(-1.4942687) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.082162) q[0];
sx q[0];
rz(-2.5100187) q[0];
sx q[0];
rz(0.53379121) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61482347) q[2];
sx q[2];
rz(-1.5346926) q[2];
sx q[2];
rz(1.6143798) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4758265) q[1];
sx q[1];
rz(-1.2428871) q[1];
sx q[1];
rz(1.2535415) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6127729) q[3];
sx q[3];
rz(-1.6838361) q[3];
sx q[3];
rz(2.6338995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5747052) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(-3.1097143) q[2];
rz(2.3632862) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(0.2955029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42302172) q[0];
sx q[0];
rz(-1.5324677) q[0];
sx q[0];
rz(1.8146429) q[0];
rz(3.0146535) q[1];
sx q[1];
rz(-2.9025684) q[1];
sx q[1];
rz(-2.92166) q[1];
rz(-3.1359966) q[2];
sx q[2];
rz(-1.4311309) q[2];
sx q[2];
rz(-2.8972814) q[2];
rz(-1.4999785) q[3];
sx q[3];
rz(-0.67588617) q[3];
sx q[3];
rz(-2.5839154) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
