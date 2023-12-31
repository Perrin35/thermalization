OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(-2.321474) q[0];
rz(2.7804873) q[1];
sx q[1];
rz(-0.63280025) q[1];
sx q[1];
rz(-0.83067218) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9554064) q[0];
sx q[0];
rz(-2.6468228) q[0];
sx q[0];
rz(0.68899378) q[0];
rz(2.6248706) q[2];
sx q[2];
rz(-0.80846723) q[2];
sx q[2];
rz(2.7441623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4315011) q[1];
sx q[1];
rz(-2.3570865) q[1];
sx q[1];
rz(2.0246519) q[1];
rz(-1.9740231) q[3];
sx q[3];
rz(-1.1551757) q[3];
sx q[3];
rz(0.88413873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(-1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(-2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-2.2297915) q[0];
rz(2.3520825) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(0.3266913) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2355736) q[0];
sx q[0];
rz(-2.2740907) q[0];
sx q[0];
rz(-2.2554382) q[0];
x q[1];
rz(-2.0147444) q[2];
sx q[2];
rz(-1.7146646) q[2];
sx q[2];
rz(-2.4614046) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82636278) q[1];
sx q[1];
rz(-2.3781207) q[1];
sx q[1];
rz(-2.1748494) q[1];
x q[2];
rz(3.0558415) q[3];
sx q[3];
rz(-2.2567856) q[3];
sx q[3];
rz(0.12967295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6313173) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(3.0316947) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(-1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.592955) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(2.8116995) q[0];
rz(0.27711162) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-1.057391) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43305106) q[0];
sx q[0];
rz(-1.9110702) q[0];
sx q[0];
rz(0.021854594) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6234987) q[2];
sx q[2];
rz(-2.7535541) q[2];
sx q[2];
rz(-0.54259091) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96506572) q[1];
sx q[1];
rz(-2.2245193) q[1];
sx q[1];
rz(-2.7214126) q[1];
rz(-pi) q[2];
rz(-2.7338545) q[3];
sx q[3];
rz(-0.45555112) q[3];
sx q[3];
rz(1.3726485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.3624297) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841004) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(0.5350565) q[0];
rz(-1.1401945) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(-2.9761956) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1918614) q[0];
sx q[0];
rz(-2.9182069) q[0];
sx q[0];
rz(0.97747691) q[0];
x q[1];
rz(2.9781614) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(-0.93726678) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0796697) q[1];
sx q[1];
rz(-1.8975782) q[1];
sx q[1];
rz(2.5667739) q[1];
rz(-pi) q[2];
rz(-2.1256251) q[3];
sx q[3];
rz(-1.8053683) q[3];
sx q[3];
rz(-0.50333422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(-0.20544927) q[2];
rz(1.127634) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(-1.9873437) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-0.23194557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1922798) q[0];
sx q[0];
rz(-0.67040196) q[0];
sx q[0];
rz(0.85286661) q[0];
x q[1];
rz(-0.69552341) q[2];
sx q[2];
rz(-1.4624274) q[2];
sx q[2];
rz(-0.34851375) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0143118) q[1];
sx q[1];
rz(-0.31186549) q[1];
sx q[1];
rz(-0.77906268) q[1];
rz(-0.36851818) q[3];
sx q[3];
rz(-1.5889865) q[3];
sx q[3];
rz(2.8858678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0014687) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(-1.2403963) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0444788) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(0.99745497) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7339242) q[0];
sx q[0];
rz(-1.7051538) q[0];
sx q[0];
rz(2.2785447) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7315797) q[2];
sx q[2];
rz(-0.7747246) q[2];
sx q[2];
rz(2.0395525) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4221229) q[1];
sx q[1];
rz(-0.49233961) q[1];
sx q[1];
rz(-1.0455529) q[1];
x q[2];
rz(-0.44316767) q[3];
sx q[3];
rz(-2.245095) q[3];
sx q[3];
rz(2.5600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.22770195) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(-0.4513936) q[2];
rz(-2.732892) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(-1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.30615) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(2.5174482) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(0.61378941) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1873338) q[0];
sx q[0];
rz(-1.7566534) q[0];
sx q[0];
rz(-0.82877393) q[0];
rz(-pi) q[1];
rz(1.5067528) q[2];
sx q[2];
rz(-1.2876858) q[2];
sx q[2];
rz(-1.3751021) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.504096) q[1];
sx q[1];
rz(-1.6747961) q[1];
sx q[1];
rz(-1.5118447) q[1];
rz(3.0244163) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(-1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(1.1748574) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90010086) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(2.6532145) q[0];
rz(1.6237367) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(2.1571295) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41868608) q[0];
sx q[0];
rz(-1.5333999) q[0];
sx q[0];
rz(-2.5953369) q[0];
rz(-pi) q[1];
rz(2.7266399) q[2];
sx q[2];
rz(-2.3256362) q[2];
sx q[2];
rz(-0.50843898) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8967646) q[1];
sx q[1];
rz(-2.7762189) q[1];
sx q[1];
rz(-3.1246964) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6650388) q[3];
sx q[3];
rz(-1.2911951) q[3];
sx q[3];
rz(-2.9001146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70665923) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(2.4556659) q[0];
rz(0.39086875) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(2.2156782) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6305011) q[0];
sx q[0];
rz(-1.7112268) q[0];
sx q[0];
rz(1.3569843) q[0];
rz(-pi) q[1];
rz(0.7943031) q[2];
sx q[2];
rz(-1.0216733) q[2];
sx q[2];
rz(-1.6910451) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3013819) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(1.128003) q[1];
rz(-pi) q[2];
rz(0.52243201) q[3];
sx q[3];
rz(-2.9070832) q[3];
sx q[3];
rz(2.0731376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-0.53608981) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(1.6292054) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(-2.1283456) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7618338) q[0];
sx q[0];
rz(-1.338431) q[0];
sx q[0];
rz(0.14302111) q[0];
rz(-pi) q[1];
rz(1.6147862) q[2];
sx q[2];
rz(-1.9242052) q[2];
sx q[2];
rz(-2.8994438) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80876795) q[1];
sx q[1];
rz(-1.1507532) q[1];
sx q[1];
rz(1.4648449) q[1];
rz(2.3771044) q[3];
sx q[3];
rz(-1.0478684) q[3];
sx q[3];
rz(0.58512277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(2.3821793) q[2];
rz(1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(-2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(-0.57327523) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(-1.4105994) q[2];
sx q[2];
rz(-0.44993958) q[2];
sx q[2];
rz(-1.1100563) q[2];
rz(0.063354062) q[3];
sx q[3];
rz(-0.8739211) q[3];
sx q[3];
rz(-2.7779761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
