OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20237246) q[0];
sx q[0];
rz(-2.7352754) q[0];
sx q[0];
rz(2.321474) q[0];
rz(2.7804873) q[1];
sx q[1];
rz(-0.63280025) q[1];
sx q[1];
rz(-0.83067218) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75738534) q[0];
sx q[0];
rz(-1.2641347) q[0];
sx q[0];
rz(0.39461179) q[0];
x q[1];
rz(2.4029998) q[2];
sx q[2];
rz(-1.9361708) q[2];
sx q[2];
rz(1.5942758) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7100916) q[1];
sx q[1];
rz(-0.7845062) q[1];
sx q[1];
rz(1.1169408) q[1];
rz(-pi) q[2];
rz(-1.1675695) q[3];
sx q[3];
rz(-1.1551757) q[3];
sx q[3];
rz(2.2574539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43710199) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(-3.0214156) q[2];
rz(1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08081089) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(-2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(0.3266913) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.806843) q[0];
sx q[0];
rz(-0.93870367) q[0];
sx q[0];
rz(0.64081162) q[0];
rz(-1.8961043) q[2];
sx q[2];
rz(-2.6763958) q[2];
sx q[2];
rz(2.543769) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5527199) q[1];
sx q[1];
rz(-2.1761804) q[1];
sx q[1];
rz(2.6436716) q[1];
rz(-0.085751199) q[3];
sx q[3];
rz(-2.2567856) q[3];
sx q[3];
rz(0.12967295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6313173) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(-0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(-1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54863769) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(-0.27711162) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(1.057391) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7085416) q[0];
sx q[0];
rz(-1.9110702) q[0];
sx q[0];
rz(-3.1197381) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6234987) q[2];
sx q[2];
rz(-0.38803852) q[2];
sx q[2];
rz(-2.5990017) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96506572) q[1];
sx q[1];
rz(-2.2245193) q[1];
sx q[1];
rz(0.42018004) q[1];
rz(-1.3789165) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(-2.2172745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.7827591) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.7791629) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3574922) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(1.1401945) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(-0.16539703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7968984) q[0];
sx q[0];
rz(-1.3860774) q[0];
sx q[0];
rz(0.12634191) q[0];
rz(2.303316) q[2];
sx q[2];
rz(-0.24176134) q[2];
sx q[2];
rz(-1.4571112) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0619229) q[1];
sx q[1];
rz(-1.2440145) q[1];
sx q[1];
rz(2.5667739) q[1];
rz(-pi) q[2];
rz(1.1449279) q[3];
sx q[3];
rz(-2.5440359) q[3];
sx q[3];
rz(-1.715341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(-2.9361434) q[2];
rz(-2.0139587) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(-1.8849461) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797392) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(-0.35476312) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(-0.23194557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35271586) q[0];
sx q[0];
rz(-2.0577601) q[0];
sx q[0];
rz(-0.48077521) q[0];
x q[1];
rz(-1.7115713) q[2];
sx q[2];
rz(-2.2614334) q[2];
sx q[2];
rz(-1.8292793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31084727) q[1];
sx q[1];
rz(-1.7880882) q[1];
sx q[1];
rz(0.22549916) q[1];
x q[2];
rz(3.0911345) q[3];
sx q[3];
rz(-2.7726463) q[3];
sx q[3];
rz(-1.7794533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0014687) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(-1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0971138) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(-0.99745497) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31842445) q[0];
sx q[0];
rz(-0.71821763) q[0];
sx q[0];
rz(1.3658001) q[0];
x q[1];
rz(-2.4100254) q[2];
sx q[2];
rz(-1.288207) q[2];
sx q[2];
rz(2.3716795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5207386) q[1];
sx q[1];
rz(-1.3315017) q[1];
sx q[1];
rz(-1.1362856) q[1];
rz(-pi) q[2];
rz(-2.2951179) q[3];
sx q[3];
rz(-1.2292362) q[3];
sx q[3];
rz(-2.4404756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9138907) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(2.690199) q[2];
rz(-2.732892) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.30615) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(2.5278032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81543024) q[0];
sx q[0];
rz(-2.3809732) q[0];
sx q[0];
rz(-1.2994231) q[0];
rz(-2.8579312) q[2];
sx q[2];
rz(-1.6322871) q[2];
sx q[2];
rz(2.9279857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9871414) q[1];
sx q[1];
rz(-0.11949355) q[1];
sx q[1];
rz(0.51388545) q[1];
rz(3.0244163) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(1.5069435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.9667352) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(-1.6237367) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(-0.98446313) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7229066) q[0];
sx q[0];
rz(-1.6081928) q[0];
sx q[0];
rz(-0.54625578) q[0];
x q[1];
rz(-1.1659053) q[2];
sx q[2];
rz(-2.3003909) q[2];
sx q[2];
rz(2.0617495) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8967646) q[1];
sx q[1];
rz(-2.7762189) q[1];
sx q[1];
rz(-0.016896292) q[1];
x q[2];
rz(1.8833313) q[3];
sx q[3];
rz(-1.1141889) q[3];
sx q[3];
rz(1.4708335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(0.53517503) q[2];
rz(2.0914071) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64131367) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(-2.4556659) q[0];
rz(2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(2.2156782) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029322421) q[0];
sx q[0];
rz(-1.3591213) q[0];
sx q[0];
rz(0.14365833) q[0];
rz(2.4326153) q[2];
sx q[2];
rz(-2.211494) q[2];
sx q[2];
rz(-2.7880653) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8030501) q[1];
sx q[1];
rz(-1.1330789) q[1];
sx q[1];
rz(-0.16191698) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6191606) q[3];
sx q[3];
rz(-2.9070832) q[3];
sx q[3];
rz(-2.0731376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9459076) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(-2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560028) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(-1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(-2.1283456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2241867) q[0];
sx q[0];
rz(-1.4316443) q[0];
sx q[0];
rz(-1.3361206) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35372325) q[2];
sx q[2];
rz(-1.612066) q[2];
sx q[2];
rz(-1.313414) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0638748) q[1];
sx q[1];
rz(-0.43242726) q[1];
sx q[1];
rz(-0.23250154) q[1];
x q[2];
rz(-0.89684422) q[3];
sx q[3];
rz(-2.2138811) q[3];
sx q[3];
rz(0.53900063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(-0.75941336) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41291819) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(2.5683174) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(-1.1258833) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
rz(3.0782386) q[3];
sx q[3];
rz(-2.2676716) q[3];
sx q[3];
rz(0.36361658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
