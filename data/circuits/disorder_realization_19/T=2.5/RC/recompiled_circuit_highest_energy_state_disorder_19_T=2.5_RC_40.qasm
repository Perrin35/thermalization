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
rz(2.1228696) q[0];
sx q[0];
rz(-2.2824204) q[0];
sx q[0];
rz(-2.3265042) q[0];
rz(1.9563142) q[1];
sx q[1];
rz(4.5524608) q[1];
sx q[1];
rz(8.3571385) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6237804) q[0];
sx q[0];
rz(-1.5945487) q[0];
sx q[0];
rz(-1.1393113) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63500603) q[2];
sx q[2];
rz(-0.62261144) q[2];
sx q[2];
rz(-0.16281637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8274535) q[1];
sx q[1];
rz(-1.8874536) q[1];
sx q[1];
rz(-0.10152557) q[1];
x q[2];
rz(1.2514011) q[3];
sx q[3];
rz(-2.4301964) q[3];
sx q[3];
rz(0.2371725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.69433576) q[2];
sx q[2];
rz(-0.7434291) q[2];
sx q[2];
rz(1.8667963) q[2];
rz(-0.46191195) q[3];
sx q[3];
rz(-0.67449823) q[3];
sx q[3];
rz(-1.1326724) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79134113) q[0];
sx q[0];
rz(-0.30838648) q[0];
sx q[0];
rz(1.8885008) q[0];
rz(-0.14532267) q[1];
sx q[1];
rz(-1.3959613) q[1];
sx q[1];
rz(-1.0911509) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4900794) q[0];
sx q[0];
rz(-1.3567748) q[0];
sx q[0];
rz(-1.1351311) q[0];
rz(-0.65204377) q[2];
sx q[2];
rz(-0.55527675) q[2];
sx q[2];
rz(-0.20737831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5877643) q[1];
sx q[1];
rz(-1.7436281) q[1];
sx q[1];
rz(0.21765222) q[1];
rz(-1.9698148) q[3];
sx q[3];
rz(-1.1091091) q[3];
sx q[3];
rz(-2.6103013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6609409) q[2];
sx q[2];
rz(-1.7512243) q[2];
sx q[2];
rz(0.52978984) q[2];
rz(-2.3482813) q[3];
sx q[3];
rz(-1.6042234) q[3];
sx q[3];
rz(0.62354273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6563501) q[0];
sx q[0];
rz(-0.95004496) q[0];
sx q[0];
rz(-0.52870885) q[0];
rz(-0.54620019) q[1];
sx q[1];
rz(-0.9587973) q[1];
sx q[1];
rz(0.34034696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84547323) q[0];
sx q[0];
rz(-1.8438135) q[0];
sx q[0];
rz(0.34015981) q[0];
x q[1];
rz(1.7163926) q[2];
sx q[2];
rz(-1.3383838) q[2];
sx q[2];
rz(0.53047859) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9941147) q[1];
sx q[1];
rz(-3.0469739) q[1];
sx q[1];
rz(3.1069744) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82967088) q[3];
sx q[3];
rz(-1.8541012) q[3];
sx q[3];
rz(-1.0782575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9693552) q[2];
sx q[2];
rz(-0.41219553) q[2];
sx q[2];
rz(-0.42984143) q[2];
rz(2.3853081) q[3];
sx q[3];
rz(-2.9679306) q[3];
sx q[3];
rz(-2.3156796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7569358) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(-1.4935619) q[0];
rz(0.76796302) q[1];
sx q[1];
rz(-2.6710644) q[1];
sx q[1];
rz(0.54642645) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1892197) q[0];
sx q[0];
rz(-1.6354939) q[0];
sx q[0];
rz(1.6420664) q[0];
rz(-pi) q[1];
rz(-2.619952) q[2];
sx q[2];
rz(-1.120479) q[2];
sx q[2];
rz(-2.3797824) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1122894) q[1];
sx q[1];
rz(-1.422907) q[1];
sx q[1];
rz(-1.6690955) q[1];
rz(-pi) q[2];
rz(-0.91937842) q[3];
sx q[3];
rz(-0.99777824) q[3];
sx q[3];
rz(-2.1185377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7118608) q[2];
sx q[2];
rz(-1.4873361) q[2];
sx q[2];
rz(-2.8374953) q[2];
rz(-1.7763304) q[3];
sx q[3];
rz(-1.920776) q[3];
sx q[3];
rz(-1.0767267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6292608) q[0];
sx q[0];
rz(-2.6973695) q[0];
sx q[0];
rz(-3.1410134) q[0];
rz(2.0594788) q[1];
sx q[1];
rz(-2.6172456) q[1];
sx q[1];
rz(-2.2023315) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2321286) q[0];
sx q[0];
rz(-2.3298023) q[0];
sx q[0];
rz(-0.62007298) q[0];
rz(-0.95920697) q[2];
sx q[2];
rz(-1.8610753) q[2];
sx q[2];
rz(2.2242351) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68785948) q[1];
sx q[1];
rz(-0.61111585) q[1];
sx q[1];
rz(1.9495717) q[1];
rz(2.9143798) q[3];
sx q[3];
rz(-1.876653) q[3];
sx q[3];
rz(1.3671041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1028563) q[2];
sx q[2];
rz(-1.880371) q[2];
sx q[2];
rz(2.8803414) q[2];
rz(2.2767565) q[3];
sx q[3];
rz(-1.4120925) q[3];
sx q[3];
rz(-2.849546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.29702) q[0];
sx q[0];
rz(-1.427587) q[0];
sx q[0];
rz(-0.20508668) q[0];
rz(1.0467485) q[1];
sx q[1];
rz(-1.7457242) q[1];
sx q[1];
rz(2.7632025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32961938) q[0];
sx q[0];
rz(-1.1866335) q[0];
sx q[0];
rz(0.20150082) q[0];
rz(2.66983) q[2];
sx q[2];
rz(-1.4233575) q[2];
sx q[2];
rz(-2.9094537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4503895) q[1];
sx q[1];
rz(-2.4892163) q[1];
sx q[1];
rz(2.1689586) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69244416) q[3];
sx q[3];
rz(-1.2298349) q[3];
sx q[3];
rz(2.8657262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4947074) q[2];
sx q[2];
rz(-1.9834221) q[2];
sx q[2];
rz(-1.0542487) q[2];
rz(0.65822893) q[3];
sx q[3];
rz(-2.4963278) q[3];
sx q[3];
rz(2.6575507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14195104) q[0];
sx q[0];
rz(-1.9331837) q[0];
sx q[0];
rz(-2.5010338) q[0];
rz(0.41796747) q[1];
sx q[1];
rz(-1.2203981) q[1];
sx q[1];
rz(-2.3366065) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98950878) q[0];
sx q[0];
rz(-1.5474657) q[0];
sx q[0];
rz(-0.70744608) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1518794) q[2];
sx q[2];
rz(-2.3275314) q[2];
sx q[2];
rz(2.4474622) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5857081) q[1];
sx q[1];
rz(-1.5991365) q[1];
sx q[1];
rz(3.1217087) q[1];
rz(-pi) q[2];
rz(0.59155699) q[3];
sx q[3];
rz(-1.322896) q[3];
sx q[3];
rz(1.3755368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54887041) q[2];
sx q[2];
rz(-0.9950811) q[2];
sx q[2];
rz(2.3804046) q[2];
rz(0.74782863) q[3];
sx q[3];
rz(-1.8239832) q[3];
sx q[3];
rz(0.67659155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6250896) q[0];
sx q[0];
rz(-0.28446063) q[0];
sx q[0];
rz(1.6647343) q[0];
rz(-0.45627108) q[1];
sx q[1];
rz(-1.7218593) q[1];
sx q[1];
rz(-2.2705073) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81646279) q[0];
sx q[0];
rz(-1.5832381) q[0];
sx q[0];
rz(-0.61412707) q[0];
x q[1];
rz(-1.4159059) q[2];
sx q[2];
rz(-1.911507) q[2];
sx q[2];
rz(2.129385) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91034094) q[1];
sx q[1];
rz(-1.7897072) q[1];
sx q[1];
rz(0.33556767) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1444451) q[3];
sx q[3];
rz(-1.6255857) q[3];
sx q[3];
rz(2.3940115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90290922) q[2];
sx q[2];
rz(-0.52672714) q[2];
sx q[2];
rz(0.61863679) q[2];
rz(-3.1213308) q[3];
sx q[3];
rz(-0.94347763) q[3];
sx q[3];
rz(-0.77997911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063754931) q[0];
sx q[0];
rz(-1.5564593) q[0];
sx q[0];
rz(-2.7177287) q[0];
rz(-0.18691143) q[1];
sx q[1];
rz(-0.82004768) q[1];
sx q[1];
rz(1.6835469) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8992726) q[0];
sx q[0];
rz(-1.5849216) q[0];
sx q[0];
rz(0.11798162) q[0];
rz(-2.0293268) q[2];
sx q[2];
rz(-0.56156172) q[2];
sx q[2];
rz(-2.4379345) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9002237) q[1];
sx q[1];
rz(-1.4321064) q[1];
sx q[1];
rz(1.5189511) q[1];
x q[2];
rz(-0.23641674) q[3];
sx q[3];
rz(-2.7763753) q[3];
sx q[3];
rz(-1.055298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75108782) q[2];
sx q[2];
rz(-1.0672528) q[2];
sx q[2];
rz(-1.0941774) q[2];
rz(-1.68082) q[3];
sx q[3];
rz(-1.3760309) q[3];
sx q[3];
rz(-1.338753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8193034) q[0];
sx q[0];
rz(-1.7604473) q[0];
sx q[0];
rz(-1.5018916) q[0];
rz(-0.27885258) q[1];
sx q[1];
rz(-2.1809705) q[1];
sx q[1];
rz(-1.0241114) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0251966) q[0];
sx q[0];
rz(-1.0828583) q[0];
sx q[0];
rz(-1.4665718) q[0];
x q[1];
rz(1.7488519) q[2];
sx q[2];
rz(-2.8626056) q[2];
sx q[2];
rz(1.6131608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1361724) q[1];
sx q[1];
rz(-1.3498303) q[1];
sx q[1];
rz(1.3653838) q[1];
rz(-2.315178) q[3];
sx q[3];
rz(-1.8262939) q[3];
sx q[3];
rz(0.57533732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9378822) q[2];
sx q[2];
rz(-1.2869765) q[2];
sx q[2];
rz(-0.53696519) q[2];
rz(-1.9133866) q[3];
sx q[3];
rz(-0.95913404) q[3];
sx q[3];
rz(0.29135191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1163597) q[0];
sx q[0];
rz(-1.3675084) q[0];
sx q[0];
rz(1.0687923) q[0];
rz(0.7069201) q[1];
sx q[1];
rz(-2.0126577) q[1];
sx q[1];
rz(-0.001002034) q[1];
rz(2.7189485) q[2];
sx q[2];
rz(-0.44036897) q[2];
sx q[2];
rz(-1.8457495) q[2];
rz(1.3523921) q[3];
sx q[3];
rz(-1.7159749) q[3];
sx q[3];
rz(0.80930474) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
