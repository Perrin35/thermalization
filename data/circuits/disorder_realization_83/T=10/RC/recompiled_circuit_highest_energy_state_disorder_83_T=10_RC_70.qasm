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
rz(2.7163765) q[0];
sx q[0];
rz(-1.334231) q[0];
sx q[0];
rz(0.38036007) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(2.9549197) q[1];
sx q[1];
rz(8.5001707) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7036669) q[0];
sx q[0];
rz(-2.3406918) q[0];
sx q[0];
rz(-0.4842224) q[0];
rz(1.0546646) q[2];
sx q[2];
rz(-1.4251815) q[2];
sx q[2];
rz(0.21805412) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8262337) q[1];
sx q[1];
rz(-2.4118638) q[1];
sx q[1];
rz(-1.6297831) q[1];
rz(-pi) q[2];
rz(3.0136631) q[3];
sx q[3];
rz(-1.5797714) q[3];
sx q[3];
rz(0.25024807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1251462) q[2];
sx q[2];
rz(-2.1852198) q[2];
sx q[2];
rz(-2.1535786) q[2];
rz(2.9347349) q[3];
sx q[3];
rz(-1.4068973) q[3];
sx q[3];
rz(-0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4619231) q[0];
sx q[0];
rz(-1.4689057) q[0];
sx q[0];
rz(-2.8894506) q[0];
rz(-0.02154669) q[1];
sx q[1];
rz(-2.4485059) q[1];
sx q[1];
rz(-0.85539877) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48777521) q[0];
sx q[0];
rz(-0.027055351) q[0];
sx q[0];
rz(1.5144801) q[0];
rz(-pi) q[1];
rz(-2.5163745) q[2];
sx q[2];
rz(-2.4897235) q[2];
sx q[2];
rz(0.21833459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8586848) q[1];
sx q[1];
rz(-2.5577684) q[1];
sx q[1];
rz(1.1633384) q[1];
rz(-0.40327252) q[3];
sx q[3];
rz(-0.72411116) q[3];
sx q[3];
rz(1.4693174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1713193) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(-0.26101905) q[2];
rz(-0.039693443) q[3];
sx q[3];
rz(-1.3986162) q[3];
sx q[3];
rz(-2.5980914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4259341) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(-0.69450992) q[0];
rz(2.014324) q[1];
sx q[1];
rz(-2.290461) q[1];
sx q[1];
rz(0.99004254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58082579) q[0];
sx q[0];
rz(-0.35835727) q[0];
sx q[0];
rz(-0.23507765) q[0];
rz(-0.19845138) q[2];
sx q[2];
rz(-0.99764148) q[2];
sx q[2];
rz(-1.0400553) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5591919) q[1];
sx q[1];
rz(-1.9291996) q[1];
sx q[1];
rz(0.47689516) q[1];
rz(-pi) q[2];
rz(-0.95616266) q[3];
sx q[3];
rz(-0.80397881) q[3];
sx q[3];
rz(-1.1756736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29828829) q[2];
sx q[2];
rz(-2.1683606) q[2];
sx q[2];
rz(-2.6514371) q[2];
rz(-0.35689029) q[3];
sx q[3];
rz(-2.7481952) q[3];
sx q[3];
rz(-1.2662158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-3.0855584) q[0];
sx q[0];
rz(-1.9558676) q[0];
sx q[0];
rz(-2.2991142) q[0];
rz(0.8017686) q[1];
sx q[1];
rz(-0.27920488) q[1];
sx q[1];
rz(2.3559949) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88670896) q[0];
sx q[0];
rz(-2.2629316) q[0];
sx q[0];
rz(-1.9002302) q[0];
rz(-pi) q[1];
rz(1.4198205) q[2];
sx q[2];
rz(-2.2822126) q[2];
sx q[2];
rz(0.47455088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.848104) q[1];
sx q[1];
rz(-1.8791208) q[1];
sx q[1];
rz(2.5749194) q[1];
rz(-pi) q[2];
rz(-2.3219452) q[3];
sx q[3];
rz(-1.7213806) q[3];
sx q[3];
rz(-0.25170799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0164612) q[2];
sx q[2];
rz(-0.74840122) q[2];
sx q[2];
rz(1.008519) q[2];
rz(-0.85159167) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.1763024) q[0];
sx q[0];
rz(-0.98133123) q[0];
sx q[0];
rz(-2.0710131) q[0];
rz(2.3640682) q[1];
sx q[1];
rz(-0.91064015) q[1];
sx q[1];
rz(-2.1308965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7532446) q[0];
sx q[0];
rz(-2.0439612) q[0];
sx q[0];
rz(1.1654668) q[0];
x q[1];
rz(-1.2161184) q[2];
sx q[2];
rz(-0.18624072) q[2];
sx q[2];
rz(-2.8493488) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8292099) q[1];
sx q[1];
rz(-2.0264033) q[1];
sx q[1];
rz(-2.4036951) q[1];
rz(-pi) q[2];
rz(0.56212496) q[3];
sx q[3];
rz(-1.486612) q[3];
sx q[3];
rz(2.4474194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8396478) q[2];
sx q[2];
rz(-0.10493111) q[2];
sx q[2];
rz(1.5555752) q[2];
rz(1.0157061) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(-1.7293845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3509336) q[0];
sx q[0];
rz(-2.9582294) q[0];
sx q[0];
rz(2.3405128) q[0];
rz(-3.0084897) q[1];
sx q[1];
rz(-1.3214654) q[1];
sx q[1];
rz(0.4020234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19045705) q[0];
sx q[0];
rz(-1.2371089) q[0];
sx q[0];
rz(-0.78898375) q[0];
rz(0.0021462321) q[2];
sx q[2];
rz(-1.8615926) q[2];
sx q[2];
rz(-2.4418497) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.525156) q[1];
sx q[1];
rz(-1.0911988) q[1];
sx q[1];
rz(0.31071556) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4834255) q[3];
sx q[3];
rz(-1.6699381) q[3];
sx q[3];
rz(2.1226573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63783995) q[2];
sx q[2];
rz(-1.736234) q[2];
sx q[2];
rz(-2.3657738) q[2];
rz(1.6860298) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(-2.1337401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044947226) q[0];
sx q[0];
rz(-2.2438887) q[0];
sx q[0];
rz(0.67895472) q[0];
rz(-1.8567765) q[1];
sx q[1];
rz(-0.7130475) q[1];
sx q[1];
rz(1.3791893) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64039907) q[0];
sx q[0];
rz(-1.3974579) q[0];
sx q[0];
rz(-2.3480183) q[0];
rz(-pi) q[1];
rz(1.8502589) q[2];
sx q[2];
rz(-1.4905635) q[2];
sx q[2];
rz(-0.81278518) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1759291) q[1];
sx q[1];
rz(-1.129985) q[1];
sx q[1];
rz(0.55056527) q[1];
x q[2];
rz(-3.0857063) q[3];
sx q[3];
rz(-1.7706141) q[3];
sx q[3];
rz(-0.36147396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3369559) q[2];
sx q[2];
rz(-0.75504428) q[2];
sx q[2];
rz(-1.9026559) q[2];
rz(-1.8094481) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(-0.24470394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234583) q[0];
sx q[0];
rz(-1.1661538) q[0];
sx q[0];
rz(-0.13846692) q[0];
rz(2.9578517) q[1];
sx q[1];
rz(-0.67444363) q[1];
sx q[1];
rz(-2.8750681) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4999102) q[0];
sx q[0];
rz(-1.2453834) q[0];
sx q[0];
rz(2.9703559) q[0];
rz(-pi) q[1];
rz(-1.0037854) q[2];
sx q[2];
rz(-1.3409753) q[2];
sx q[2];
rz(-2.9858495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9723818) q[1];
sx q[1];
rz(-1.5145497) q[1];
sx q[1];
rz(-1.1509658) q[1];
x q[2];
rz(-2.4146904) q[3];
sx q[3];
rz(-1.57734) q[3];
sx q[3];
rz(0.18101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0515685) q[2];
sx q[2];
rz(-1.2568306) q[2];
sx q[2];
rz(0.23452342) q[2];
rz(-1.4759493) q[3];
sx q[3];
rz(-1.6022976) q[3];
sx q[3];
rz(-0.75638151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.804857) q[0];
sx q[0];
rz(-0.8661626) q[0];
sx q[0];
rz(2.3418703) q[0];
rz(2.5526478) q[1];
sx q[1];
rz(-0.84732333) q[1];
sx q[1];
rz(-2.934093) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070415592) q[0];
sx q[0];
rz(-1.9635597) q[0];
sx q[0];
rz(-1.429274) q[0];
x q[1];
rz(-2.8325808) q[2];
sx q[2];
rz(-2.2591619) q[2];
sx q[2];
rz(2.7613044) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.55940699) q[1];
sx q[1];
rz(-1.9684125) q[1];
sx q[1];
rz(2.29672) q[1];
rz(1.6274243) q[3];
sx q[3];
rz(-0.47041962) q[3];
sx q[3];
rz(-2.3868167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3346682) q[2];
sx q[2];
rz(-2.25756) q[2];
sx q[2];
rz(2.7743288) q[2];
rz(1.7043097) q[3];
sx q[3];
rz(-1.9673012) q[3];
sx q[3];
rz(-1.1737163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.2631898) q[0];
sx q[0];
rz(-2.1387565) q[0];
sx q[0];
rz(-0.81018418) q[0];
rz(2.0267678) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(-2.5883163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081757717) q[0];
sx q[0];
rz(-0.6404658) q[0];
sx q[0];
rz(-1.4778796) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3478339) q[2];
sx q[2];
rz(-2.0362034) q[2];
sx q[2];
rz(-2.2528354) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3685963) q[1];
sx q[1];
rz(-0.34075156) q[1];
sx q[1];
rz(0.98253886) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7891879) q[3];
sx q[3];
rz(-2.3239072) q[3];
sx q[3];
rz(-1.529983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1222003) q[2];
sx q[2];
rz(-2.4330752) q[2];
sx q[2];
rz(2.9019287) q[2];
rz(-1.3278809) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(-0.46752587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079530579) q[0];
sx q[0];
rz(-1.7876328) q[0];
sx q[0];
rz(-2.7476516) q[0];
rz(-2.486034) q[1];
sx q[1];
rz(-1.1759023) q[1];
sx q[1];
rz(-2.1942153) q[1];
rz(-1.2389567) q[2];
sx q[2];
rz(-1.9658026) q[2];
sx q[2];
rz(2.6323242) q[2];
rz(2.5611193) q[3];
sx q[3];
rz(-1.7844641) q[3];
sx q[3];
rz(-2.6337558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
