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
rz(1.3820833) q[0];
sx q[0];
rz(-0.74892646) q[0];
sx q[0];
rz(1.7480667) q[0];
rz(0.8028318) q[1];
sx q[1];
rz(-0.79610151) q[1];
sx q[1];
rz(0.5308477) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33541629) q[0];
sx q[0];
rz(-2.5160976) q[0];
sx q[0];
rz(-0.92910398) q[0];
rz(-1.5535001) q[2];
sx q[2];
rz(-1.4040385) q[2];
sx q[2];
rz(-1.1823428) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4179791) q[1];
sx q[1];
rz(-2.6043677) q[1];
sx q[1];
rz(1.4239271) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10194709) q[3];
sx q[3];
rz(-1.0539712) q[3];
sx q[3];
rz(2.068813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.406245) q[2];
sx q[2];
rz(-2.2082059) q[2];
sx q[2];
rz(2.1437342) q[2];
rz(2.2570299) q[3];
sx q[3];
rz(-1.9622784) q[3];
sx q[3];
rz(-2.4295889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74580055) q[0];
sx q[0];
rz(-0.53676787) q[0];
sx q[0];
rz(-1.0774379) q[0];
rz(2.5224345) q[1];
sx q[1];
rz(-1.8749323) q[1];
sx q[1];
rz(-1.9177297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2015468) q[0];
sx q[0];
rz(-2.8000961) q[0];
sx q[0];
rz(-1.0251371) q[0];
rz(3.139758) q[2];
sx q[2];
rz(-1.6009838) q[2];
sx q[2];
rz(-1.0890397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.465816) q[1];
sx q[1];
rz(-0.73938939) q[1];
sx q[1];
rz(-0.63111102) q[1];
rz(-pi) q[2];
rz(2.4183351) q[3];
sx q[3];
rz(-0.92953909) q[3];
sx q[3];
rz(2.0632405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9549442) q[2];
sx q[2];
rz(-0.84543219) q[2];
sx q[2];
rz(2.1050982) q[2];
rz(1.9519818) q[3];
sx q[3];
rz(-1.2396953) q[3];
sx q[3];
rz(0.069124393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67967296) q[0];
sx q[0];
rz(-2.9000977) q[0];
sx q[0];
rz(-1.4659721) q[0];
rz(1.2566603) q[1];
sx q[1];
rz(-0.64777056) q[1];
sx q[1];
rz(-0.57674903) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7490596) q[0];
sx q[0];
rz(-2.1196094) q[0];
sx q[0];
rz(-0.69290036) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45879634) q[2];
sx q[2];
rz(-1.5374628) q[2];
sx q[2];
rz(-1.3622487) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27611342) q[1];
sx q[1];
rz(-0.41655585) q[1];
sx q[1];
rz(1.4426484) q[1];
rz(3.0528487) q[3];
sx q[3];
rz(-1.7341494) q[3];
sx q[3];
rz(0.08069144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3999148) q[2];
sx q[2];
rz(-2.8027813) q[2];
sx q[2];
rz(2.3728288) q[2];
rz(0.1712884) q[3];
sx q[3];
rz(-1.4408709) q[3];
sx q[3];
rz(2.7870074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6866368) q[0];
sx q[0];
rz(-0.88382116) q[0];
sx q[0];
rz(0.53632847) q[0];
rz(2.3388011) q[1];
sx q[1];
rz(-1.2840459) q[1];
sx q[1];
rz(0.098043052) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8756384) q[0];
sx q[0];
rz(-1.2438626) q[0];
sx q[0];
rz(0.083766706) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78196193) q[2];
sx q[2];
rz(-1.1267917) q[2];
sx q[2];
rz(0.16829106) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9598011) q[1];
sx q[1];
rz(-1.0611294) q[1];
sx q[1];
rz(0.23284973) q[1];
x q[2];
rz(2.732787) q[3];
sx q[3];
rz(-1.9681566) q[3];
sx q[3];
rz(-1.8026343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1596277) q[2];
sx q[2];
rz(-2.1106796) q[2];
sx q[2];
rz(-0.31943303) q[2];
rz(2.5751513) q[3];
sx q[3];
rz(-2.3947377) q[3];
sx q[3];
rz(-1.1613065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40120688) q[0];
sx q[0];
rz(-1.3196608) q[0];
sx q[0];
rz(-2.5886986) q[0];
rz(2.3623908) q[1];
sx q[1];
rz(-0.82256493) q[1];
sx q[1];
rz(-0.78831569) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8126683) q[0];
sx q[0];
rz(-1.5018437) q[0];
sx q[0];
rz(-2.8806995) q[0];
x q[1];
rz(-2.9328038) q[2];
sx q[2];
rz(-1.5940856) q[2];
sx q[2];
rz(-1.4643027) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31438875) q[1];
sx q[1];
rz(-2.6174859) q[1];
sx q[1];
rz(-0.57669611) q[1];
rz(0.29946976) q[3];
sx q[3];
rz(-1.3336998) q[3];
sx q[3];
rz(0.32215986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7439482) q[2];
sx q[2];
rz(-1.1279736) q[2];
sx q[2];
rz(2.9359342) q[2];
rz(-1.7913943) q[3];
sx q[3];
rz(-0.74068991) q[3];
sx q[3];
rz(-3.1309483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81269294) q[0];
sx q[0];
rz(-2.6943272) q[0];
sx q[0];
rz(1.6269667) q[0];
rz(-2.336592) q[1];
sx q[1];
rz(-1.4470419) q[1];
sx q[1];
rz(2.8256493) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7132063) q[0];
sx q[0];
rz(-1.180219) q[0];
sx q[0];
rz(-1.9414563) q[0];
rz(-1.9419901) q[2];
sx q[2];
rz(-0.99200574) q[2];
sx q[2];
rz(0.82497342) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1263663) q[1];
sx q[1];
rz(-2.7607007) q[1];
sx q[1];
rz(0.81896255) q[1];
x q[2];
rz(0.31031761) q[3];
sx q[3];
rz(-1.7162706) q[3];
sx q[3];
rz(-0.014257243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.463795) q[2];
sx q[2];
rz(-2.4424876) q[2];
sx q[2];
rz(-0.15764906) q[2];
rz(-0.53358233) q[3];
sx q[3];
rz(-1.491051) q[3];
sx q[3];
rz(1.5177479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59876281) q[0];
sx q[0];
rz(-2.6451126) q[0];
sx q[0];
rz(0.98094034) q[0];
rz(-0.44278231) q[1];
sx q[1];
rz(-1.9552224) q[1];
sx q[1];
rz(2.669899) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1991292) q[0];
sx q[0];
rz(-0.53265306) q[0];
sx q[0];
rz(-2.1739302) q[0];
rz(-pi) q[1];
rz(-1.4002789) q[2];
sx q[2];
rz(-1.8461175) q[2];
sx q[2];
rz(-1.4845276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.065475796) q[1];
sx q[1];
rz(-3.0055176) q[1];
sx q[1];
rz(-3.0398366) q[1];
rz(-pi) q[2];
rz(-2.6568233) q[3];
sx q[3];
rz(-2.7607684) q[3];
sx q[3];
rz(0.72431475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.20785759) q[2];
sx q[2];
rz(-0.8693049) q[2];
sx q[2];
rz(3.0233439) q[2];
rz(-2.5737428) q[3];
sx q[3];
rz(-1.7852781) q[3];
sx q[3];
rz(-2.3387199) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0589013) q[0];
sx q[0];
rz(-1.7401798) q[0];
sx q[0];
rz(-0.68728224) q[0];
rz(1.8742689) q[1];
sx q[1];
rz(-1.9272389) q[1];
sx q[1];
rz(-0.6257239) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83118781) q[0];
sx q[0];
rz(-3.1200231) q[0];
sx q[0];
rz(1.789272) q[0];
x q[1];
rz(-0.86144336) q[2];
sx q[2];
rz(-0.77327912) q[2];
sx q[2];
rz(2.803162) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3290112) q[1];
sx q[1];
rz(-2.1454403) q[1];
sx q[1];
rz(-3.126815) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1504907) q[3];
sx q[3];
rz(-1.8324346) q[3];
sx q[3];
rz(-1.8707448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.22831336) q[2];
sx q[2];
rz(-1.2484756) q[2];
sx q[2];
rz(1.0949562) q[2];
rz(2.6878808) q[3];
sx q[3];
rz(-0.79115051) q[3];
sx q[3];
rz(-2.9554844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5472645) q[0];
sx q[0];
rz(-0.493395) q[0];
sx q[0];
rz(0.067597978) q[0];
rz(-1.0293845) q[1];
sx q[1];
rz(-1.8600347) q[1];
sx q[1];
rz(-3.0060815) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96582039) q[0];
sx q[0];
rz(-2.1012573) q[0];
sx q[0];
rz(2.4932753) q[0];
rz(-pi) q[1];
rz(2.4004647) q[2];
sx q[2];
rz(-1.6911881) q[2];
sx q[2];
rz(-0.84726221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77733946) q[1];
sx q[1];
rz(-2.4852264) q[1];
sx q[1];
rz(-2.6692764) q[1];
rz(-pi) q[2];
rz(-0.71262757) q[3];
sx q[3];
rz(-0.36667675) q[3];
sx q[3];
rz(0.7798051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6873261) q[2];
sx q[2];
rz(-0.25298515) q[2];
sx q[2];
rz(2.6542286) q[2];
rz(0.38960114) q[3];
sx q[3];
rz(-0.95249683) q[3];
sx q[3];
rz(0.065936955) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366632) q[0];
sx q[0];
rz(-1.0856029) q[0];
sx q[0];
rz(1.3847463) q[0];
rz(2.0866277) q[1];
sx q[1];
rz(-2.2672548) q[1];
sx q[1];
rz(-0.25996444) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6523165) q[0];
sx q[0];
rz(-0.92149177) q[0];
sx q[0];
rz(-0.32353521) q[0];
rz(-pi) q[1];
rz(2.2222338) q[2];
sx q[2];
rz(-2.9769889) q[2];
sx q[2];
rz(-2.2155264) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7284026) q[1];
sx q[1];
rz(-1.5356881) q[1];
sx q[1];
rz(2.1988695) q[1];
rz(1.7694951) q[3];
sx q[3];
rz(-1.3532146) q[3];
sx q[3];
rz(-2.8189903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7197623) q[2];
sx q[2];
rz(-2.7037342) q[2];
sx q[2];
rz(-1.7913294) q[2];
rz(2.0566025) q[3];
sx q[3];
rz(-1.9271873) q[3];
sx q[3];
rz(-2.0124281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.0706901) q[0];
sx q[0];
rz(-0.68690837) q[0];
sx q[0];
rz(2.6223781) q[0];
rz(1.6593973) q[1];
sx q[1];
rz(-1.8435602) q[1];
sx q[1];
rz(-1.6866121) q[1];
rz(2.5605911) q[2];
sx q[2];
rz(-1.1270564) q[2];
sx q[2];
rz(-1.8869274) q[2];
rz(-1.0709892) q[3];
sx q[3];
rz(-1.3060112) q[3];
sx q[3];
rz(-1.3597091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
