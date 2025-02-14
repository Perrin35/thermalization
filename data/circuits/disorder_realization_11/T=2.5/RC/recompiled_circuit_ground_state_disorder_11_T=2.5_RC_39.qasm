OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1144855) q[0];
sx q[0];
rz(3.0740102) q[0];
sx q[0];
rz(10.013784) q[0];
rz(2.0677805) q[1];
sx q[1];
rz(-1.685073) q[1];
sx q[1];
rz(-2.9340802) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1197101) q[0];
sx q[0];
rz(-1.3346905) q[0];
sx q[0];
rz(-0.2050928) q[0];
x q[1];
rz(-2.95822) q[2];
sx q[2];
rz(-1.4916821) q[2];
sx q[2];
rz(-0.26814207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3564975) q[1];
sx q[1];
rz(-1.5672088) q[1];
sx q[1];
rz(1.5990348) q[1];
x q[2];
rz(-0.21353586) q[3];
sx q[3];
rz(-2.6693404) q[3];
sx q[3];
rz(-0.47073281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.33240685) q[2];
sx q[2];
rz(-3.1250592) q[2];
sx q[2];
rz(-0.22745505) q[2];
rz(-2.4849232) q[3];
sx q[3];
rz(-0.69461099) q[3];
sx q[3];
rz(-0.53546661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.4559795) q[0];
sx q[0];
rz(-3.1094636) q[0];
sx q[0];
rz(-0.44422126) q[0];
rz(-1.1086858) q[1];
sx q[1];
rz(-1.3408835) q[1];
sx q[1];
rz(-1.1681555) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5614119) q[0];
sx q[0];
rz(-1.7636443) q[0];
sx q[0];
rz(0.69985244) q[0];
rz(-1.6158478) q[2];
sx q[2];
rz(-2.0241996) q[2];
sx q[2];
rz(-0.13870961) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8634974) q[1];
sx q[1];
rz(-1.4953164) q[1];
sx q[1];
rz(-1.1202034) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2222667) q[3];
sx q[3];
rz(-0.72685234) q[3];
sx q[3];
rz(-0.60692274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.68219677) q[2];
sx q[2];
rz(-2.0708059) q[2];
sx q[2];
rz(3.0219141) q[2];
rz(0.97673544) q[3];
sx q[3];
rz(-3.1141545) q[3];
sx q[3];
rz(-1.9500835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-1.479849) q[0];
sx q[0];
rz(-2.9763344) q[0];
sx q[0];
rz(1.4953493) q[0];
rz(-0.34548512) q[1];
sx q[1];
rz(-2.5472239) q[1];
sx q[1];
rz(2.3631309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0986106) q[0];
sx q[0];
rz(-1.5152644) q[0];
sx q[0];
rz(1.5776724) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2092152) q[2];
sx q[2];
rz(-0.86019963) q[2];
sx q[2];
rz(-1.4833409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.567019) q[1];
sx q[1];
rz(-1.4664093) q[1];
sx q[1];
rz(2.0135572) q[1];
rz(-pi) q[2];
rz(2.4069294) q[3];
sx q[3];
rz(-1.0473932) q[3];
sx q[3];
rz(1.911834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2848795) q[2];
sx q[2];
rz(-3.1046125) q[2];
sx q[2];
rz(-1.3899089) q[2];
rz(-3.1189647) q[3];
sx q[3];
rz(-2.8152864) q[3];
sx q[3];
rz(-2.7271395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6445899) q[0];
sx q[0];
rz(-0.039660064) q[0];
sx q[0];
rz(0.52763754) q[0];
rz(2.7854994) q[1];
sx q[1];
rz(-1.5047319) q[1];
sx q[1];
rz(1.5783232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8601807) q[0];
sx q[0];
rz(-1.1426272) q[0];
sx q[0];
rz(2.5606945) q[0];
rz(-pi) q[1];
rz(-3.061297) q[2];
sx q[2];
rz(-2.4622637) q[2];
sx q[2];
rz(-1.7905362) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8265022) q[1];
sx q[1];
rz(-2.0036526) q[1];
sx q[1];
rz(2.9587246) q[1];
rz(1.4107693) q[3];
sx q[3];
rz(-0.48086325) q[3];
sx q[3];
rz(-0.72070011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57033527) q[2];
sx q[2];
rz(-0.0084849914) q[2];
sx q[2];
rz(-0.045489475) q[2];
rz(2.3174543) q[3];
sx q[3];
rz(-2.6233311) q[3];
sx q[3];
rz(0.99388188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5387886) q[0];
sx q[0];
rz(-2.0863918) q[0];
sx q[0];
rz(1.7428727) q[0];
rz(2.973373) q[1];
sx q[1];
rz(-1.8054447) q[1];
sx q[1];
rz(-0.22163637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0960064) q[0];
sx q[0];
rz(-2.0915338) q[0];
sx q[0];
rz(-2.8686499) q[0];
x q[1];
rz(3.0732642) q[2];
sx q[2];
rz(-1.3573651) q[2];
sx q[2];
rz(2.1377856) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0757743) q[1];
sx q[1];
rz(-1.5464968) q[1];
sx q[1];
rz(1.1323117) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4517253) q[3];
sx q[3];
rz(-0.97056164) q[3];
sx q[3];
rz(-2.8056895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1571656) q[2];
sx q[2];
rz(-0.027313622) q[2];
sx q[2];
rz(-2.1402764) q[2];
rz(-2.2973513) q[3];
sx q[3];
rz(-3.0735569) q[3];
sx q[3];
rz(-0.74725738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35185128) q[0];
sx q[0];
rz(-2.8838367) q[0];
sx q[0];
rz(-1.3699654) q[0];
rz(2.9657956) q[1];
sx q[1];
rz(-1.628123) q[1];
sx q[1];
rz(1.0584077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31613126) q[0];
sx q[0];
rz(-2.2062345) q[0];
sx q[0];
rz(0.35346377) q[0];
rz(-1.5537698) q[2];
sx q[2];
rz(-2.5527708) q[2];
sx q[2];
rz(0.39521171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5762944) q[1];
sx q[1];
rz(-0.75913069) q[1];
sx q[1];
rz(-0.36328333) q[1];
x q[2];
rz(-0.36096548) q[3];
sx q[3];
rz(-2.1885893) q[3];
sx q[3];
rz(-1.82919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83219641) q[2];
sx q[2];
rz(-3.1377073) q[2];
sx q[2];
rz(-2.2957809) q[2];
rz(0.63515615) q[3];
sx q[3];
rz(-0.35717765) q[3];
sx q[3];
rz(-0.14491189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41088596) q[0];
sx q[0];
rz(-2.9697953) q[0];
sx q[0];
rz(-2.8806277) q[0];
rz(-1.7102526) q[1];
sx q[1];
rz(-0.15161082) q[1];
sx q[1];
rz(1.3014911) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1430014) q[0];
sx q[0];
rz(-0.81969417) q[0];
sx q[0];
rz(-0.9350594) q[0];
rz(3.030483) q[2];
sx q[2];
rz(-1.5859563) q[2];
sx q[2];
rz(2.3738101) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72922414) q[1];
sx q[1];
rz(-0.43392402) q[1];
sx q[1];
rz(0.093355066) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5742794) q[3];
sx q[3];
rz(-0.26393587) q[3];
sx q[3];
rz(1.4802295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5362376) q[2];
sx q[2];
rz(-1.321512) q[2];
sx q[2];
rz(3.1129254) q[2];
rz(0.043449314) q[3];
sx q[3];
rz(-2.9048007) q[3];
sx q[3];
rz(-1.4437599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0726149) q[0];
sx q[0];
rz(-0.34729877) q[0];
sx q[0];
rz(-0.29796991) q[0];
rz(-1.7295674) q[1];
sx q[1];
rz(-0.35930082) q[1];
sx q[1];
rz(1.5922458) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51247728) q[0];
sx q[0];
rz(-2.7393259) q[0];
sx q[0];
rz(0.8921807) q[0];
rz(-pi) q[1];
rz(1.9321597) q[2];
sx q[2];
rz(-3.123715) q[2];
sx q[2];
rz(1.0635384) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9142278) q[1];
sx q[1];
rz(-2.4691205) q[1];
sx q[1];
rz(-2.6537031) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9202254) q[3];
sx q[3];
rz(-1.2949756) q[3];
sx q[3];
rz(0.779895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90184244) q[2];
sx q[2];
rz(-0.032278927) q[2];
sx q[2];
rz(-0.65995222) q[2];
rz(-0.7800855) q[3];
sx q[3];
rz(-1.2939021) q[3];
sx q[3];
rz(1.9735146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6338947) q[0];
sx q[0];
rz(-0.036490353) q[0];
sx q[0];
rz(-2.3376035) q[0];
rz(-0.2313624) q[1];
sx q[1];
rz(-1.7295126) q[1];
sx q[1];
rz(3.0661327) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36614019) q[0];
sx q[0];
rz(-2.7260927) q[0];
sx q[0];
rz(0.52966161) q[0];
rz(-pi) q[1];
rz(-1.5701887) q[2];
sx q[2];
rz(-1.570833) q[2];
sx q[2];
rz(-1.4609877) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14743155) q[1];
sx q[1];
rz(-2.2975337) q[1];
sx q[1];
rz(-2.1925111) q[1];
rz(-pi) q[2];
rz(1.9583804) q[3];
sx q[3];
rz(-1.5389666) q[3];
sx q[3];
rz(1.4393161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1239473) q[2];
sx q[2];
rz(-0.46111527) q[2];
sx q[2];
rz(-0.41575113) q[2];
rz(1.4990643) q[3];
sx q[3];
rz(-3.1406904) q[3];
sx q[3];
rz(-0.7846964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7163664) q[0];
sx q[0];
rz(-3.0854736) q[0];
sx q[0];
rz(0.29384336) q[0];
rz(1.6055239) q[1];
sx q[1];
rz(-0.87296456) q[1];
sx q[1];
rz(-1.6764838) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7820202) q[0];
sx q[0];
rz(-1.5609419) q[0];
sx q[0];
rz(0.1043679) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56581409) q[2];
sx q[2];
rz(-2.4560258) q[2];
sx q[2];
rz(-3.1367347) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6584545) q[1];
sx q[1];
rz(-0.92441974) q[1];
sx q[1];
rz(-1.5366244) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5663537) q[3];
sx q[3];
rz(-1.5779788) q[3];
sx q[3];
rz(1.3346163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6426223) q[2];
sx q[2];
rz(-0.62752807) q[2];
sx q[2];
rz(1.9210531) q[2];
rz(-0.39310655) q[3];
sx q[3];
rz(-0.016963907) q[3];
sx q[3];
rz(0.19256798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4569693) q[0];
sx q[0];
rz(-1.4263117) q[0];
sx q[0];
rz(-0.34761467) q[0];
rz(2.5034703) q[1];
sx q[1];
rz(-2.3633524) q[1];
sx q[1];
rz(-0.6755158) q[1];
rz(0.38370321) q[2];
sx q[2];
rz(-2.9331238) q[2];
sx q[2];
rz(2.5274656) q[2];
rz(1.5777231) q[3];
sx q[3];
rz(-1.6835811) q[3];
sx q[3];
rz(-1.2006105) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
