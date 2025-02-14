OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.423288) q[0];
sx q[0];
rz(-2.365132) q[0];
sx q[0];
rz(-2.7756696) q[0];
rz(-2.794682) q[1];
sx q[1];
rz(3.4263098) q[1];
sx q[1];
rz(10.813536) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66807085) q[0];
sx q[0];
rz(-1.6374303) q[0];
sx q[0];
rz(0.20247831) q[0];
rz(0.69268815) q[2];
sx q[2];
rz(-1.8028304) q[2];
sx q[2];
rz(1.2218066) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5386913) q[1];
sx q[1];
rz(-0.60990342) q[1];
sx q[1];
rz(-0.84714386) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4284989) q[3];
sx q[3];
rz(-2.4131107) q[3];
sx q[3];
rz(-0.99458867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1846788) q[2];
sx q[2];
rz(-0.92014402) q[2];
sx q[2];
rz(-2.206395) q[2];
rz(-0.30185559) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(-2.7479318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(-2.1644156) q[0];
sx q[0];
rz(-1.2657607) q[0];
sx q[0];
rz(-2.6432977) q[0];
rz(1.7138819) q[1];
sx q[1];
rz(-1.9352813) q[1];
sx q[1];
rz(2.0548342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3918864) q[0];
sx q[0];
rz(-1.8072309) q[0];
sx q[0];
rz(-2.4357585) q[0];
rz(2.1939414) q[2];
sx q[2];
rz(-0.32440475) q[2];
sx q[2];
rz(2.0632921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47269184) q[1];
sx q[1];
rz(-1.7311061) q[1];
sx q[1];
rz(0.44071381) q[1];
rz(-pi) q[2];
rz(3.1051136) q[3];
sx q[3];
rz(-2.8231797) q[3];
sx q[3];
rz(-0.94389254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.19865092) q[2];
sx q[2];
rz(-2.8309839) q[2];
sx q[2];
rz(-1.3322213) q[2];
rz(1.1474991) q[3];
sx q[3];
rz(-1.56286) q[3];
sx q[3];
rz(1.3326299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9320817) q[0];
sx q[0];
rz(-1.4028343) q[0];
sx q[0];
rz(2.984356) q[0];
rz(1.2278185) q[1];
sx q[1];
rz(-2.9324052) q[1];
sx q[1];
rz(-0.42207178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76772375) q[0];
sx q[0];
rz(-1.6607303) q[0];
sx q[0];
rz(-1.7492245) q[0];
x q[1];
rz(-2.8639925) q[2];
sx q[2];
rz(-1.7381845) q[2];
sx q[2];
rz(-0.041942747) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2873309) q[1];
sx q[1];
rz(-0.39188436) q[1];
sx q[1];
rz(-1.803483) q[1];
rz(-pi) q[2];
rz(-1.8881599) q[3];
sx q[3];
rz(-0.79447132) q[3];
sx q[3];
rz(-1.0004071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.32855365) q[2];
sx q[2];
rz(-0.96031323) q[2];
sx q[2];
rz(-0.73406827) q[2];
rz(-1.0775393) q[3];
sx q[3];
rz(-2.4534093) q[3];
sx q[3];
rz(0.075210007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5828736) q[0];
sx q[0];
rz(-2.131077) q[0];
sx q[0];
rz(0.016481312) q[0];
rz(3.0670498) q[1];
sx q[1];
rz(-2.6838979) q[1];
sx q[1];
rz(0.25513908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82127386) q[0];
sx q[0];
rz(-0.50423586) q[0];
sx q[0];
rz(-2.8404499) q[0];
rz(-pi) q[1];
rz(-2.8066977) q[2];
sx q[2];
rz(-1.0303632) q[2];
sx q[2];
rz(-1.188736) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8447579) q[1];
sx q[1];
rz(-2.1898309) q[1];
sx q[1];
rz(-2.1059787) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11394088) q[3];
sx q[3];
rz(-1.2234395) q[3];
sx q[3];
rz(-1.8983656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5269346) q[2];
sx q[2];
rz(-0.69710985) q[2];
sx q[2];
rz(-2.441791) q[2];
rz(-3.0497293) q[3];
sx q[3];
rz(-1.2173165) q[3];
sx q[3];
rz(-2.6868668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62234539) q[0];
sx q[0];
rz(-0.82038251) q[0];
sx q[0];
rz(-2.0297594) q[0];
rz(2.9513997) q[1];
sx q[1];
rz(-0.88514248) q[1];
sx q[1];
rz(-1.1598738) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5525813) q[0];
sx q[0];
rz(-1.9144626) q[0];
sx q[0];
rz(3.1049743) q[0];
x q[1];
rz(2.2091523) q[2];
sx q[2];
rz(-1.6886504) q[2];
sx q[2];
rz(-1.7913851) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.24007356) q[1];
sx q[1];
rz(-2.1131971) q[1];
sx q[1];
rz(-1.8611004) q[1];
rz(-pi) q[2];
rz(-2.0926863) q[3];
sx q[3];
rz(-1.5318222) q[3];
sx q[3];
rz(-0.029595395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5708892) q[2];
sx q[2];
rz(-2.3184226) q[2];
sx q[2];
rz(-0.23045753) q[2];
rz(2.0795836) q[3];
sx q[3];
rz(-1.8657203) q[3];
sx q[3];
rz(-1.6331858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.688711) q[0];
sx q[0];
rz(-2.0827561) q[0];
sx q[0];
rz(-1.7857312) q[0];
rz(-2.611825) q[1];
sx q[1];
rz(-0.7822839) q[1];
sx q[1];
rz(-1.1518325) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065581948) q[0];
sx q[0];
rz(-1.9261203) q[0];
sx q[0];
rz(1.4655345) q[0];
rz(-0.095345796) q[2];
sx q[2];
rz(-1.0311677) q[2];
sx q[2];
rz(-2.6971872) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7499654) q[1];
sx q[1];
rz(-0.4164975) q[1];
sx q[1];
rz(-2.6451664) q[1];
rz(-pi) q[2];
rz(-1.4532667) q[3];
sx q[3];
rz(-1.3580048) q[3];
sx q[3];
rz(0.26373395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4749703) q[2];
sx q[2];
rz(-0.64133659) q[2];
sx q[2];
rz(-0.46710157) q[2];
rz(-2.1304255) q[3];
sx q[3];
rz(-2.7979388) q[3];
sx q[3];
rz(1.7721734) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2924627) q[0];
sx q[0];
rz(-2.9865773) q[0];
sx q[0];
rz(2.9169061) q[0];
rz(-2.977773) q[1];
sx q[1];
rz(-0.33912173) q[1];
sx q[1];
rz(2.4952369) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77262473) q[0];
sx q[0];
rz(-0.59212055) q[0];
sx q[0];
rz(2.0298241) q[0];
rz(-pi) q[1];
rz(-0.21804131) q[2];
sx q[2];
rz(-1.8865117) q[2];
sx q[2];
rz(1.9670847) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4187635) q[1];
sx q[1];
rz(-1.495528) q[1];
sx q[1];
rz(-2.597517) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4348174) q[3];
sx q[3];
rz(-1.1866015) q[3];
sx q[3];
rz(-0.7766436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9859163) q[2];
sx q[2];
rz(-1.6935657) q[2];
sx q[2];
rz(-0.4099561) q[2];
rz(2.3112467) q[3];
sx q[3];
rz(-2.1928619) q[3];
sx q[3];
rz(-1.4478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53772563) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(0.24630462) q[0];
rz(-0.82659563) q[1];
sx q[1];
rz(-1.7431424) q[1];
sx q[1];
rz(-2.8776317) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9519099) q[0];
sx q[0];
rz(-0.095746843) q[0];
sx q[0];
rz(-1.4950947) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3578431) q[2];
sx q[2];
rz(-1.2211868) q[2];
sx q[2];
rz(-2.8658681) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3726793) q[1];
sx q[1];
rz(-0.9039592) q[1];
sx q[1];
rz(-0.36075488) q[1];
x q[2];
rz(1.1353536) q[3];
sx q[3];
rz(-2.1135984) q[3];
sx q[3];
rz(-0.28771985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0366514) q[2];
sx q[2];
rz(-1.2678601) q[2];
sx q[2];
rz(0.71511739) q[2];
rz(-3.0702843) q[3];
sx q[3];
rz(-1.5569867) q[3];
sx q[3];
rz(-2.3947072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1320553) q[0];
sx q[0];
rz(-1.4893463) q[0];
sx q[0];
rz(0.43553964) q[0];
rz(-2.9938193) q[1];
sx q[1];
rz(-1.4316033) q[1];
sx q[1];
rz(-1.4685644) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72318059) q[0];
sx q[0];
rz(-1.3910798) q[0];
sx q[0];
rz(-0.43850684) q[0];
rz(-pi) q[1];
rz(3.0473273) q[2];
sx q[2];
rz(-0.69957367) q[2];
sx q[2];
rz(2.6852644) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0316511) q[1];
sx q[1];
rz(-1.9579971) q[1];
sx q[1];
rz(-0.18020682) q[1];
x q[2];
rz(0.40364175) q[3];
sx q[3];
rz(-1.7701245) q[3];
sx q[3];
rz(-2.2011873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74464166) q[2];
sx q[2];
rz(-1.3775185) q[2];
sx q[2];
rz(0.89293876) q[2];
rz(-1.8630155) q[3];
sx q[3];
rz(-1.9847001) q[3];
sx q[3];
rz(-1.6781835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8466012) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(0.60605979) q[0];
rz(3.1312969) q[1];
sx q[1];
rz(-2.1538487) q[1];
sx q[1];
rz(3.0616679) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91658501) q[0];
sx q[0];
rz(-1.3189303) q[0];
sx q[0];
rz(-0.99117898) q[0];
rz(-2.8541508) q[2];
sx q[2];
rz(-0.86311695) q[2];
sx q[2];
rz(-1.8302994) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2458107) q[1];
sx q[1];
rz(-1.5485829) q[1];
sx q[1];
rz(0.40708812) q[1];
x q[2];
rz(1.8013838) q[3];
sx q[3];
rz(-1.7713431) q[3];
sx q[3];
rz(-1.6236562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.14110485) q[2];
sx q[2];
rz(-0.92685574) q[2];
sx q[2];
rz(2.1263988) q[2];
rz(2.5403533) q[3];
sx q[3];
rz(-0.70178086) q[3];
sx q[3];
rz(2.0950441) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64423185) q[0];
sx q[0];
rz(-1.624122) q[0];
sx q[0];
rz(-2.4540785) q[0];
rz(-0.58912206) q[1];
sx q[1];
rz(-2.4644869) q[1];
sx q[1];
rz(-0.45225515) q[1];
rz(2.9420058) q[2];
sx q[2];
rz(-2.6885586) q[2];
sx q[2];
rz(2.7282245) q[2];
rz(-1.4429055) q[3];
sx q[3];
rz(-2.3913132) q[3];
sx q[3];
rz(-2.0146418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
