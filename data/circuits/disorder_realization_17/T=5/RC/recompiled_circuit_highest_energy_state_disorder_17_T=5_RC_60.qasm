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
rz(-1.0626592) q[0];
sx q[0];
rz(-0.81544977) q[0];
sx q[0];
rz(0.71001473) q[0];
rz(-0.31399909) q[1];
sx q[1];
rz(-0.93500885) q[1];
sx q[1];
rz(11.234565) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3478617) q[0];
sx q[0];
rz(-1.8457185) q[0];
sx q[0];
rz(1.3375651) q[0];
rz(-pi) q[1];
rz(-0.62628905) q[2];
sx q[2];
rz(-0.41538996) q[2];
sx q[2];
rz(0.20520575) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5422022) q[1];
sx q[1];
rz(-0.97394651) q[1];
sx q[1];
rz(-0.29920293) q[1];
x q[2];
rz(2.2043742) q[3];
sx q[3];
rz(-1.9301747) q[3];
sx q[3];
rz(-2.0546953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4189202) q[2];
sx q[2];
rz(-1.924943) q[2];
sx q[2];
rz(2.326272) q[2];
rz(0.80960649) q[3];
sx q[3];
rz(-1.4987192) q[3];
sx q[3];
rz(-1.0838375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0631436) q[0];
sx q[0];
rz(-2.0922631) q[0];
sx q[0];
rz(0.11058841) q[0];
rz(0.94509205) q[1];
sx q[1];
rz(-0.4393622) q[1];
sx q[1];
rz(1.5623215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62842436) q[0];
sx q[0];
rz(-1.1187176) q[0];
sx q[0];
rz(-2.3356209) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3198691) q[2];
sx q[2];
rz(-2.5603449) q[2];
sx q[2];
rz(1.7591214) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.72772273) q[1];
sx q[1];
rz(-1.6752104) q[1];
sx q[1];
rz(2.6563679) q[1];
rz(-2.6591127) q[3];
sx q[3];
rz(-2.2594995) q[3];
sx q[3];
rz(2.3107993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9940146) q[2];
sx q[2];
rz(-0.14573228) q[2];
sx q[2];
rz(0.4501403) q[2];
rz(-2.0077997) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(-2.1639737) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56871539) q[0];
sx q[0];
rz(-2.1935538) q[0];
sx q[0];
rz(-2.8602588) q[0];
rz(1.3866792) q[1];
sx q[1];
rz(-1.7117056) q[1];
sx q[1];
rz(2.5414355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48334405) q[0];
sx q[0];
rz(-1.5243825) q[0];
sx q[0];
rz(-0.29057403) q[0];
x q[1];
rz(0.96397206) q[2];
sx q[2];
rz(-1.6168878) q[2];
sx q[2];
rz(2.5986236) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9439796) q[1];
sx q[1];
rz(-0.67786067) q[1];
sx q[1];
rz(0.53450905) q[1];
x q[2];
rz(1.7417817) q[3];
sx q[3];
rz(-2.1372843) q[3];
sx q[3];
rz(-1.2681707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0448138) q[2];
sx q[2];
rz(-1.1324977) q[2];
sx q[2];
rz(2.0951994) q[2];
rz(-2.7438296) q[3];
sx q[3];
rz(-2.4989276) q[3];
sx q[3];
rz(0.1325632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1716877) q[0];
sx q[0];
rz(-0.02709087) q[0];
sx q[0];
rz(-1.0525674) q[0];
rz(-1.196208) q[1];
sx q[1];
rz(-1.9320678) q[1];
sx q[1];
rz(0.41950163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013166817) q[0];
sx q[0];
rz(-1.3072398) q[0];
sx q[0];
rz(0.58851425) q[0];
rz(-1.2673032) q[2];
sx q[2];
rz(-2.4848273) q[2];
sx q[2];
rz(-3.0820897) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0191325) q[1];
sx q[1];
rz(-1.1943218) q[1];
sx q[1];
rz(-2.2458162) q[1];
x q[2];
rz(-0.49286263) q[3];
sx q[3];
rz(-2.5774084) q[3];
sx q[3];
rz(-1.197937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18998751) q[2];
sx q[2];
rz(-1.1185458) q[2];
sx q[2];
rz(2.0513963) q[2];
rz(-0.53127855) q[3];
sx q[3];
rz(-0.71610206) q[3];
sx q[3];
rz(-1.6147015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7200318) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(1.3979727) q[0];
rz(0.13294237) q[1];
sx q[1];
rz(-1.839919) q[1];
sx q[1];
rz(3.1403819) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5842297) q[0];
sx q[0];
rz(-1.5546636) q[0];
sx q[0];
rz(0.38059142) q[0];
x q[1];
rz(1.8563849) q[2];
sx q[2];
rz(-1.4175709) q[2];
sx q[2];
rz(-1.1040579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9922793) q[1];
sx q[1];
rz(-2.1985801) q[1];
sx q[1];
rz(1.9491484) q[1];
rz(-pi) q[2];
rz(1.7204956) q[3];
sx q[3];
rz(-1.5970917) q[3];
sx q[3];
rz(-2.5291075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2289537) q[2];
sx q[2];
rz(-1.3372083) q[2];
sx q[2];
rz(0.21052989) q[2];
rz(-2.1052965) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(1.2150631) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9482816) q[0];
sx q[0];
rz(-1.223215) q[0];
sx q[0];
rz(2.7358828) q[0];
rz(-1.7871208) q[1];
sx q[1];
rz(-2.3190119) q[1];
sx q[1];
rz(-2.2192661) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2887696) q[0];
sx q[0];
rz(-1.4643702) q[0];
sx q[0];
rz(-1.861051) q[0];
rz(-2.2375003) q[2];
sx q[2];
rz(-1.4897926) q[2];
sx q[2];
rz(1.1912277) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2353224) q[1];
sx q[1];
rz(-2.1791237) q[1];
sx q[1];
rz(3.0023605) q[1];
x q[2];
rz(-2.5821463) q[3];
sx q[3];
rz(-1.9044975) q[3];
sx q[3];
rz(2.4041686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7102082) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(-2.7396438) q[2];
rz(-1.2881783) q[3];
sx q[3];
rz(-2.2737019) q[3];
sx q[3];
rz(-1.8624381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61889082) q[0];
sx q[0];
rz(-2.0382477) q[0];
sx q[0];
rz(-3.0772305) q[0];
rz(-0.83802682) q[1];
sx q[1];
rz(-2.3152654) q[1];
sx q[1];
rz(-1.8109969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143599) q[0];
sx q[0];
rz(-1.4454675) q[0];
sx q[0];
rz(-0.57526875) q[0];
rz(-pi) q[1];
rz(0.60637252) q[2];
sx q[2];
rz(-1.9909711) q[2];
sx q[2];
rz(0.5768896) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2423252) q[1];
sx q[1];
rz(-2.0914255) q[1];
sx q[1];
rz(-1.8249056) q[1];
rz(-pi) q[2];
rz(-0.95860062) q[3];
sx q[3];
rz(-1.7122088) q[3];
sx q[3];
rz(2.8116016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7360721) q[2];
sx q[2];
rz(-1.6330999) q[2];
sx q[2];
rz(-2.4580477) q[2];
rz(2.4077967) q[3];
sx q[3];
rz(-1.7283231) q[3];
sx q[3];
rz(0.84764135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34506327) q[0];
sx q[0];
rz(-2.0579484) q[0];
sx q[0];
rz(1.1731359) q[0];
rz(-0.37551156) q[1];
sx q[1];
rz(-0.48986062) q[1];
sx q[1];
rz(1.0367905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79085892) q[0];
sx q[0];
rz(-3.1111801) q[0];
sx q[0];
rz(-2.4615672) q[0];
x q[1];
rz(-1.0388184) q[2];
sx q[2];
rz(-2.8237282) q[2];
sx q[2];
rz(1.7828072) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6822676) q[1];
sx q[1];
rz(-0.81574355) q[1];
sx q[1];
rz(1.3774648) q[1];
x q[2];
rz(2.100714) q[3];
sx q[3];
rz(-1.138759) q[3];
sx q[3];
rz(-1.0505067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5658687) q[2];
sx q[2];
rz(-2.2278991) q[2];
sx q[2];
rz(2.4110528) q[2];
rz(-2.9324487) q[3];
sx q[3];
rz(-2.6181965) q[3];
sx q[3];
rz(1.8714347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65649477) q[0];
sx q[0];
rz(-0.42335835) q[0];
sx q[0];
rz(0.04743162) q[0];
rz(-2.9196396) q[1];
sx q[1];
rz(-1.0739948) q[1];
sx q[1];
rz(2.2304631) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61924261) q[0];
sx q[0];
rz(-1.5372835) q[0];
sx q[0];
rz(-1.5977809) q[0];
rz(-pi) q[1];
rz(3.1268551) q[2];
sx q[2];
rz(-0.92916766) q[2];
sx q[2];
rz(-0.82889885) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.641727) q[1];
sx q[1];
rz(-1.2033101) q[1];
sx q[1];
rz(-1.5843452) q[1];
x q[2];
rz(-1.0044596) q[3];
sx q[3];
rz(-1.7438423) q[3];
sx q[3];
rz(2.6030948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86965108) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(0.16858777) q[2];
rz(2.9224959) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(-1.8144089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.99437) q[0];
sx q[0];
rz(-0.34039012) q[0];
sx q[0];
rz(-0.45743531) q[0];
rz(-1.9153197) q[1];
sx q[1];
rz(-2.8044082) q[1];
sx q[1];
rz(-1.3630684) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4922614) q[0];
sx q[0];
rz(-2.1773585) q[0];
sx q[0];
rz(-1.6442293) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8149557) q[2];
sx q[2];
rz(-0.95521046) q[2];
sx q[2];
rz(-0.89060874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.05039617) q[1];
sx q[1];
rz(-1.947787) q[1];
sx q[1];
rz(-1.7524377) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2533204) q[3];
sx q[3];
rz(-1.6366881) q[3];
sx q[3];
rz(2.1988188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3739796) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(0.45200959) q[2];
rz(2.7376145) q[3];
sx q[3];
rz(-1.890506) q[3];
sx q[3];
rz(2.0154791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6954738) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(0.52234621) q[1];
sx q[1];
rz(-2.6112687) q[1];
sx q[1];
rz(-1.912259) q[1];
rz(-1.8277373) q[2];
sx q[2];
rz(-1.2972472) q[2];
sx q[2];
rz(-2.2157584) q[2];
rz(1.4215076) q[3];
sx q[3];
rz(-1.310077) q[3];
sx q[3];
rz(2.6826912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
