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
rz(2.4438357) q[0];
sx q[0];
rz(-1.9480167) q[0];
sx q[0];
rz(0.20456631) q[0];
rz(-0.76072955) q[1];
sx q[1];
rz(-1.3890356) q[1];
sx q[1];
rz(-1.7420446) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4142864) q[0];
sx q[0];
rz(-2.901863) q[0];
sx q[0];
rz(-1.1778465) q[0];
x q[1];
rz(1.1209773) q[2];
sx q[2];
rz(-1.9391141) q[2];
sx q[2];
rz(0.30003795) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0006932) q[1];
sx q[1];
rz(-1.4834036) q[1];
sx q[1];
rz(0.42229514) q[1];
x q[2];
rz(-1.8379962) q[3];
sx q[3];
rz(-2.2624216) q[3];
sx q[3];
rz(-1.2822145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.43510398) q[2];
sx q[2];
rz(-0.032328345) q[2];
sx q[2];
rz(-0.48933634) q[2];
rz(-1.0232183) q[3];
sx q[3];
rz(-0.018298572) q[3];
sx q[3];
rz(-2.0160915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045687549) q[0];
sx q[0];
rz(-2.4850595) q[0];
sx q[0];
rz(-0.80279654) q[0];
rz(-3.070201) q[1];
sx q[1];
rz(-0.26669058) q[1];
sx q[1];
rz(-3.0838222) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5934385) q[0];
sx q[0];
rz(-1.3134985) q[0];
sx q[0];
rz(2.1307039) q[0];
rz(-0.34681706) q[2];
sx q[2];
rz(-2.8379734) q[2];
sx q[2];
rz(0.70250073) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.36684588) q[1];
sx q[1];
rz(-0.54301942) q[1];
sx q[1];
rz(-2.0319035) q[1];
rz(-pi) q[2];
rz(-0.85912786) q[3];
sx q[3];
rz(-1.6831846) q[3];
sx q[3];
rz(2.4939362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1959261) q[2];
sx q[2];
rz(-2.0464996) q[2];
sx q[2];
rz(1.2954953) q[2];
rz(0.96674353) q[3];
sx q[3];
rz(-0.77015489) q[3];
sx q[3];
rz(-2.3841592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029723786) q[0];
sx q[0];
rz(-1.3628549) q[0];
sx q[0];
rz(1.7080074) q[0];
rz(-3.0729821) q[1];
sx q[1];
rz(-1.5674633) q[1];
sx q[1];
rz(-2.5624018) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96087181) q[0];
sx q[0];
rz(-1.6287043) q[0];
sx q[0];
rz(-1.4790034) q[0];
rz(-pi) q[1];
rz(-3.1379329) q[2];
sx q[2];
rz(-0.47572593) q[2];
sx q[2];
rz(1.2086679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.9186185) q[1];
sx q[1];
rz(-0.22906216) q[1];
sx q[1];
rz(-0.15840662) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3287401) q[3];
sx q[3];
rz(-0.62004706) q[3];
sx q[3];
rz(-0.64921415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.27089831) q[2];
sx q[2];
rz(-1.8879994) q[2];
sx q[2];
rz(2.9581621) q[2];
rz(-0.80426788) q[3];
sx q[3];
rz(-0.99884123) q[3];
sx q[3];
rz(0.53406322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428228) q[0];
sx q[0];
rz(-0.11480055) q[0];
sx q[0];
rz(2.542069) q[0];
rz(2.7291258) q[1];
sx q[1];
rz(-3.1209374) q[1];
sx q[1];
rz(-1.0106769) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3314047) q[0];
sx q[0];
rz(-1.6972901) q[0];
sx q[0];
rz(-2.7888915) q[0];
rz(0.050881906) q[2];
sx q[2];
rz(-2.2023099) q[2];
sx q[2];
rz(-2.8040407) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4356723) q[1];
sx q[1];
rz(-1.3002035) q[1];
sx q[1];
rz(-0.95939221) q[1];
rz(-pi) q[2];
rz(-1.4323727) q[3];
sx q[3];
rz(-2.3218621) q[3];
sx q[3];
rz(2.8734796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7708873) q[2];
sx q[2];
rz(-2.7949896) q[2];
sx q[2];
rz(2.8240805) q[2];
rz(-0.62234771) q[3];
sx q[3];
rz(-0.93572664) q[3];
sx q[3];
rz(-2.5573825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7790826) q[0];
sx q[0];
rz(-2.2267987) q[0];
sx q[0];
rz(-2.1146178) q[0];
rz(2.5912071) q[1];
sx q[1];
rz(-3.0774979) q[1];
sx q[1];
rz(1.9245573) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4452222) q[0];
sx q[0];
rz(-1.0125986) q[0];
sx q[0];
rz(-0.77409111) q[0];
rz(1.4929068) q[2];
sx q[2];
rz(-2.1078034) q[2];
sx q[2];
rz(-0.80918771) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0266307) q[1];
sx q[1];
rz(-1.3991881) q[1];
sx q[1];
rz(2.2725355) q[1];
rz(-1.9868136) q[3];
sx q[3];
rz(-0.75473329) q[3];
sx q[3];
rz(-2.5861135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7050742) q[2];
sx q[2];
rz(-1.7784092) q[2];
sx q[2];
rz(-0.77318937) q[2];
rz(0.1117205) q[3];
sx q[3];
rz(-1.819928) q[3];
sx q[3];
rz(-0.93938655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47950995) q[0];
sx q[0];
rz(-0.22308068) q[0];
sx q[0];
rz(0.42698419) q[0];
rz(2.2110979) q[1];
sx q[1];
rz(-0.016914802) q[1];
sx q[1];
rz(-2.6771136) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31373608) q[0];
sx q[0];
rz(-1.858874) q[0];
sx q[0];
rz(-1.7749857) q[0];
rz(-0.33916766) q[2];
sx q[2];
rz(-2.1438823) q[2];
sx q[2];
rz(-3.1011875) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1982616) q[1];
sx q[1];
rz(-1.2020018) q[1];
sx q[1];
rz(-1.6095227) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6397912) q[3];
sx q[3];
rz(-0.42519266) q[3];
sx q[3];
rz(-2.0824279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6102607) q[2];
sx q[2];
rz(-1.8208296) q[2];
sx q[2];
rz(-2.8533234) q[2];
rz(-2.098295) q[3];
sx q[3];
rz(-0.61560029) q[3];
sx q[3];
rz(2.402795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4814602) q[0];
sx q[0];
rz(-1.2876502) q[0];
sx q[0];
rz(0.75755358) q[0];
rz(0.07269147) q[1];
sx q[1];
rz(-3.1158267) q[1];
sx q[1];
rz(-3.0933948) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5023247) q[0];
sx q[0];
rz(-1.6002965) q[0];
sx q[0];
rz(0.98639368) q[0];
rz(-pi) q[1];
rz(0.9725857) q[2];
sx q[2];
rz(-1.3763104) q[2];
sx q[2];
rz(1.1301646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1485702) q[1];
sx q[1];
rz(-2.3766368) q[1];
sx q[1];
rz(0.92924849) q[1];
x q[2];
rz(0.052560135) q[3];
sx q[3];
rz(-1.1014928) q[3];
sx q[3];
rz(-1.8183501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66182071) q[2];
sx q[2];
rz(-1.6419819) q[2];
sx q[2];
rz(3.0721967) q[2];
rz(-1.5875459) q[3];
sx q[3];
rz(-0.78726751) q[3];
sx q[3];
rz(2.9230996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3875535) q[0];
sx q[0];
rz(-2.0856922) q[0];
sx q[0];
rz(-1.4132502) q[0];
rz(-2.3130401) q[1];
sx q[1];
rz(-0.041752432) q[1];
sx q[1];
rz(0.54263306) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56930243) q[0];
sx q[0];
rz(-2.9539032) q[0];
sx q[0];
rz(2.45298) q[0];
rz(-pi) q[1];
rz(-1.2130402) q[2];
sx q[2];
rz(-1.2160436) q[2];
sx q[2];
rz(-0.32217978) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5396351) q[1];
sx q[1];
rz(-1.4550721) q[1];
sx q[1];
rz(-1.5166111) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7067634) q[3];
sx q[3];
rz(-1.2246338) q[3];
sx q[3];
rz(1.0258254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1594306) q[2];
sx q[2];
rz(-2.7320778) q[2];
sx q[2];
rz(2.8816667) q[2];
rz(-1.0204756) q[3];
sx q[3];
rz(-2.8838938) q[3];
sx q[3];
rz(-0.79403383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.4644153) q[0];
sx q[0];
rz(-2.9554415) q[0];
sx q[0];
rz(-1.4790685) q[0];
rz(-1.6089449) q[1];
sx q[1];
rz(-1.0391935) q[1];
sx q[1];
rz(-2.398568) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047516454) q[0];
sx q[0];
rz(-1.0595982) q[0];
sx q[0];
rz(-0.1316977) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7588927) q[2];
sx q[2];
rz(-2.0819132) q[2];
sx q[2];
rz(-1.2148884) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8919395) q[1];
sx q[1];
rz(-1.6257451) q[1];
sx q[1];
rz(1.6366538) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0366692) q[3];
sx q[3];
rz(-1.6015953) q[3];
sx q[3];
rz(-1.2237751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.6741901) q[2];
sx q[2];
rz(-0.82618606) q[2];
sx q[2];
rz(0.80545938) q[2];
rz(-1.6921267) q[3];
sx q[3];
rz(-1.2286681) q[3];
sx q[3];
rz(0.8031556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2529124) q[0];
sx q[0];
rz(-2.5997933) q[0];
sx q[0];
rz(-2.3895277) q[0];
rz(-1.9883142) q[1];
sx q[1];
rz(-2.2572932) q[1];
sx q[1];
rz(2.8582252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3842938) q[0];
sx q[0];
rz(-1.481825) q[0];
sx q[0];
rz(-0.53953895) q[0];
rz(-1.7287909) q[2];
sx q[2];
rz(-0.758095) q[2];
sx q[2];
rz(-1.5811063) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99774573) q[1];
sx q[1];
rz(-2.4334868) q[1];
sx q[1];
rz(2.6216402) q[1];
rz(-pi) q[2];
rz(-2.4872736) q[3];
sx q[3];
rz(-1.2551184) q[3];
sx q[3];
rz(1.1512427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2470384) q[2];
sx q[2];
rz(-0.082823195) q[2];
sx q[2];
rz(-1.7029597) q[2];
rz(-2.8476207) q[3];
sx q[3];
rz(-3.1271264) q[3];
sx q[3];
rz(1.0283874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033584874) q[0];
sx q[0];
rz(-1.7243732) q[0];
sx q[0];
rz(1.617817) q[0];
rz(-2.602018) q[1];
sx q[1];
rz(-0.78782606) q[1];
sx q[1];
rz(0.16001564) q[1];
rz(2.9931184) q[2];
sx q[2];
rz(-1.21429) q[2];
sx q[2];
rz(0.16333632) q[2];
rz(-0.60565518) q[3];
sx q[3];
rz(-2.6867821) q[3];
sx q[3];
rz(2.3985779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
