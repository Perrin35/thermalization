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
rz(0.053057916) q[0];
sx q[0];
rz(-2.7795656) q[0];
sx q[0];
rz(-1.1658573) q[0];
rz(-0.8968269) q[1];
sx q[1];
rz(-1.4520175) q[1];
sx q[1];
rz(-1.7133763) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13442909) q[0];
sx q[0];
rz(-2.0822133) q[0];
sx q[0];
rz(-1.5017139) q[0];
rz(-pi) q[1];
rz(2.8642902) q[2];
sx q[2];
rz(-2.8805974) q[2];
sx q[2];
rz(2.9948611) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.00014405202) q[1];
sx q[1];
rz(-0.23113092) q[1];
sx q[1];
rz(-0.96047251) q[1];
rz(-pi) q[2];
rz(1.0827829) q[3];
sx q[3];
rz(-1.3329771) q[3];
sx q[3];
rz(2.2453515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6894138) q[2];
sx q[2];
rz(-3.1248326) q[2];
sx q[2];
rz(2.9407799) q[2];
rz(2.9928442) q[3];
sx q[3];
rz(-0.0047618682) q[3];
sx q[3];
rz(2.8301921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0018604) q[0];
sx q[0];
rz(-2.5480324) q[0];
sx q[0];
rz(1.0396022) q[0];
rz(3.1270341) q[1];
sx q[1];
rz(-1.2332375) q[1];
sx q[1];
rz(-1.5537517) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7717465) q[0];
sx q[0];
rz(-2.2622007) q[0];
sx q[0];
rz(-0.86172523) q[0];
x q[1];
rz(1.3755685) q[2];
sx q[2];
rz(-0.075610925) q[2];
sx q[2];
rz(1.4812673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5841516) q[1];
sx q[1];
rz(-1.6049892) q[1];
sx q[1];
rz(1.8299915) q[1];
x q[2];
rz(0.513345) q[3];
sx q[3];
rz(-2.7637595) q[3];
sx q[3];
rz(-0.53841089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62850922) q[2];
sx q[2];
rz(-1.5336978) q[2];
sx q[2];
rz(1.7536564) q[2];
rz(1.7657109) q[3];
sx q[3];
rz(-1.0466156) q[3];
sx q[3];
rz(-2.8468813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3917711) q[0];
sx q[0];
rz(-0.23673683) q[0];
sx q[0];
rz(-2.5340875) q[0];
rz(1.5433743) q[1];
sx q[1];
rz(-2.9608455) q[1];
sx q[1];
rz(-2.1764596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9133234) q[0];
sx q[0];
rz(-1.5748657) q[0];
sx q[0];
rz(-1.5913056) q[0];
rz(-pi) q[1];
rz(2.7982462) q[2];
sx q[2];
rz(-2.008524) q[2];
sx q[2];
rz(-0.9870607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7399921) q[1];
sx q[1];
rz(-1.5349746) q[1];
sx q[1];
rz(-1.4261246) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4749281) q[3];
sx q[3];
rz(-1.6536668) q[3];
sx q[3];
rz(1.7516983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33660108) q[2];
sx q[2];
rz(-0.6707297) q[2];
sx q[2];
rz(2.2750308) q[2];
rz(-2.0349515) q[3];
sx q[3];
rz(-1.589078) q[3];
sx q[3];
rz(1.470587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.5690145) q[0];
sx q[0];
rz(-2.4622279) q[0];
sx q[0];
rz(1.5255852) q[0];
rz(-3.1309879) q[1];
sx q[1];
rz(-0.0037825982) q[1];
sx q[1];
rz(-2.3912281) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7667267) q[0];
sx q[0];
rz(-0.9095053) q[0];
sx q[0];
rz(1.4505211) q[0];
rz(2.5360944) q[2];
sx q[2];
rz(-1.7586305) q[2];
sx q[2];
rz(-0.8643291) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.039783104) q[1];
sx q[1];
rz(-1.0652115) q[1];
sx q[1];
rz(1.4391446) q[1];
rz(-pi) q[2];
rz(-1.697465) q[3];
sx q[3];
rz(-1.2614246) q[3];
sx q[3];
rz(-0.29361967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7287207) q[2];
sx q[2];
rz(-2.0527288) q[2];
sx q[2];
rz(1.9000165) q[2];
rz(-0.0056754644) q[3];
sx q[3];
rz(-2.3311876) q[3];
sx q[3];
rz(-0.84398213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64549696) q[0];
sx q[0];
rz(-3.0912919) q[0];
sx q[0];
rz(-2.2182933) q[0];
rz(0.80054379) q[1];
sx q[1];
rz(-3.1381331) q[1];
sx q[1];
rz(-2.961535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6430017) q[0];
sx q[0];
rz(-2.8992436) q[0];
sx q[0];
rz(1.7447628) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6512958) q[2];
sx q[2];
rz(-1.668612) q[2];
sx q[2];
rz(-2.7076633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.99111588) q[1];
sx q[1];
rz(-1.2186945) q[1];
sx q[1];
rz(-1.8262509) q[1];
rz(-pi) q[2];
rz(1.9365129) q[3];
sx q[3];
rz(-1.4673181) q[3];
sx q[3];
rz(2.8998855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8397612) q[2];
sx q[2];
rz(-1.8415035) q[2];
sx q[2];
rz(1.7054455) q[2];
rz(1.2561717) q[3];
sx q[3];
rz(-1.4830282) q[3];
sx q[3];
rz(0.095631599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489478) q[0];
sx q[0];
rz(-2.5697932) q[0];
sx q[0];
rz(-0.44334626) q[0];
rz(1.1706785) q[1];
sx q[1];
rz(-3.1405293) q[1];
sx q[1];
rz(0.43177691) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1373119) q[0];
sx q[0];
rz(-2.3385424) q[0];
sx q[0];
rz(-2.7025168) q[0];
rz(-pi) q[1];
rz(-3.058039) q[2];
sx q[2];
rz(-1.3992893) q[2];
sx q[2];
rz(0.87860859) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1237632) q[1];
sx q[1];
rz(-0.69849724) q[1];
sx q[1];
rz(2.5269954) q[1];
rz(0.59817846) q[3];
sx q[3];
rz(-2.7477187) q[3];
sx q[3];
rz(-2.1109714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40054896) q[2];
sx q[2];
rz(-0.86047518) q[2];
sx q[2];
rz(-1.384037) q[2];
rz(0.69796973) q[3];
sx q[3];
rz(-0.75708404) q[3];
sx q[3];
rz(2.82011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8292002) q[0];
sx q[0];
rz(-1.7990524) q[0];
sx q[0];
rz(2.7685352) q[0];
rz(2.864605) q[1];
sx q[1];
rz(-3.1412558) q[1];
sx q[1];
rz(-2.3854947) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33462229) q[0];
sx q[0];
rz(-0.96617132) q[0];
sx q[0];
rz(1.104825) q[0];
rz(-pi) q[1];
rz(2.7268488) q[2];
sx q[2];
rz(-1.9090134) q[2];
sx q[2];
rz(1.2646706) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0758772) q[1];
sx q[1];
rz(-2.1098032) q[1];
sx q[1];
rz(-0.4735986) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30310615) q[3];
sx q[3];
rz(-1.9744996) q[3];
sx q[3];
rz(1.2991326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21195352) q[2];
sx q[2];
rz(-2.5888207) q[2];
sx q[2];
rz(0.92145222) q[2];
rz(0.067342162) q[3];
sx q[3];
rz(-1.9332956) q[3];
sx q[3];
rz(-1.6578081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15143722) q[0];
sx q[0];
rz(-0.28011265) q[0];
sx q[0];
rz(-0.13023278) q[0];
rz(2.3479334) q[1];
sx q[1];
rz(-0.001611324) q[1];
sx q[1];
rz(-2.8614955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02592382) q[0];
sx q[0];
rz(-1.6243132) q[0];
sx q[0];
rz(-1.491886) q[0];
rz(-pi) q[1];
rz(-2.75861) q[2];
sx q[2];
rz(-2.9071701) q[2];
sx q[2];
rz(-1.2146666) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4520766) q[1];
sx q[1];
rz(-2.048564) q[1];
sx q[1];
rz(0.46124752) q[1];
x q[2];
rz(2.4305268) q[3];
sx q[3];
rz(-1.1991522) q[3];
sx q[3];
rz(-1.7106595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.085999504) q[2];
sx q[2];
rz(-1.4693825) q[2];
sx q[2];
rz(0.80297339) q[2];
rz(-1.5351852) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(2.3101961) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11237535) q[0];
sx q[0];
rz(-3.1387098) q[0];
sx q[0];
rz(3.032384) q[0];
rz(2.7681007) q[1];
sx q[1];
rz(-1.9544574) q[1];
sx q[1];
rz(2.5901897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38574252) q[0];
sx q[0];
rz(-1.047121) q[0];
sx q[0];
rz(0.68501212) q[0];
rz(-pi) q[1];
rz(0.59342758) q[2];
sx q[2];
rz(-2.9205236) q[2];
sx q[2];
rz(-1.9367557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99433339) q[1];
sx q[1];
rz(-0.87169391) q[1];
sx q[1];
rz(-0.84529133) q[1];
x q[2];
rz(-2.738036) q[3];
sx q[3];
rz(-1.7975032) q[3];
sx q[3];
rz(-1.3317396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6138844) q[2];
sx q[2];
rz(-1.8324499) q[2];
sx q[2];
rz(1.8180397) q[2];
rz(1.8771133) q[3];
sx q[3];
rz(-1.2885965) q[3];
sx q[3];
rz(0.0059787353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5962113) q[0];
sx q[0];
rz(-2.5061506) q[0];
sx q[0];
rz(0.7363466) q[0];
rz(2.9296854) q[1];
sx q[1];
rz(-1.0286464) q[1];
sx q[1];
rz(1.5440936) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0422393) q[0];
sx q[0];
rz(-2.8618715) q[0];
sx q[0];
rz(3.0763402) q[0];
x q[1];
rz(1.5681463) q[2];
sx q[2];
rz(-1.6011597) q[2];
sx q[2];
rz(-1.4043491) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5527344) q[1];
sx q[1];
rz(-0.97179669) q[1];
sx q[1];
rz(-2.1373279) q[1];
x q[2];
rz(2.6968217) q[3];
sx q[3];
rz(-1.8666779) q[3];
sx q[3];
rz(2.5887161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.836901) q[2];
sx q[2];
rz(-2.3062314) q[2];
sx q[2];
rz(-1.2738073) q[2];
rz(1.6972208) q[3];
sx q[3];
rz(-3.0675409) q[3];
sx q[3];
rz(-1.4950289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.973751) q[0];
sx q[0];
rz(-1.583562) q[0];
sx q[0];
rz(-1.2927443) q[0];
rz(-1.6043067) q[1];
sx q[1];
rz(-0.91265408) q[1];
sx q[1];
rz(0.18462054) q[1];
rz(0.22278788) q[2];
sx q[2];
rz(-3.1004799) q[2];
sx q[2];
rz(-2.0903175) q[2];
rz(2.1783185) q[3];
sx q[3];
rz(-1.4701075) q[3];
sx q[3];
rz(-0.95095271) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
