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
rz(2.8528557) q[0];
sx q[0];
rz(-0.69536916) q[0];
sx q[0];
rz(-2.8737972) q[0];
rz(0.42203045) q[1];
sx q[1];
rz(4.0623436) q[1];
sx q[1];
rz(10.698591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0604297) q[0];
sx q[0];
rz(-1.3295242) q[0];
sx q[0];
rz(0.046242899) q[0];
x q[1];
rz(-2.535475) q[2];
sx q[2];
rz(-0.96783468) q[2];
sx q[2];
rz(-2.9204592) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2886658) q[1];
sx q[1];
rz(-1.5024606) q[1];
sx q[1];
rz(2.4803376) q[1];
rz(-pi) q[2];
rz(1.7142322) q[3];
sx q[3];
rz(-1.6063476) q[3];
sx q[3];
rz(-1.8604904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6196809) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(3.0989975) q[2];
rz(2.8604782) q[3];
sx q[3];
rz(-1.5646076) q[3];
sx q[3];
rz(-0.59578305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5113145) q[0];
sx q[0];
rz(-1.1742641) q[0];
sx q[0];
rz(-1.0761155) q[0];
rz(-2.1108744) q[1];
sx q[1];
rz(-2.0298256) q[1];
sx q[1];
rz(1.3105185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1210107) q[0];
sx q[0];
rz(-0.91269976) q[0];
sx q[0];
rz(-0.026348635) q[0];
rz(2.0115701) q[2];
sx q[2];
rz(-0.78579575) q[2];
sx q[2];
rz(0.94968972) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6678929) q[1];
sx q[1];
rz(-1.4034499) q[1];
sx q[1];
rz(-1.6909864) q[1];
x q[2];
rz(0.63046534) q[3];
sx q[3];
rz(-2.5098636) q[3];
sx q[3];
rz(-0.28030685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.217546) q[2];
sx q[2];
rz(-0.7520389) q[2];
sx q[2];
rz(0.95620608) q[2];
rz(-1.8337967) q[3];
sx q[3];
rz(-0.81495133) q[3];
sx q[3];
rz(-0.1213049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8420551) q[0];
sx q[0];
rz(-0.50899035) q[0];
sx q[0];
rz(-0.036238413) q[0];
rz(-2.3402975) q[1];
sx q[1];
rz(-1.5472629) q[1];
sx q[1];
rz(1.3005728) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42388445) q[0];
sx q[0];
rz(-0.1926271) q[0];
sx q[0];
rz(-0.9813375) q[0];
x q[1];
rz(0.043134281) q[2];
sx q[2];
rz(-0.97402527) q[2];
sx q[2];
rz(-1.5297001) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76884586) q[1];
sx q[1];
rz(-2.0018491) q[1];
sx q[1];
rz(-0.94668862) q[1];
rz(2.7058296) q[3];
sx q[3];
rz(-0.73832694) q[3];
sx q[3];
rz(-2.3070525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7654045) q[2];
sx q[2];
rz(-1.7690965) q[2];
sx q[2];
rz(-0.21793951) q[2];
rz(-0.91207063) q[3];
sx q[3];
rz(-1.4927161) q[3];
sx q[3];
rz(2.2899341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7243778) q[0];
sx q[0];
rz(-1.0715002) q[0];
sx q[0];
rz(0.81556129) q[0];
rz(2.0888445) q[1];
sx q[1];
rz(-0.88200724) q[1];
sx q[1];
rz(1.7900593) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049168436) q[0];
sx q[0];
rz(-1.4935078) q[0];
sx q[0];
rz(-0.46991761) q[0];
rz(-0.10064023) q[2];
sx q[2];
rz(-0.71906861) q[2];
sx q[2];
rz(1.6883862) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38511577) q[1];
sx q[1];
rz(-2.6653892) q[1];
sx q[1];
rz(-1.508273) q[1];
x q[2];
rz(-0.71664401) q[3];
sx q[3];
rz(-1.1782681) q[3];
sx q[3];
rz(1.6300622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1286596) q[2];
sx q[2];
rz(-1.9482875) q[2];
sx q[2];
rz(0.43221691) q[2];
rz(-0.84248078) q[3];
sx q[3];
rz(-1.0336927) q[3];
sx q[3];
rz(-2.1459818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1603482) q[0];
sx q[0];
rz(-0.46537414) q[0];
sx q[0];
rz(2.5904742) q[0];
rz(-2.5235858) q[1];
sx q[1];
rz(-1.8564686) q[1];
sx q[1];
rz(-1.0689703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81449303) q[0];
sx q[0];
rz(-0.90225039) q[0];
sx q[0];
rz(-1.296968) q[0];
rz(-1.1177914) q[2];
sx q[2];
rz(-2.0108622) q[2];
sx q[2];
rz(-2.5534782) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5139937) q[1];
sx q[1];
rz(-0.74666903) q[1];
sx q[1];
rz(-1.7563266) q[1];
x q[2];
rz(-1.7154416) q[3];
sx q[3];
rz(-1.0624043) q[3];
sx q[3];
rz(1.4163464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7908287) q[2];
sx q[2];
rz(-0.5841693) q[2];
sx q[2];
rz(-0.9551777) q[2];
rz(-0.59349924) q[3];
sx q[3];
rz(-2.1953526) q[3];
sx q[3];
rz(1.745863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7668358) q[0];
sx q[0];
rz(-3.0992442) q[0];
sx q[0];
rz(-1.7497077) q[0];
rz(-1.1426686) q[1];
sx q[1];
rz(-1.3418158) q[1];
sx q[1];
rz(1.8291738) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4785193) q[0];
sx q[0];
rz(-1.0940625) q[0];
sx q[0];
rz(-0.8404503) q[0];
rz(-0.80760132) q[2];
sx q[2];
rz(-0.86059216) q[2];
sx q[2];
rz(1.0654895) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23199546) q[1];
sx q[1];
rz(-1.0568406) q[1];
sx q[1];
rz(-1.6327052) q[1];
rz(-pi) q[2];
x q[2];
rz(3.083701) q[3];
sx q[3];
rz(-1.8857919) q[3];
sx q[3];
rz(2.0712899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.352508) q[2];
sx q[2];
rz(-0.75560537) q[2];
sx q[2];
rz(1.4136723) q[2];
rz(0.46755725) q[3];
sx q[3];
rz(-0.65842015) q[3];
sx q[3];
rz(-0.55268923) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308052) q[0];
sx q[0];
rz(-0.063022114) q[0];
sx q[0];
rz(-2.4600273) q[0];
rz(-1.956578) q[1];
sx q[1];
rz(-2.3763035) q[1];
sx q[1];
rz(-2.3550745) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2241572) q[0];
sx q[0];
rz(-2.6154714) q[0];
sx q[0];
rz(0.71147646) q[0];
rz(-pi) q[1];
rz(-0.94633905) q[2];
sx q[2];
rz(-1.1558487) q[2];
sx q[2];
rz(-3.0111661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1934172) q[1];
sx q[1];
rz(-1.6778291) q[1];
sx q[1];
rz(-0.89148895) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6166572) q[3];
sx q[3];
rz(-0.73165441) q[3];
sx q[3];
rz(1.589244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53175348) q[2];
sx q[2];
rz(-2.9015151) q[2];
sx q[2];
rz(1.5698203) q[2];
rz(-1.8543367) q[3];
sx q[3];
rz(-1.5232892) q[3];
sx q[3];
rz(2.1985445) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59931961) q[0];
sx q[0];
rz(-2.4241408) q[0];
sx q[0];
rz(1.9532816) q[0];
rz(3.1207454) q[1];
sx q[1];
rz(-2.3884845) q[1];
sx q[1];
rz(2.5426224) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5064829) q[0];
sx q[0];
rz(-2.0744893) q[0];
sx q[0];
rz(3.0290108) q[0];
rz(2.2856667) q[2];
sx q[2];
rz(-1.2675261) q[2];
sx q[2];
rz(-2.0963469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4074583) q[1];
sx q[1];
rz(-2.2770513) q[1];
sx q[1];
rz(-0.64532537) q[1];
x q[2];
rz(0.4407626) q[3];
sx q[3];
rz(-2.2513933) q[3];
sx q[3];
rz(-1.1632533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.282436) q[2];
sx q[2];
rz(-1.2500637) q[2];
sx q[2];
rz(0.51521987) q[2];
rz(0.53269261) q[3];
sx q[3];
rz(-1.0849413) q[3];
sx q[3];
rz(0.93713078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027997967) q[0];
sx q[0];
rz(-0.6830712) q[0];
sx q[0];
rz(-2.7711476) q[0];
rz(2.0619552) q[1];
sx q[1];
rz(-2.7307983) q[1];
sx q[1];
rz(-3.0830141) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1616042) q[0];
sx q[0];
rz(-1.7287041) q[0];
sx q[0];
rz(-0.61451332) q[0];
x q[1];
rz(0.04444261) q[2];
sx q[2];
rz(-1.343633) q[2];
sx q[2];
rz(-0.73611605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1778922) q[1];
sx q[1];
rz(-0.78989115) q[1];
sx q[1];
rz(-3.0067985) q[1];
x q[2];
rz(-0.94335572) q[3];
sx q[3];
rz(-0.84057144) q[3];
sx q[3];
rz(0.7920023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8687245) q[2];
sx q[2];
rz(-0.803002) q[2];
sx q[2];
rz(0.33099428) q[2];
rz(3.1281779) q[3];
sx q[3];
rz(-0.9534854) q[3];
sx q[3];
rz(-2.0851871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57824221) q[0];
sx q[0];
rz(-0.42911068) q[0];
sx q[0];
rz(0.44813928) q[0];
rz(0.093712417) q[1];
sx q[1];
rz(-2.8401076) q[1];
sx q[1];
rz(1.9261446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5543723) q[0];
sx q[0];
rz(-2.0036812) q[0];
sx q[0];
rz(3.109193) q[0];
rz(-pi) q[1];
rz(2.3289069) q[2];
sx q[2];
rz(-1.2764837) q[2];
sx q[2];
rz(3.091623) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3443904) q[1];
sx q[1];
rz(-2.0995525) q[1];
sx q[1];
rz(0.31086287) q[1];
x q[2];
rz(-1.2566725) q[3];
sx q[3];
rz(-0.41659714) q[3];
sx q[3];
rz(-0.18942094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1756246) q[2];
sx q[2];
rz(-1.5886687) q[2];
sx q[2];
rz(2.442181) q[2];
rz(-1.2753963) q[3];
sx q[3];
rz(-1.1558775) q[3];
sx q[3];
rz(1.2709966) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70915837) q[0];
sx q[0];
rz(-0.86031886) q[0];
sx q[0];
rz(1.1501089) q[0];
rz(1.746183) q[1];
sx q[1];
rz(-1.2184873) q[1];
sx q[1];
rz(1.7582735) q[1];
rz(-1.4906314) q[2];
sx q[2];
rz(-0.90654984) q[2];
sx q[2];
rz(-0.64150099) q[2];
rz(1.4518573) q[3];
sx q[3];
rz(-1.107286) q[3];
sx q[3];
rz(-2.7881691) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
