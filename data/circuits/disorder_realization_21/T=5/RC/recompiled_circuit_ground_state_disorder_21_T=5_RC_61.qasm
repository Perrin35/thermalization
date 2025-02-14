OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.8747099) q[0];
sx q[0];
rz(4.2031718) q[0];
sx q[0];
rz(9.940552) q[0];
rz(0.90142673) q[1];
sx q[1];
rz(-0.20144784) q[1];
sx q[1];
rz(2.1405061) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7667215) q[0];
sx q[0];
rz(-2.3618638) q[0];
sx q[0];
rz(-2.4089902) q[0];
rz(-pi) q[1];
rz(0.41462173) q[2];
sx q[2];
rz(-1.5650563) q[2];
sx q[2];
rz(1.8513377) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.6902509) q[1];
sx q[1];
rz(-0.81762839) q[1];
sx q[1];
rz(0.62599384) q[1];
rz(-2.9941213) q[3];
sx q[3];
rz(-1.6912973) q[3];
sx q[3];
rz(-1.7013753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72656816) q[2];
sx q[2];
rz(-0.51076204) q[2];
sx q[2];
rz(-1.4600352) q[2];
rz(-2.7728752) q[3];
sx q[3];
rz(-1.5548778) q[3];
sx q[3];
rz(-1.4812428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52861315) q[0];
sx q[0];
rz(-3.0280085) q[0];
sx q[0];
rz(-3.0533277) q[0];
rz(-2.2093692) q[1];
sx q[1];
rz(-2.014522) q[1];
sx q[1];
rz(0.50311911) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4934782) q[0];
sx q[0];
rz(-2.0556306) q[0];
sx q[0];
rz(-1.6877112) q[0];
x q[1];
rz(0.18791942) q[2];
sx q[2];
rz(-2.2320679) q[2];
sx q[2];
rz(2.5235944) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4331665) q[1];
sx q[1];
rz(-1.462877) q[1];
sx q[1];
rz(-0.34924653) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1386836) q[3];
sx q[3];
rz(-2.3140597) q[3];
sx q[3];
rz(1.7167909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3162389) q[2];
sx q[2];
rz(-2.1672858) q[2];
sx q[2];
rz(0.45315722) q[2];
rz(3.1406) q[3];
sx q[3];
rz(-1.5627292) q[3];
sx q[3];
rz(-2.4003975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2082625) q[0];
sx q[0];
rz(-0.18284155) q[0];
sx q[0];
rz(-0.46874794) q[0];
rz(-2.5543429) q[1];
sx q[1];
rz(-1.2882371) q[1];
sx q[1];
rz(-0.25150484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5875781) q[0];
sx q[0];
rz(-1.16246) q[0];
sx q[0];
rz(0.46550444) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38893803) q[2];
sx q[2];
rz(-1.374482) q[2];
sx q[2];
rz(1.0599355) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.1826102) q[1];
sx q[1];
rz(-1.6291367) q[1];
sx q[1];
rz(1.6946409) q[1];
rz(-0.44106828) q[3];
sx q[3];
rz(-1.8905283) q[3];
sx q[3];
rz(-2.5351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59715366) q[2];
sx q[2];
rz(-0.73828283) q[2];
sx q[2];
rz(-3.0944589) q[2];
rz(-0.86826396) q[3];
sx q[3];
rz(-1.9009813) q[3];
sx q[3];
rz(2.6422083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78859162) q[0];
sx q[0];
rz(-0.40591875) q[0];
sx q[0];
rz(1.3141919) q[0];
rz(-0.81002533) q[1];
sx q[1];
rz(-1.6255197) q[1];
sx q[1];
rz(-0.31585082) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0994713) q[0];
sx q[0];
rz(-1.0183304) q[0];
sx q[0];
rz(0.61037678) q[0];
rz(-pi) q[1];
rz(-0.48607488) q[2];
sx q[2];
rz(-2.7856084) q[2];
sx q[2];
rz(-0.69185585) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1882945) q[1];
sx q[1];
rz(-0.85894924) q[1];
sx q[1];
rz(-1.1813335) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9258435) q[3];
sx q[3];
rz(-2.5255754) q[3];
sx q[3];
rz(2.8293138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4874068) q[2];
sx q[2];
rz(-0.6260637) q[2];
sx q[2];
rz(-0.89228863) q[2];
rz(-1.6728632) q[3];
sx q[3];
rz(-1.059831) q[3];
sx q[3];
rz(-2.0784126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47657087) q[0];
sx q[0];
rz(-0.63711089) q[0];
sx q[0];
rz(2.2549905) q[0];
rz(2.2992112) q[1];
sx q[1];
rz(-1.8319172) q[1];
sx q[1];
rz(-1.1096035) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32217596) q[0];
sx q[0];
rz(-2.3482394) q[0];
sx q[0];
rz(-0.94950907) q[0];
rz(0.20528741) q[2];
sx q[2];
rz(-0.73371202) q[2];
sx q[2];
rz(3.044341) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.230727) q[1];
sx q[1];
rz(-0.80692569) q[1];
sx q[1];
rz(0.74300503) q[1];
rz(2.4585633) q[3];
sx q[3];
rz(-2.0583526) q[3];
sx q[3];
rz(-1.4157996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2716918) q[2];
sx q[2];
rz(-0.39125189) q[2];
sx q[2];
rz(2.5539577) q[2];
rz(0.21688004) q[3];
sx q[3];
rz(-1.8606595) q[3];
sx q[3];
rz(1.7873526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0103535) q[0];
sx q[0];
rz(-2.6709747) q[0];
sx q[0];
rz(1.7577897) q[0];
rz(2.9673987) q[1];
sx q[1];
rz(-1.1015588) q[1];
sx q[1];
rz(-1.1879638) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0452703) q[0];
sx q[0];
rz(-0.29106566) q[0];
sx q[0];
rz(1.8610659) q[0];
rz(-pi) q[1];
rz(1.401004) q[2];
sx q[2];
rz(-2.0706958) q[2];
sx q[2];
rz(-0.57725805) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9932258) q[1];
sx q[1];
rz(-1.6794551) q[1];
sx q[1];
rz(0.0309561) q[1];
rz(2.8219696) q[3];
sx q[3];
rz(-2.351056) q[3];
sx q[3];
rz(-2.9755693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4478628) q[2];
sx q[2];
rz(-1.2749981) q[2];
sx q[2];
rz(1.0542144) q[2];
rz(-0.55274719) q[3];
sx q[3];
rz(-0.25944969) q[3];
sx q[3];
rz(-2.0235846) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066910557) q[0];
sx q[0];
rz(-0.56147611) q[0];
sx q[0];
rz(-3.1232324) q[0];
rz(-1.8439937) q[1];
sx q[1];
rz(-1.8074139) q[1];
sx q[1];
rz(2.0265354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0783065) q[0];
sx q[0];
rz(-3.0111599) q[0];
sx q[0];
rz(3.0060769) q[0];
rz(-pi) q[1];
rz(-1.1184022) q[2];
sx q[2];
rz(-1.2873189) q[2];
sx q[2];
rz(2.7615135) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0167078) q[1];
sx q[1];
rz(-0.77335268) q[1];
sx q[1];
rz(-0.70752899) q[1];
rz(-2.2998718) q[3];
sx q[3];
rz(-0.63243619) q[3];
sx q[3];
rz(0.034733437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4516419) q[2];
sx q[2];
rz(-0.74863282) q[2];
sx q[2];
rz(2.3484223) q[2];
rz(2.4801109) q[3];
sx q[3];
rz(-1.7659148) q[3];
sx q[3];
rz(-0.64600265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20872214) q[0];
sx q[0];
rz(-2.0168309) q[0];
sx q[0];
rz(2.9710508) q[0];
rz(2.1216682) q[1];
sx q[1];
rz(-2.7248757) q[1];
sx q[1];
rz(0.65933093) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19457283) q[0];
sx q[0];
rz(-1.6259369) q[0];
sx q[0];
rz(1.5558262) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.061672) q[2];
sx q[2];
rz(-1.1069555) q[2];
sx q[2];
rz(2.1042175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0773191) q[1];
sx q[1];
rz(-0.96076316) q[1];
sx q[1];
rz(1.2011106) q[1];
x q[2];
rz(-0.59472707) q[3];
sx q[3];
rz(-2.2910614) q[3];
sx q[3];
rz(1.3095462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3063804) q[2];
sx q[2];
rz(-0.39432085) q[2];
sx q[2];
rz(-0.96768704) q[2];
rz(-2.5505998) q[3];
sx q[3];
rz(-1.3700181) q[3];
sx q[3];
rz(-1.437457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10385926) q[0];
sx q[0];
rz(-0.93872207) q[0];
sx q[0];
rz(-0.18873225) q[0];
rz(1.0602779) q[1];
sx q[1];
rz(-2.4791368) q[1];
sx q[1];
rz(2.2417384) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8322139) q[0];
sx q[0];
rz(-0.77118528) q[0];
sx q[0];
rz(-2.3126649) q[0];
rz(3.138928) q[2];
sx q[2];
rz(-0.84328534) q[2];
sx q[2];
rz(-1.8553986) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2303992) q[1];
sx q[1];
rz(-2.491103) q[1];
sx q[1];
rz(0.002928811) q[1];
rz(-1.7873618) q[3];
sx q[3];
rz(-0.73288554) q[3];
sx q[3];
rz(-0.93001825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4047644) q[2];
sx q[2];
rz(-1.9177723) q[2];
sx q[2];
rz(-1.169211) q[2];
rz(2.2660008) q[3];
sx q[3];
rz(-2.4647522) q[3];
sx q[3];
rz(-3.0524047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5065696) q[0];
sx q[0];
rz(-1.4590141) q[0];
sx q[0];
rz(-0.42612472) q[0];
rz(-1.2395073) q[1];
sx q[1];
rz(-2.111221) q[1];
sx q[1];
rz(-2.820231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6792435) q[0];
sx q[0];
rz(-2.5658742) q[0];
sx q[0];
rz(-1.4690983) q[0];
rz(-pi) q[1];
rz(-1.1280766) q[2];
sx q[2];
rz(-1.6196847) q[2];
sx q[2];
rz(0.14059925) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1275275) q[1];
sx q[1];
rz(-2.7726538) q[1];
sx q[1];
rz(-1.9762726) q[1];
rz(1.0635183) q[3];
sx q[3];
rz(-0.92417704) q[3];
sx q[3];
rz(2.3648928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1539803) q[2];
sx q[2];
rz(-2.2025547) q[2];
sx q[2];
rz(3.0737446) q[2];
rz(2.1765354) q[3];
sx q[3];
rz(-2.8128251) q[3];
sx q[3];
rz(-1.4062101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52473749) q[0];
sx q[0];
rz(-2.6829834) q[0];
sx q[0];
rz(0.99880698) q[0];
rz(2.2899992) q[1];
sx q[1];
rz(-2.0709745) q[1];
sx q[1];
rz(-2.4877683) q[1];
rz(-2.5851696) q[2];
sx q[2];
rz(-1.6243373) q[2];
sx q[2];
rz(-2.7173964) q[2];
rz(1.5461934) q[3];
sx q[3];
rz(-0.8529803) q[3];
sx q[3];
rz(2.148694) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
