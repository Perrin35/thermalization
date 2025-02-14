OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7425793) q[0];
sx q[0];
rz(-1.9785545) q[0];
sx q[0];
rz(1.1385588) q[0];
rz(-0.63602716) q[1];
sx q[1];
rz(-0.50995246) q[1];
sx q[1];
rz(0.62503254) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63955414) q[0];
sx q[0];
rz(-2.9939277) q[0];
sx q[0];
rz(2.8390711) q[0];
rz(-3.1209645) q[2];
sx q[2];
rz(-1.0307023) q[2];
sx q[2];
rz(-2.2079225) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4506065) q[1];
sx q[1];
rz(-2.4391973) q[1];
sx q[1];
rz(-2.2390963) q[1];
rz(-2.9431683) q[3];
sx q[3];
rz(-2.3188667) q[3];
sx q[3];
rz(0.89694512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2193489) q[2];
sx q[2];
rz(-2.2897828) q[2];
sx q[2];
rz(0.40199486) q[2];
rz(-0.829202) q[3];
sx q[3];
rz(-1.821937) q[3];
sx q[3];
rz(1.7792262) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40792337) q[0];
sx q[0];
rz(-2.5897554) q[0];
sx q[0];
rz(-2.0358987) q[0];
rz(2.3528174) q[1];
sx q[1];
rz(-2.5194247) q[1];
sx q[1];
rz(0.50599352) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087990847) q[0];
sx q[0];
rz(-1.8058386) q[0];
sx q[0];
rz(-1.4925692) q[0];
rz(-2.3621735) q[2];
sx q[2];
rz(-0.96923087) q[2];
sx q[2];
rz(1.0600519) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11720235) q[1];
sx q[1];
rz(-1.3474696) q[1];
sx q[1];
rz(-1.9322371) q[1];
rz(-pi) q[2];
rz(-2.6873029) q[3];
sx q[3];
rz(-1.3338106) q[3];
sx q[3];
rz(1.2302356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.904423) q[2];
sx q[2];
rz(-1.0639031) q[2];
sx q[2];
rz(2.2696631) q[2];
rz(2.0875841) q[3];
sx q[3];
rz(-2.0619312) q[3];
sx q[3];
rz(1.7006629) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0525518) q[0];
sx q[0];
rz(-2.4854923) q[0];
sx q[0];
rz(-1.22714) q[0];
rz(-0.50651208) q[1];
sx q[1];
rz(-2.3498693) q[1];
sx q[1];
rz(-0.92481771) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8657579) q[0];
sx q[0];
rz(-1.468549) q[0];
sx q[0];
rz(-1.7195549) q[0];
rz(-1.6396937) q[2];
sx q[2];
rz(-2.0259948) q[2];
sx q[2];
rz(-0.26959878) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9021685) q[1];
sx q[1];
rz(-1.268506) q[1];
sx q[1];
rz(-2.2151674) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66092068) q[3];
sx q[3];
rz(-1.875289) q[3];
sx q[3];
rz(-0.083947649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0623124) q[2];
sx q[2];
rz(-1.541905) q[2];
sx q[2];
rz(2.6463032) q[2];
rz(1.7700206) q[3];
sx q[3];
rz(-0.74426952) q[3];
sx q[3];
rz(1.8194958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1247193) q[0];
sx q[0];
rz(-1.158241) q[0];
sx q[0];
rz(-0.88791263) q[0];
rz(1.3814231) q[1];
sx q[1];
rz(-1.0180749) q[1];
sx q[1];
rz(-0.51753128) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4221188) q[0];
sx q[0];
rz(-1.5848012) q[0];
sx q[0];
rz(1.8318569) q[0];
rz(2.5597455) q[2];
sx q[2];
rz(-2.3334425) q[2];
sx q[2];
rz(0.86607546) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.62152109) q[1];
sx q[1];
rz(-1.5806131) q[1];
sx q[1];
rz(-1.6689368) q[1];
rz(-2.6721084) q[3];
sx q[3];
rz(-1.9133798) q[3];
sx q[3];
rz(-2.1367578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.069328221) q[2];
sx q[2];
rz(-1.1257409) q[2];
sx q[2];
rz(-0.98108712) q[2];
rz(0.4782933) q[3];
sx q[3];
rz(-0.35224733) q[3];
sx q[3];
rz(-0.52811629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1961871) q[0];
sx q[0];
rz(-1.7178752) q[0];
sx q[0];
rz(3.1086573) q[0];
rz(-1.6163274) q[1];
sx q[1];
rz(-0.5739637) q[1];
sx q[1];
rz(2.2059435) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5870635) q[0];
sx q[0];
rz(-1.1008223) q[0];
sx q[0];
rz(1.5857459) q[0];
rz(-1.0887371) q[2];
sx q[2];
rz(-2.3080359) q[2];
sx q[2];
rz(-0.6635467) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5838353) q[1];
sx q[1];
rz(-1.6331065) q[1];
sx q[1];
rz(0.77737096) q[1];
x q[2];
rz(-2.557002) q[3];
sx q[3];
rz(-1.7742312) q[3];
sx q[3];
rz(-0.13336059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4719438) q[2];
sx q[2];
rz(-2.7574597) q[2];
sx q[2];
rz(1.6432537) q[2];
rz(0.60254997) q[3];
sx q[3];
rz(-2.3157401) q[3];
sx q[3];
rz(-1.92441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54685408) q[0];
sx q[0];
rz(-2.8475519) q[0];
sx q[0];
rz(-3.102741) q[0];
rz(-0.36000577) q[1];
sx q[1];
rz(-1.729676) q[1];
sx q[1];
rz(-2.9927599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9827844) q[0];
sx q[0];
rz(-0.17915922) q[0];
sx q[0];
rz(2.3590703) q[0];
rz(-2.5720949) q[2];
sx q[2];
rz(-1.6844771) q[2];
sx q[2];
rz(1.4194429) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56797261) q[1];
sx q[1];
rz(-0.96155969) q[1];
sx q[1];
rz(-0.9602169) q[1];
rz(-pi) q[2];
rz(-1.8735165) q[3];
sx q[3];
rz(-1.0060272) q[3];
sx q[3];
rz(0.64498745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0845324) q[2];
sx q[2];
rz(-0.82814211) q[2];
sx q[2];
rz(2.7866936) q[2];
rz(2.0830294) q[3];
sx q[3];
rz(-0.94616977) q[3];
sx q[3];
rz(-2.6088349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8232089) q[0];
sx q[0];
rz(-1.7192778) q[0];
sx q[0];
rz(0.80867714) q[0];
rz(-3.0478364) q[1];
sx q[1];
rz(-0.84051991) q[1];
sx q[1];
rz(2.7309928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5275972) q[0];
sx q[0];
rz(-1.9638914) q[0];
sx q[0];
rz(-0.40531203) q[0];
rz(-0.65356853) q[2];
sx q[2];
rz(-1.7796345) q[2];
sx q[2];
rz(-2.3318568) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4098874) q[1];
sx q[1];
rz(-0.52919878) q[1];
sx q[1];
rz(-0.7284109) q[1];
rz(0.15488829) q[3];
sx q[3];
rz(-0.67017503) q[3];
sx q[3];
rz(2.6763499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14313301) q[2];
sx q[2];
rz(-2.748558) q[2];
sx q[2];
rz(-0.03579363) q[2];
rz(1.3068457) q[3];
sx q[3];
rz(-1.9948317) q[3];
sx q[3];
rz(-2.3329363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.82658139) q[0];
sx q[0];
rz(-2.1666574) q[0];
sx q[0];
rz(-0.63631979) q[0];
rz(2.4782205) q[1];
sx q[1];
rz(-2.0320818) q[1];
sx q[1];
rz(-2.5136307) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5945744) q[0];
sx q[0];
rz(-1.5994521) q[0];
sx q[0];
rz(0.076731971) q[0];
x q[1];
rz(-2.9344995) q[2];
sx q[2];
rz(-0.40488714) q[2];
sx q[2];
rz(1.5804497) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9414569) q[1];
sx q[1];
rz(-1.0147328) q[1];
sx q[1];
rz(2.0651613) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5359587) q[3];
sx q[3];
rz(-1.3834586) q[3];
sx q[3];
rz(-3.0288937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5593354) q[2];
sx q[2];
rz(-0.11056837) q[2];
sx q[2];
rz(1.4165233) q[2];
rz(-1.0995809) q[3];
sx q[3];
rz(-2.1018335) q[3];
sx q[3];
rz(-0.85272461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13903862) q[0];
sx q[0];
rz(-2.250493) q[0];
sx q[0];
rz(0.6148327) q[0];
rz(-0.50827208) q[1];
sx q[1];
rz(-2.1891687) q[1];
sx q[1];
rz(-0.27613786) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3373703) q[0];
sx q[0];
rz(-1.7003248) q[0];
sx q[0];
rz(-1.1418377) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.040784) q[2];
sx q[2];
rz(-1.4333087) q[2];
sx q[2];
rz(1.9078209) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8135173) q[1];
sx q[1];
rz(-2.349252) q[1];
sx q[1];
rz(-0.66285073) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91137852) q[3];
sx q[3];
rz(-1.2967392) q[3];
sx q[3];
rz(-0.94737939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.189956) q[2];
sx q[2];
rz(-1.0389682) q[2];
sx q[2];
rz(-3.0432126) q[2];
rz(-1.9010057) q[3];
sx q[3];
rz(-1.0616579) q[3];
sx q[3];
rz(-3.0978751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28546178) q[0];
sx q[0];
rz(-2.0933445) q[0];
sx q[0];
rz(-2.7833126) q[0];
rz(2.1367392) q[1];
sx q[1];
rz(-0.88468164) q[1];
sx q[1];
rz(-0.34214941) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2712667) q[0];
sx q[0];
rz(-1.7881601) q[0];
sx q[0];
rz(-2.5408265) q[0];
x q[1];
rz(-1.0295632) q[2];
sx q[2];
rz(-2.4839253) q[2];
sx q[2];
rz(-0.56305158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4917716) q[1];
sx q[1];
rz(-2.559548) q[1];
sx q[1];
rz(-0.19606049) q[1];
x q[2];
rz(-1.8588641) q[3];
sx q[3];
rz(-2.7917842) q[3];
sx q[3];
rz(-1.6988848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3442012) q[2];
sx q[2];
rz(-0.48365334) q[2];
sx q[2];
rz(-2.7698216) q[2];
rz(0.19112912) q[3];
sx q[3];
rz(-1.976795) q[3];
sx q[3];
rz(0.42310664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525748) q[0];
sx q[0];
rz(-0.49674635) q[0];
sx q[0];
rz(2.4865271) q[0];
rz(2.7816506) q[1];
sx q[1];
rz(-1.5690201) q[1];
sx q[1];
rz(1.5108861) q[1];
rz(-1.7689479) q[2];
sx q[2];
rz(-2.8336278) q[2];
sx q[2];
rz(-2.8358493) q[2];
rz(1.5549079) q[3];
sx q[3];
rz(-0.75805981) q[3];
sx q[3];
rz(-0.49657878) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
