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
rz(0.19652551) q[0];
sx q[0];
rz(2.3152469) q[0];
sx q[0];
rz(8.5066975) q[0];
rz(3.0300568) q[1];
sx q[1];
rz(-1.7938951) q[1];
sx q[1];
rz(-1.5674051) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67016168) q[0];
sx q[0];
rz(-1.5613371) q[0];
sx q[0];
rz(-0.0034250101) q[0];
x q[1];
rz(0.55337535) q[2];
sx q[2];
rz(-0.63734431) q[2];
sx q[2];
rz(3.1269249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1861856) q[1];
sx q[1];
rz(-1.2689841) q[1];
sx q[1];
rz(0.11037453) q[1];
rz(0.49310103) q[3];
sx q[3];
rz(-1.6776909) q[3];
sx q[3];
rz(1.3809539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.427318) q[2];
sx q[2];
rz(-1.486342) q[2];
sx q[2];
rz(1.8753258) q[2];
rz(-2.0990939) q[3];
sx q[3];
rz(-1.6008585) q[3];
sx q[3];
rz(-1.3955759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.105044) q[0];
sx q[0];
rz(-2.5140913) q[0];
sx q[0];
rz(-3.0446766) q[0];
rz(-0.80822432) q[1];
sx q[1];
rz(-0.40882912) q[1];
sx q[1];
rz(-1.854863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42303434) q[0];
sx q[0];
rz(-2.0204394) q[0];
sx q[0];
rz(3.1021523) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1551844) q[2];
sx q[2];
rz(-0.93743582) q[2];
sx q[2];
rz(1.4937166) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60139537) q[1];
sx q[1];
rz(-0.88078618) q[1];
sx q[1];
rz(-0.50443919) q[1];
rz(-pi) q[2];
rz(2.3318491) q[3];
sx q[3];
rz(-2.6698723) q[3];
sx q[3];
rz(-2.7641486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5010995) q[2];
sx q[2];
rz(-0.88397908) q[2];
sx q[2];
rz(-0.38197771) q[2];
rz(2.6169418) q[3];
sx q[3];
rz(-1.7347521) q[3];
sx q[3];
rz(2.9978571) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76661888) q[0];
sx q[0];
rz(-1.5559649) q[0];
sx q[0];
rz(2.4554456) q[0];
rz(2.0740017) q[1];
sx q[1];
rz(-2.144404) q[1];
sx q[1];
rz(3.0608665) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9322745) q[0];
sx q[0];
rz(-0.22695146) q[0];
sx q[0];
rz(1.7448241) q[0];
rz(-pi) q[1];
rz(2.2171634) q[2];
sx q[2];
rz(-1.0253128) q[2];
sx q[2];
rz(1.1875064) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7551832) q[1];
sx q[1];
rz(-1.5940136) q[1];
sx q[1];
rz(-0.23020102) q[1];
rz(2.3537797) q[3];
sx q[3];
rz(-1.2255049) q[3];
sx q[3];
rz(2.6097176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2637691) q[2];
sx q[2];
rz(-2.4407385) q[2];
sx q[2];
rz(-0.43886718) q[2];
rz(1.8435439) q[3];
sx q[3];
rz(-2.7545007) q[3];
sx q[3];
rz(-2.7264061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2271659) q[0];
sx q[0];
rz(-0.54410797) q[0];
sx q[0];
rz(0.67657226) q[0];
rz(-0.1380955) q[1];
sx q[1];
rz(-1.0905677) q[1];
sx q[1];
rz(0.20060435) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5745068) q[0];
sx q[0];
rz(-2.4100523) q[0];
sx q[0];
rz(1.0163496) q[0];
rz(-0.62941636) q[2];
sx q[2];
rz(-1.417629) q[2];
sx q[2];
rz(-1.6352194) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.296941) q[1];
sx q[1];
rz(-1.703023) q[1];
sx q[1];
rz(2.5297715) q[1];
rz(-pi) q[2];
rz(-1.7142606) q[3];
sx q[3];
rz(-2.1273566) q[3];
sx q[3];
rz(-1.0866764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6382008) q[2];
sx q[2];
rz(-1.8663422) q[2];
sx q[2];
rz(-1.4208043) q[2];
rz(0.11451379) q[3];
sx q[3];
rz(-1.6679461) q[3];
sx q[3];
rz(-2.1421471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223709) q[0];
sx q[0];
rz(-2.350816) q[0];
sx q[0];
rz(2.1694699) q[0];
rz(-0.32360336) q[1];
sx q[1];
rz(-1.4515896) q[1];
sx q[1];
rz(-1.9591029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2401449) q[0];
sx q[0];
rz(-1.5306382) q[0];
sx q[0];
rz(1.7661278) q[0];
rz(-pi) q[1];
rz(-2.5713872) q[2];
sx q[2];
rz(-1.5657164) q[2];
sx q[2];
rz(-3.0722116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4515498) q[1];
sx q[1];
rz(-1.4937972) q[1];
sx q[1];
rz(1.4119215) q[1];
rz(-2.5706815) q[3];
sx q[3];
rz(-1.1215253) q[3];
sx q[3];
rz(-0.8415287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16367308) q[2];
sx q[2];
rz(-1.7719496) q[2];
sx q[2];
rz(-2.6046806) q[2];
rz(-0.72530693) q[3];
sx q[3];
rz(-0.0497497) q[3];
sx q[3];
rz(-0.79881001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2742915) q[0];
sx q[0];
rz(-2.0766356) q[0];
sx q[0];
rz(-1.4428447) q[0];
rz(0.2105712) q[1];
sx q[1];
rz(-1.4434394) q[1];
sx q[1];
rz(-0.64705667) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4763586) q[0];
sx q[0];
rz(-1.3148576) q[0];
sx q[0];
rz(2.7546935) q[0];
rz(-2.5850614) q[2];
sx q[2];
rz(-1.663563) q[2];
sx q[2];
rz(-1.1586729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.9155026) q[1];
sx q[1];
rz(-1.9880297) q[1];
sx q[1];
rz(-1.1933865) q[1];
rz(-2.4652867) q[3];
sx q[3];
rz(-0.42255536) q[3];
sx q[3];
rz(0.12862118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1232542) q[2];
sx q[2];
rz(-1.0945357) q[2];
sx q[2];
rz(-1.9617762) q[2];
rz(1.2849464) q[3];
sx q[3];
rz(-0.40950567) q[3];
sx q[3];
rz(-3.0783317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0354075) q[0];
sx q[0];
rz(-1.4466865) q[0];
sx q[0];
rz(0.89757288) q[0];
rz(1.8073742) q[1];
sx q[1];
rz(-1.370627) q[1];
sx q[1];
rz(-2.1468377) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.008381) q[0];
sx q[0];
rz(-1.8803839) q[0];
sx q[0];
rz(-0.82010834) q[0];
rz(0.95218539) q[2];
sx q[2];
rz(-2.4403009) q[2];
sx q[2];
rz(1.889297) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23235944) q[1];
sx q[1];
rz(-0.5482175) q[1];
sx q[1];
rz(-2.1562804) q[1];
x q[2];
rz(3.1174421) q[3];
sx q[3];
rz(-1.2901879) q[3];
sx q[3];
rz(2.5606972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0576386) q[2];
sx q[2];
rz(-1.3292162) q[2];
sx q[2];
rz(2.8482021) q[2];
rz(3.0883664) q[3];
sx q[3];
rz(-2.2727727) q[3];
sx q[3];
rz(-1.7085541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5695802) q[0];
sx q[0];
rz(-1.7520289) q[0];
sx q[0];
rz(-0.30297512) q[0];
rz(1.6527893) q[1];
sx q[1];
rz(-1.3158512) q[1];
sx q[1];
rz(-2.4724919) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8527864) q[0];
sx q[0];
rz(-1.1256721) q[0];
sx q[0];
rz(1.9381802) q[0];
x q[1];
rz(2.8995048) q[2];
sx q[2];
rz(-1.228516) q[2];
sx q[2];
rz(0.9108327) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4035621) q[1];
sx q[1];
rz(-0.68519095) q[1];
sx q[1];
rz(-2.2735647) q[1];
rz(-pi) q[2];
rz(-0.49578285) q[3];
sx q[3];
rz(-2.6675329) q[3];
sx q[3];
rz(2.1366675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1423433) q[2];
sx q[2];
rz(-1.3331022) q[2];
sx q[2];
rz(2.2744501) q[2];
rz(-1.123318) q[3];
sx q[3];
rz(-2.2893548) q[3];
sx q[3];
rz(-1.1423906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.860723) q[0];
sx q[0];
rz(-2.6779802) q[0];
sx q[0];
rz(-3.0238357) q[0];
rz(0.87751687) q[1];
sx q[1];
rz(-2.120647) q[1];
sx q[1];
rz(2.1868475) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68761008) q[0];
sx q[0];
rz(-1.1907401) q[0];
sx q[0];
rz(2.4847772) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0659201) q[2];
sx q[2];
rz(-1.1469764) q[2];
sx q[2];
rz(-1.1413121) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92312557) q[1];
sx q[1];
rz(-0.72281853) q[1];
sx q[1];
rz(-1.2413526) q[1];
rz(-pi) q[2];
rz(-0.84857984) q[3];
sx q[3];
rz(-0.89451087) q[3];
sx q[3];
rz(0.05300314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8942326) q[2];
sx q[2];
rz(-2.8017513) q[2];
sx q[2];
rz(1.6575238) q[2];
rz(-1.6190716) q[3];
sx q[3];
rz(-1.7584453) q[3];
sx q[3];
rz(1.1744261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5840983) q[0];
sx q[0];
rz(-1.7099986) q[0];
sx q[0];
rz(0.75310055) q[0];
rz(1.151471) q[1];
sx q[1];
rz(-0.52422062) q[1];
sx q[1];
rz(-1.6023191) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26254567) q[0];
sx q[0];
rz(-0.48948727) q[0];
sx q[0];
rz(-1.9681843) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1588155) q[2];
sx q[2];
rz(-0.78938198) q[2];
sx q[2];
rz(-0.090580926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6021612) q[1];
sx q[1];
rz(-1.185158) q[1];
sx q[1];
rz(-1.3630609) q[1];
rz(-0.65493802) q[3];
sx q[3];
rz(-1.7823151) q[3];
sx q[3];
rz(-1.7251273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7842385) q[2];
sx q[2];
rz(-0.18111649) q[2];
sx q[2];
rz(-0.46585807) q[2];
rz(-2.0458131) q[3];
sx q[3];
rz(-0.992479) q[3];
sx q[3];
rz(-0.54212681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096238484) q[0];
sx q[0];
rz(-1.5662554) q[0];
sx q[0];
rz(1.5731496) q[0];
rz(1.0678328) q[1];
sx q[1];
rz(-0.44034958) q[1];
sx q[1];
rz(-0.82028295) q[1];
rz(-2.8779262) q[2];
sx q[2];
rz(-1.7963526) q[2];
sx q[2];
rz(0.11033478) q[2];
rz(2.0910083) q[3];
sx q[3];
rz(-1.0903842) q[3];
sx q[3];
rz(0.36538418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
