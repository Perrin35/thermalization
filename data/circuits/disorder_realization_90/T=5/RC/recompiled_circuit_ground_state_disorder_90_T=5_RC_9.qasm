OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8425771) q[0];
sx q[0];
rz(-0.17248532) q[0];
sx q[0];
rz(0.01699288) q[0];
rz(0.51141557) q[1];
sx q[1];
rz(-0.49632448) q[1];
sx q[1];
rz(-2.6561148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4776909) q[0];
sx q[0];
rz(-1.6833651) q[0];
sx q[0];
rz(1.6825874) q[0];
rz(-pi) q[1];
rz(-2.4713712) q[2];
sx q[2];
rz(-2.1677758) q[2];
sx q[2];
rz(2.1917997) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.770931) q[1];
sx q[1];
rz(-1.071613) q[1];
sx q[1];
rz(0.73727495) q[1];
rz(-pi) q[2];
rz(1.2608246) q[3];
sx q[3];
rz(-1.4616218) q[3];
sx q[3];
rz(-2.7742164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7444676) q[2];
sx q[2];
rz(-3.0391389) q[2];
sx q[2];
rz(2.3490119) q[2];
rz(-2.1553195) q[3];
sx q[3];
rz(-1.6542185) q[3];
sx q[3];
rz(0.13315323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88696402) q[0];
sx q[0];
rz(-1.5381085) q[0];
sx q[0];
rz(-3.0358553) q[0];
rz(-0.75791439) q[1];
sx q[1];
rz(-2.4287903) q[1];
sx q[1];
rz(-2.6656718) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3041046) q[0];
sx q[0];
rz(-1.6679576) q[0];
sx q[0];
rz(2.9062953) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49293002) q[2];
sx q[2];
rz(-1.4441737) q[2];
sx q[2];
rz(-2.7585518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8393499) q[1];
sx q[1];
rz(-2.1363597) q[1];
sx q[1];
rz(-2.9981705) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1844238) q[3];
sx q[3];
rz(-2.4605721) q[3];
sx q[3];
rz(-1.0134361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.14625749) q[2];
sx q[2];
rz(-0.47440752) q[2];
sx q[2];
rz(-1.8246626) q[2];
rz(2.3138192) q[3];
sx q[3];
rz(-1.0931949) q[3];
sx q[3];
rz(-2.304346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4227609) q[0];
sx q[0];
rz(-1.7976924) q[0];
sx q[0];
rz(-1.9047009) q[0];
rz(1.3924567) q[1];
sx q[1];
rz(-1.720865) q[1];
sx q[1];
rz(-2.4512591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25094098) q[0];
sx q[0];
rz(-1.1162045) q[0];
sx q[0];
rz(0.62418749) q[0];
rz(1.1074781) q[2];
sx q[2];
rz(-1.5743739) q[2];
sx q[2];
rz(1.6301375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.954309) q[1];
sx q[1];
rz(-2.885703) q[1];
sx q[1];
rz(2.4268716) q[1];
rz(-pi) q[2];
rz(-0.66392939) q[3];
sx q[3];
rz(-2.2888043) q[3];
sx q[3];
rz(1.688886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51848015) q[2];
sx q[2];
rz(-1.5223576) q[2];
sx q[2];
rz(0.073866455) q[2];
rz(-1.1582003) q[3];
sx q[3];
rz(-1.2169714) q[3];
sx q[3];
rz(-1.7346252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.912792) q[0];
sx q[0];
rz(-0.055963628) q[0];
sx q[0];
rz(-0.19293109) q[0];
rz(1.2436766) q[1];
sx q[1];
rz(-1.7453777) q[1];
sx q[1];
rz(2.9085433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6911963) q[0];
sx q[0];
rz(-1.6653227) q[0];
sx q[0];
rz(-0.14733845) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2795882) q[2];
sx q[2];
rz(-1.5995912) q[2];
sx q[2];
rz(-0.51143247) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8687081) q[1];
sx q[1];
rz(-1.2603425) q[1];
sx q[1];
rz(-1.5116879) q[1];
x q[2];
rz(-0.95110622) q[3];
sx q[3];
rz(-2.4220805) q[3];
sx q[3];
rz(0.48966416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21726808) q[2];
sx q[2];
rz(-1.3793719) q[2];
sx q[2];
rz(0.064083727) q[2];
rz(-0.25121769) q[3];
sx q[3];
rz(-2.678674) q[3];
sx q[3];
rz(2.8502051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578882) q[0];
sx q[0];
rz(-1.2936445) q[0];
sx q[0];
rz(-1.2842913) q[0];
rz(-0.0074370782) q[1];
sx q[1];
rz(-1.9556655) q[1];
sx q[1];
rz(-2.1844905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.590132) q[0];
sx q[0];
rz(-0.55319417) q[0];
sx q[0];
rz(-0.93450089) q[0];
rz(-pi) q[1];
rz(0.77218036) q[2];
sx q[2];
rz(-2.3640458) q[2];
sx q[2];
rz(-0.61138703) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6308208) q[1];
sx q[1];
rz(-1.242852) q[1];
sx q[1];
rz(-1.4664709) q[1];
x q[2];
rz(2.4864205) q[3];
sx q[3];
rz(-1.7196894) q[3];
sx q[3];
rz(-2.6096491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0127504) q[2];
sx q[2];
rz(-2.4350171) q[2];
sx q[2];
rz(2.1112704) q[2];
rz(2.4230867) q[3];
sx q[3];
rz(-1.0822783) q[3];
sx q[3];
rz(2.4350186) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4265823) q[0];
sx q[0];
rz(-0.74786818) q[0];
sx q[0];
rz(2.7744875) q[0];
rz(2.7857065) q[1];
sx q[1];
rz(-2.0242736) q[1];
sx q[1];
rz(1.5497367) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132431) q[0];
sx q[0];
rz(-1.5322881) q[0];
sx q[0];
rz(-2.4000077) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0075304) q[2];
sx q[2];
rz(-1.4685837) q[2];
sx q[2];
rz(1.619316) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8198135) q[1];
sx q[1];
rz(-1.1604619) q[1];
sx q[1];
rz(1.3458662) q[1];
rz(2.6713085) q[3];
sx q[3];
rz(-1.6374644) q[3];
sx q[3];
rz(2.0955412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1648569) q[2];
sx q[2];
rz(-1.2931436) q[2];
sx q[2];
rz(3.1381651) q[2];
rz(2.2554743) q[3];
sx q[3];
rz(-2.0485853) q[3];
sx q[3];
rz(1.4440943) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43585983) q[0];
sx q[0];
rz(-1.4737361) q[0];
sx q[0];
rz(-0.71520299) q[0];
rz(-0.12229478) q[1];
sx q[1];
rz(-1.0194174) q[1];
sx q[1];
rz(0.14762793) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24623309) q[0];
sx q[0];
rz(-2.6608855) q[0];
sx q[0];
rz(2.3395545) q[0];
x q[1];
rz(2.4164957) q[2];
sx q[2];
rz(-2.5864961) q[2];
sx q[2];
rz(-1.8738418) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9761889) q[1];
sx q[1];
rz(-2.2641085) q[1];
sx q[1];
rz(2.90392) q[1];
rz(-pi) q[2];
x q[2];
rz(1.857599) q[3];
sx q[3];
rz(-1.6383324) q[3];
sx q[3];
rz(0.90112858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3132396) q[2];
sx q[2];
rz(-0.30892631) q[2];
sx q[2];
rz(0.1114791) q[2];
rz(-2.3765423) q[3];
sx q[3];
rz(-1.7181516) q[3];
sx q[3];
rz(-0.033871977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2400804) q[0];
sx q[0];
rz(-2.1724048) q[0];
sx q[0];
rz(-1.0028268) q[0];
rz(-0.78549939) q[1];
sx q[1];
rz(-1.5248884) q[1];
sx q[1];
rz(1.2202107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5400235) q[0];
sx q[0];
rz(-2.6271554) q[0];
sx q[0];
rz(3.0995447) q[0];
rz(-pi) q[1];
rz(-1.6613879) q[2];
sx q[2];
rz(-1.633051) q[2];
sx q[2];
rz(-1.7287776) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.073382811) q[1];
sx q[1];
rz(-2.073895) q[1];
sx q[1];
rz(-0.36243172) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16943361) q[3];
sx q[3];
rz(-0.9421351) q[3];
sx q[3];
rz(-0.43600142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23585524) q[2];
sx q[2];
rz(-1.7451655) q[2];
sx q[2];
rz(0.033795707) q[2];
rz(0.77364051) q[3];
sx q[3];
rz(-1.528911) q[3];
sx q[3];
rz(-2.0788367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77699023) q[0];
sx q[0];
rz(-0.62043014) q[0];
sx q[0];
rz(0.34238368) q[0];
rz(-1.7204334) q[1];
sx q[1];
rz(-0.8232638) q[1];
sx q[1];
rz(-1.8470496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6362444) q[0];
sx q[0];
rz(-0.3925792) q[0];
sx q[0];
rz(2.7468203) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3632107) q[2];
sx q[2];
rz(-2.2515544) q[2];
sx q[2];
rz(-1.7964586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0659938) q[1];
sx q[1];
rz(-2.3477049) q[1];
sx q[1];
rz(2.9692279) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7429084) q[3];
sx q[3];
rz(-1.2381319) q[3];
sx q[3];
rz(1.1928967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2932059) q[2];
sx q[2];
rz(-1.7561971) q[2];
sx q[2];
rz(-2.6050341) q[2];
rz(-2.7287591) q[3];
sx q[3];
rz(-0.53240132) q[3];
sx q[3];
rz(-2.0280973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8263016) q[0];
sx q[0];
rz(-1.2696126) q[0];
sx q[0];
rz(1.4971365) q[0];
rz(-0.81028384) q[1];
sx q[1];
rz(-1.5850001) q[1];
sx q[1];
rz(-2.2946045) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6893456) q[0];
sx q[0];
rz(-0.23437491) q[0];
sx q[0];
rz(-2.330392) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55436418) q[2];
sx q[2];
rz(-2.1424166) q[2];
sx q[2];
rz(-1.8272561) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3851953) q[1];
sx q[1];
rz(-2.1902731) q[1];
sx q[1];
rz(1.2379012) q[1];
rz(-1.296792) q[3];
sx q[3];
rz(-1.4782259) q[3];
sx q[3];
rz(-2.0333729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0392796) q[2];
sx q[2];
rz(-1.3156834) q[2];
sx q[2];
rz(2.5073591) q[2];
rz(-0.63747326) q[3];
sx q[3];
rz(-2.0135148) q[3];
sx q[3];
rz(-2.9243961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98916003) q[0];
sx q[0];
rz(-1.9518873) q[0];
sx q[0];
rz(1.963203) q[0];
rz(0.70869008) q[1];
sx q[1];
rz(-1.3460174) q[1];
sx q[1];
rz(-1.8585471) q[1];
rz(-1.218956) q[2];
sx q[2];
rz(-1.9250122) q[2];
sx q[2];
rz(-0.4831947) q[2];
rz(-1.1905963) q[3];
sx q[3];
rz(-2.1722542) q[3];
sx q[3];
rz(0.62361591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
