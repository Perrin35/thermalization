OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7862608) q[0];
sx q[0];
rz(-0.064602764) q[0];
sx q[0];
rz(0.021615418) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(1.8099161) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0776805) q[0];
sx q[0];
rz(-0.67718107) q[0];
sx q[0];
rz(1.022524) q[0];
rz(-pi) q[1];
rz(1.3433427) q[2];
sx q[2];
rz(-1.6709575) q[2];
sx q[2];
rz(1.7286466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6380438) q[1];
sx q[1];
rz(-2.330144) q[1];
sx q[1];
rz(-1.7502977) q[1];
rz(1.1159775) q[3];
sx q[3];
rz(-2.4664719) q[3];
sx q[3];
rz(1.3779373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(0.60418207) q[2];
rz(2.1172681) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169096) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(2.0781793) q[0];
rz(-2.2564607) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(-3.1399472) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7749274) q[0];
sx q[0];
rz(-1.6158551) q[0];
sx q[0];
rz(-0.010618322) q[0];
x q[1];
rz(2.8029289) q[2];
sx q[2];
rz(-2.4709457) q[2];
sx q[2];
rz(1.8530958) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8811223) q[1];
sx q[1];
rz(-2.3350041) q[1];
sx q[1];
rz(1.2299728) q[1];
rz(0.9886338) q[3];
sx q[3];
rz(-0.30246624) q[3];
sx q[3];
rz(-2.8305588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1559747) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22096069) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(0.0050841252) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(-1.089383) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1339061) q[0];
sx q[0];
rz(-2.0984762) q[0];
sx q[0];
rz(0.49467996) q[0];
x q[1];
rz(-0.58615746) q[2];
sx q[2];
rz(-0.68550368) q[2];
sx q[2];
rz(-0.11867487) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5685801) q[1];
sx q[1];
rz(-1.9553361) q[1];
sx q[1];
rz(-0.09210715) q[1];
rz(-pi) q[2];
rz(-0.79077625) q[3];
sx q[3];
rz(-1.1547935) q[3];
sx q[3];
rz(-2.5115867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0638782) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(-0.88469488) q[2];
rz(1.9583154) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(-2.4147721) q[0];
rz(-2.3379393) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(2.7817536) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4041876) q[0];
sx q[0];
rz(-0.92659896) q[0];
sx q[0];
rz(2.9106211) q[0];
x q[1];
rz(0.6511351) q[2];
sx q[2];
rz(-1.480181) q[2];
sx q[2];
rz(0.091094253) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6105146) q[1];
sx q[1];
rz(-0.58137608) q[1];
sx q[1];
rz(2.4060712) q[1];
x q[2];
rz(2.2877341) q[3];
sx q[3];
rz(-0.44970185) q[3];
sx q[3];
rz(-2.2886697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56746733) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(2.3763669) q[2];
rz(-2.3848173) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(-0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1446447) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(1.416052) q[0];
rz(-2.7744746) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(-1.0353154) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9177592) q[0];
sx q[0];
rz(-0.98063722) q[0];
sx q[0];
rz(-2.607164) q[0];
rz(3.0683124) q[2];
sx q[2];
rz(-2.6490232) q[2];
sx q[2];
rz(2.5626593) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6742179) q[1];
sx q[1];
rz(-1.4754703) q[1];
sx q[1];
rz(1.6256888) q[1];
rz(-0.64621332) q[3];
sx q[3];
rz(-2.8013902) q[3];
sx q[3];
rz(0.37341213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7845903) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(0.33603493) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(1.1451716) q[0];
rz(-1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-2.9343658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475537) q[0];
sx q[0];
rz(-2.5045966) q[0];
sx q[0];
rz(0.23547049) q[0];
rz(2.5099498) q[2];
sx q[2];
rz(-0.70438671) q[2];
sx q[2];
rz(0.44142516) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.44240272) q[1];
sx q[1];
rz(-1.8132205) q[1];
sx q[1];
rz(1.8276023) q[1];
rz(-2.3239922) q[3];
sx q[3];
rz(-1.451965) q[3];
sx q[3];
rz(2.6314051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(-0.72171372) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24467829) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(2.4095643) q[0];
rz(-3.1320944) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(0.2917372) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0193298) q[0];
sx q[0];
rz(-0.42207345) q[0];
sx q[0];
rz(-3.0164369) q[0];
rz(1.8740011) q[2];
sx q[2];
rz(-1.6723987) q[2];
sx q[2];
rz(2.6333957) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3953494) q[1];
sx q[1];
rz(-2.8031073) q[1];
sx q[1];
rz(-2.6836718) q[1];
rz(-pi) q[2];
rz(2.6351356) q[3];
sx q[3];
rz(-2.5589057) q[3];
sx q[3];
rz(-1.2014233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(2.3042802) q[2];
rz(1.9705747) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0432805) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(0.6638546) q[0];
rz(0.10617667) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(0.95867872) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52866919) q[0];
sx q[0];
rz(-2.0818424) q[0];
sx q[0];
rz(-2.5445166) q[0];
x q[1];
rz(-2.8616222) q[2];
sx q[2];
rz(-2.0915871) q[2];
sx q[2];
rz(-2.1998646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6764076) q[1];
sx q[1];
rz(-1.0470111) q[1];
sx q[1];
rz(-0.1233867) q[1];
x q[2];
rz(1.6931375) q[3];
sx q[3];
rz(-0.99007505) q[3];
sx q[3];
rz(-1.3337222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0083996) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(0.91910249) q[2];
rz(1.3778) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(-1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063426) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(0.43689716) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(2.6760496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8252194) q[0];
sx q[0];
rz(-1.4857978) q[0];
sx q[0];
rz(-2.4337016) q[0];
rz(-pi) q[1];
rz(3.0131857) q[2];
sx q[2];
rz(-1.1616716) q[2];
sx q[2];
rz(-0.94305925) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25632875) q[1];
sx q[1];
rz(-2.0205106) q[1];
sx q[1];
rz(-1.1880258) q[1];
rz(-1.9302619) q[3];
sx q[3];
rz(-1.8065479) q[3];
sx q[3];
rz(0.93709968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.643606) q[2];
rz(2.9368029) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(-2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64514226) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(2.0196594) q[0];
rz(-2.3902068) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(-0.5823935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2975785) q[0];
sx q[0];
rz(-1.7739002) q[0];
sx q[0];
rz(0.88589478) q[0];
rz(-pi) q[1];
rz(-0.15162823) q[2];
sx q[2];
rz(-2.0564338) q[2];
sx q[2];
rz(3.1079353) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6813587) q[1];
sx q[1];
rz(-1.6885307) q[1];
sx q[1];
rz(1.387499) q[1];
rz(-2.9243745) q[3];
sx q[3];
rz(-1.5434885) q[3];
sx q[3];
rz(-2.581493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1311243) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(2.3790512) q[2];
rz(-1.4108346) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0762155) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(1.3394042) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(-0.81417685) q[2];
sx q[2];
rz(-1.7780751) q[2];
sx q[2];
rz(1.8572394) q[2];
rz(0.4776095) q[3];
sx q[3];
rz(-1.1160679) q[3];
sx q[3];
rz(3.06649) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
