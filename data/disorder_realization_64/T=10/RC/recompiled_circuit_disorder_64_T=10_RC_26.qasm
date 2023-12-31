OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(0.064602764) q[0];
sx q[0];
rz(6.3048007) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(-1.3316766) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0927147) q[0];
sx q[0];
rz(-1.9034916) q[0];
sx q[0];
rz(-0.96941745) q[0];
rz(-1.7982499) q[2];
sx q[2];
rz(-1.4706352) q[2];
sx q[2];
rz(1.4129461) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.761258) q[1];
sx q[1];
rz(-2.3654656) q[1];
sx q[1];
rz(2.9556729) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8033923) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(-0.81830922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8109479) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(-0.60418207) q[2];
rz(1.0243246) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(-1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.8169096) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(-1.0634134) q[0];
rz(0.88513199) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(-3.1399472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5981818) q[0];
sx q[0];
rz(-3.0953005) q[0];
sx q[0];
rz(1.8020736) q[0];
rz(2.4992141) q[2];
sx q[2];
rz(-1.7787691) q[2];
sx q[2];
rz(0.013052879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8811223) q[1];
sx q[1];
rz(-0.80658856) q[1];
sx q[1];
rz(-1.9116198) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8257636) q[3];
sx q[3];
rz(-1.4062738) q[3];
sx q[3];
rz(-0.69873519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1559747) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-0.70811159) q[2];
rz(2.0478785) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(-2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22096069) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(-3.1365085) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(1.089383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1339061) q[0];
sx q[0];
rz(-2.0984762) q[0];
sx q[0];
rz(2.6469127) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.543534) q[2];
sx q[2];
rz(-1.2130249) q[2];
sx q[2];
rz(-1.9269112) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1091724) q[1];
sx q[1];
rz(-1.65616) q[1];
sx q[1];
rz(-1.1847772) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5855519) q[3];
sx q[3];
rz(-2.2696113) q[3];
sx q[3];
rz(-1.3211105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0777145) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(2.2568978) q[2];
rz(1.9583154) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(-1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(-0.72682056) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(2.7817536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.115288) q[0];
sx q[0];
rz(-1.7548772) q[0];
sx q[0];
rz(2.227965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6511351) q[2];
sx q[2];
rz(-1.480181) q[2];
sx q[2];
rz(3.0504984) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29405669) q[1];
sx q[1];
rz(-1.9901853) q[1];
sx q[1];
rz(1.9860752) q[1];
x q[2];
rz(2.8344645) q[3];
sx q[3];
rz(-1.2369452) q[3];
sx q[3];
rz(-3.0577554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5741253) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(-2.3763669) q[2];
rz(-2.3848173) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(-2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1446447) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(-1.7255406) q[0];
rz(2.7744746) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(2.1062772) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9177592) q[0];
sx q[0];
rz(-2.1609554) q[0];
sx q[0];
rz(2.607164) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0683124) q[2];
sx q[2];
rz(-0.49256941) q[2];
sx q[2];
rz(-2.5626593) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1977735) q[1];
sx q[1];
rz(-3.0316331) q[1];
sx q[1];
rz(0.52093671) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27541311) q[3];
sx q[3];
rz(-1.3684891) q[3];
sx q[3];
rz(1.3261258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7845903) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(-2.3727097) q[2];
rz(0.33603493) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.4679573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.1795905) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(1.1451716) q[0];
rz(2.0369453) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(0.2072269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6672872) q[0];
sx q[0];
rz(-1.4315839) q[0];
sx q[0];
rz(2.5179203) q[0];
rz(-pi) q[1];
rz(-2.5404262) q[2];
sx q[2];
rz(-1.1784369) q[2];
sx q[2];
rz(0.62077921) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6991899) q[1];
sx q[1];
rz(-1.8132205) q[1];
sx q[1];
rz(1.3139903) q[1];
rz(2.9793671) q[3];
sx q[3];
rz(-2.3174006) q[3];
sx q[3];
rz(-1.970286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3383639) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(-0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8969144) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-2.4095643) q[0];
rz(0.0094982068) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(-2.8498555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4371944) q[0];
sx q[0];
rz(-1.5196374) q[0];
sx q[0];
rz(2.722446) q[0];
rz(0.10642274) q[2];
sx q[2];
rz(-1.8723882) q[2];
sx q[2];
rz(2.0472722) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8818672) q[1];
sx q[1];
rz(-1.7181267) q[1];
sx q[1];
rz(-0.30585652) q[1];
rz(-1.2613867) q[3];
sx q[3];
rz(-1.0687807) q[3];
sx q[3];
rz(-1.3537784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(-2.3042802) q[2];
rz(1.1710179) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432805) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(2.4777381) q[0];
rz(-3.035416) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(0.95867872) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7784568) q[0];
sx q[0];
rz(-2.083337) q[0];
sx q[0];
rz(2.1666359) q[0];
rz(-pi) q[1];
rz(2.0197228) q[2];
sx q[2];
rz(-0.58510963) q[2];
sx q[2];
rz(1.4657071) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2221335) q[1];
sx q[1];
rz(-0.53680116) q[1];
sx q[1];
rz(-1.7807351) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4484552) q[3];
sx q[3];
rz(-2.1515176) q[3];
sx q[3];
rz(-1.8078705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0083996) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(-0.91910249) q[2];
rz(-1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8063426) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(-0.43689716) q[0];
rz(-0.70029744) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(2.6760496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31637329) q[0];
sx q[0];
rz(-1.6557949) q[0];
sx q[0];
rz(0.70789106) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1586458) q[2];
sx q[2];
rz(-1.4530384) q[2];
sx q[2];
rz(-0.6790557) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8852639) q[1];
sx q[1];
rz(-1.1210821) q[1];
sx q[1];
rz(-1.9535669) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97165473) q[3];
sx q[3];
rz(-2.7145436) q[3];
sx q[3];
rz(-3.0640102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(-1.4979866) q[2];
rz(-2.9368029) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(-0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4964504) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(2.0196594) q[0];
rz(0.75138584) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(-2.5591992) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84401417) q[0];
sx q[0];
rz(-1.7739002) q[0];
sx q[0];
rz(-0.88589478) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2920612) q[2];
sx q[2];
rz(-0.50694743) q[2];
sx q[2];
rz(-2.7915733) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6873684) q[1];
sx q[1];
rz(-0.2174938) q[1];
sx q[1];
rz(0.99517676) q[1];
rz(-0.21721812) q[3];
sx q[3];
rz(-1.5981042) q[3];
sx q[3];
rz(0.56009968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1311243) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(2.3790512) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(-1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0762155) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(-1.3394042) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(-2.8600678) q[2];
sx q[2];
rz(-0.83419656) q[2];
sx q[2];
rz(0.4783334) q[2];
rz(-2.3253757) q[3];
sx q[3];
rz(-0.64707884) q[3];
sx q[3];
rz(-0.94221471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
