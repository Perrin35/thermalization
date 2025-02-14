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
rz(0.084963381) q[0];
sx q[0];
rz(3.4439937) q[0];
sx q[0];
rz(6.2511282) q[0];
rz(-0.0078460296) q[1];
sx q[1];
rz(-0.50385952) q[1];
sx q[1];
rz(2.388968) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0467211) q[0];
sx q[0];
rz(-1.3699475) q[0];
sx q[0];
rz(-2.8152391) q[0];
rz(0.80208366) q[2];
sx q[2];
rz(-0.49979106) q[2];
sx q[2];
rz(-2.0715203) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51576534) q[1];
sx q[1];
rz(-1.9099565) q[1];
sx q[1];
rz(1.8372507) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21423046) q[3];
sx q[3];
rz(-1.2204683) q[3];
sx q[3];
rz(0.080834724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.96693119) q[2];
sx q[2];
rz(-1.175468) q[2];
sx q[2];
rz(-0.82007972) q[2];
rz(2.7005633) q[3];
sx q[3];
rz(-2.0377772) q[3];
sx q[3];
rz(-2.5681514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012861982) q[0];
sx q[0];
rz(-2.2282889) q[0];
sx q[0];
rz(2.7313857) q[0];
rz(-2.7606616) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(-1.7832696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8976645) q[0];
sx q[0];
rz(-1.9635734) q[0];
sx q[0];
rz(2.200109) q[0];
x q[1];
rz(-0.83186457) q[2];
sx q[2];
rz(-0.68507776) q[2];
sx q[2];
rz(-2.5545504) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0002901) q[1];
sx q[1];
rz(-1.7385671) q[1];
sx q[1];
rz(-3.0649158) q[1];
x q[2];
rz(2.9866657) q[3];
sx q[3];
rz(-2.2048031) q[3];
sx q[3];
rz(2.4661494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.35839781) q[2];
sx q[2];
rz(-0.40893778) q[2];
sx q[2];
rz(-2.6386293) q[2];
rz(1.8424235) q[3];
sx q[3];
rz(-2.3830569) q[3];
sx q[3];
rz(0.23183307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7314887) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(-0.9758392) q[0];
rz(-2.2052235) q[1];
sx q[1];
rz(-1.5872637) q[1];
sx q[1];
rz(3.0036614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5463211) q[0];
sx q[0];
rz(-1.0965276) q[0];
sx q[0];
rz(2.7300937) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7439026) q[2];
sx q[2];
rz(-1.4639336) q[2];
sx q[2];
rz(-1.7025089) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5974849) q[1];
sx q[1];
rz(-0.710809) q[1];
sx q[1];
rz(-1.7729458) q[1];
x q[2];
rz(1.8432328) q[3];
sx q[3];
rz(-2.0217293) q[3];
sx q[3];
rz(0.48905269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5522573) q[2];
sx q[2];
rz(-1.5312803) q[2];
sx q[2];
rz(1.4790347) q[2];
rz(2.3432483) q[3];
sx q[3];
rz(-0.84404293) q[3];
sx q[3];
rz(0.020309694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8524356) q[0];
sx q[0];
rz(-0.11141369) q[0];
sx q[0];
rz(1.0668466) q[0];
rz(-1.5929219) q[1];
sx q[1];
rz(-1.8916062) q[1];
sx q[1];
rz(2.7222395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28946149) q[0];
sx q[0];
rz(-1.7999819) q[0];
sx q[0];
rz(0.16332345) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52813645) q[2];
sx q[2];
rz(-1.4924876) q[2];
sx q[2];
rz(-2.9700043) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49653445) q[1];
sx q[1];
rz(-2.7784111) q[1];
sx q[1];
rz(-1.3892322) q[1];
rz(-pi) q[2];
rz(0.97030117) q[3];
sx q[3];
rz(-2.5733549) q[3];
sx q[3];
rz(-1.1201536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7833917) q[2];
sx q[2];
rz(-0.50938598) q[2];
sx q[2];
rz(1.5979213) q[2];
rz(1.915043) q[3];
sx q[3];
rz(-1.2164601) q[3];
sx q[3];
rz(-1.8515057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.876038) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(2.3853886) q[0];
rz(-1.2528231) q[1];
sx q[1];
rz(-0.93106657) q[1];
sx q[1];
rz(1.7255712) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9276816) q[0];
sx q[0];
rz(-2.1577303) q[0];
sx q[0];
rz(-1.8326493) q[0];
rz(-pi) q[1];
rz(1.6966693) q[2];
sx q[2];
rz(-1.2294266) q[2];
sx q[2];
rz(1.9370796) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6292343) q[1];
sx q[1];
rz(-2.6027789) q[1];
sx q[1];
rz(3.1278361) q[1];
x q[2];
rz(-2.6719499) q[3];
sx q[3];
rz(-0.40294632) q[3];
sx q[3];
rz(-2.6279891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34117928) q[2];
sx q[2];
rz(-2.4004553) q[2];
sx q[2];
rz(-2.9608012) q[2];
rz(1.4888658) q[3];
sx q[3];
rz(-1.7427665) q[3];
sx q[3];
rz(1.6256049) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7086696) q[0];
sx q[0];
rz(-1.0754508) q[0];
sx q[0];
rz(2.3413626) q[0];
rz(-0.284614) q[1];
sx q[1];
rz(-2.0814643) q[1];
sx q[1];
rz(3.0175993) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6033961) q[0];
sx q[0];
rz(-1.5992224) q[0];
sx q[0];
rz(1.70065) q[0];
rz(1.8058067) q[2];
sx q[2];
rz(-1.5127276) q[2];
sx q[2];
rz(-3.0259759) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.65035076) q[1];
sx q[1];
rz(-1.4891152) q[1];
sx q[1];
rz(-0.30568576) q[1];
rz(-0.37185566) q[3];
sx q[3];
rz(-1.0069435) q[3];
sx q[3];
rz(-1.2363889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2064712) q[2];
sx q[2];
rz(-1.7266885) q[2];
sx q[2];
rz(-0.07621152) q[2];
rz(1.8501836) q[3];
sx q[3];
rz(-2.0468057) q[3];
sx q[3];
rz(-0.61856234) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16977075) q[0];
sx q[0];
rz(-1.562028) q[0];
sx q[0];
rz(-0.75702697) q[0];
rz(0.79611671) q[1];
sx q[1];
rz(-1.7346104) q[1];
sx q[1];
rz(-0.98181358) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039092075) q[0];
sx q[0];
rz(-1.287724) q[0];
sx q[0];
rz(2.3209178) q[0];
x q[1];
rz(0.070008833) q[2];
sx q[2];
rz(-1.4638136) q[2];
sx q[2];
rz(1.9221523) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10985366) q[1];
sx q[1];
rz(-0.64424911) q[1];
sx q[1];
rz(2.2883313) q[1];
x q[2];
rz(2.997274) q[3];
sx q[3];
rz(-1.8486406) q[3];
sx q[3];
rz(-1.2793737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41165274) q[2];
sx q[2];
rz(-2.3306658) q[2];
sx q[2];
rz(0.5438424) q[2];
rz(-0.020921556) q[3];
sx q[3];
rz(-1.7048416) q[3];
sx q[3];
rz(2.3232443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.92886096) q[0];
sx q[0];
rz(-2.2305771) q[0];
sx q[0];
rz(-2.5182356) q[0];
rz(2.8202672) q[1];
sx q[1];
rz(-0.61505452) q[1];
sx q[1];
rz(-2.8772433) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0118222) q[0];
sx q[0];
rz(-1.7826102) q[0];
sx q[0];
rz(2.2925633) q[0];
x q[1];
rz(-1.4521763) q[2];
sx q[2];
rz(-0.33782321) q[2];
sx q[2];
rz(2.0488957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26786727) q[1];
sx q[1];
rz(-1.8494938) q[1];
sx q[1];
rz(0.81588094) q[1];
rz(-0.13295138) q[3];
sx q[3];
rz(-1.8859409) q[3];
sx q[3];
rz(-0.12115762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.136772) q[2];
sx q[2];
rz(-0.68237582) q[2];
sx q[2];
rz(-0.47478673) q[2];
rz(-2.7759806) q[3];
sx q[3];
rz(-2.6545299) q[3];
sx q[3];
rz(2.9584208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30597618) q[0];
sx q[0];
rz(-0.37537471) q[0];
sx q[0];
rz(-0.48582745) q[0];
rz(-2.0756508) q[1];
sx q[1];
rz(-1.7168047) q[1];
sx q[1];
rz(-0.33445439) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6897278) q[0];
sx q[0];
rz(-2.7948927) q[0];
sx q[0];
rz(-2.6348216) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29269258) q[2];
sx q[2];
rz(-2.6862157) q[2];
sx q[2];
rz(0.51560452) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.124467) q[1];
sx q[1];
rz(-1.7199868) q[1];
sx q[1];
rz(-2.6890432) q[1];
rz(2.3452737) q[3];
sx q[3];
rz(-1.2346754) q[3];
sx q[3];
rz(1.4848061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.36755422) q[2];
sx q[2];
rz(-0.84711051) q[2];
sx q[2];
rz(0.96662194) q[2];
rz(1.4878081) q[3];
sx q[3];
rz(-0.32854587) q[3];
sx q[3];
rz(0.070913471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3286572) q[0];
sx q[0];
rz(-2.2870977) q[0];
sx q[0];
rz(-0.27035126) q[0];
rz(-1.5125754) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(-0.62634748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4975166) q[0];
sx q[0];
rz(-0.77267161) q[0];
sx q[0];
rz(2.0324043) q[0];
rz(-2.7105408) q[2];
sx q[2];
rz(-2.1630175) q[2];
sx q[2];
rz(0.86041245) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8737265) q[1];
sx q[1];
rz(-1.6609816) q[1];
sx q[1];
rz(-1.6117881) q[1];
rz(-pi) q[2];
rz(-0.6932859) q[3];
sx q[3];
rz(-0.48637128) q[3];
sx q[3];
rz(2.6457172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7865929) q[2];
sx q[2];
rz(-0.91675106) q[2];
sx q[2];
rz(1.5724486) q[2];
rz(0.7507503) q[3];
sx q[3];
rz(-0.34199491) q[3];
sx q[3];
rz(2.2641838) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56562051) q[0];
sx q[0];
rz(-1.8436057) q[0];
sx q[0];
rz(2.5634503) q[0];
rz(-1.7987953) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(-0.12231355) q[2];
sx q[2];
rz(-1.5708455) q[2];
sx q[2];
rz(1.4932426) q[2];
rz(-0.20082898) q[3];
sx q[3];
rz(-0.26293892) q[3];
sx q[3];
rz(0.50504167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
