OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(1.8811037) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0215065) q[0];
sx q[0];
rz(-1.9061631) q[0];
sx q[0];
rz(1.0943227) q[0];
rz(-pi) q[1];
x q[1];
rz(1.853763) q[2];
sx q[2];
rz(-1.4322865) q[2];
sx q[2];
rz(-1.1793009) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0337861) q[1];
sx q[1];
rz(-2.4596655) q[1];
sx q[1];
rz(0.71180196) q[1];
rz(-2.5039623) q[3];
sx q[3];
rz(-0.59083592) q[3];
sx q[3];
rz(1.5096111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91360056) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(0.16201924) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0733923) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(-2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.686036) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992422) q[0];
sx q[0];
rz(-2.143321) q[0];
sx q[0];
rz(-0.19897977) q[0];
rz(2.935264) q[2];
sx q[2];
rz(-0.38197877) q[2];
sx q[2];
rz(-1.0155592) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.96869722) q[1];
sx q[1];
rz(-2.6186133) q[1];
sx q[1];
rz(-1.5978659) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7299269) q[3];
sx q[3];
rz(-1.5951612) q[3];
sx q[3];
rz(-1.0893694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.3519752) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(-2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882554) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(0.80048168) q[0];
rz(0.02877409) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(1.172539) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4370255) q[0];
sx q[0];
rz(-0.66172681) q[0];
sx q[0];
rz(0.5163124) q[0];
rz(-1.2433979) q[2];
sx q[2];
rz(-0.84257579) q[2];
sx q[2];
rz(2.3936405) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.4634358) q[1];
sx q[1];
rz(-2.0883745) q[1];
sx q[1];
rz(2.8606158) q[1];
x q[2];
rz(0.058733744) q[3];
sx q[3];
rz(-1.8211094) q[3];
sx q[3];
rz(-0.55610031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0744434) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(-2.2303936) q[2];
rz(2.1905812) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(2.2495911) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(0.12810853) q[0];
rz(-3.065486) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(0.52350837) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2098171) q[0];
sx q[0];
rz(-2.1486001) q[0];
sx q[0];
rz(2.0156167) q[0];
rz(2.8977215) q[2];
sx q[2];
rz(-2.802231) q[2];
sx q[2];
rz(1.3432168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.254863) q[1];
sx q[1];
rz(-2.2543395) q[1];
sx q[1];
rz(0.32943326) q[1];
rz(-pi) q[2];
rz(0.97949667) q[3];
sx q[3];
rz(-1.6846091) q[3];
sx q[3];
rz(2.1554961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(0.28856746) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-0.99194828) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6447727) q[0];
sx q[0];
rz(-0.20475514) q[0];
sx q[0];
rz(0.55069189) q[0];
x q[1];
rz(-3.132658) q[2];
sx q[2];
rz(-1.2679456) q[2];
sx q[2];
rz(0.47700275) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3757513) q[1];
sx q[1];
rz(-1.5305133) q[1];
sx q[1];
rz(2.8326616) q[1];
rz(1.4335853) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(-1.8828132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(0.27080718) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4145684) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(-2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(1.7165002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3544918) q[0];
sx q[0];
rz(-1.0757425) q[0];
sx q[0];
rz(-1.3160734) q[0];
rz(-pi) q[1];
rz(-1.9082597) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(0.5459107) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7879494) q[1];
sx q[1];
rz(-2.0197581) q[1];
sx q[1];
rz(-3.06762) q[1];
rz(-2.607843) q[3];
sx q[3];
rz(-2.1042049) q[3];
sx q[3];
rz(0.48983869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(2.9610736) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5086223) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-2.3197876) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99188995) q[0];
sx q[0];
rz(-0.8677965) q[0];
sx q[0];
rz(0.71233149) q[0];
x q[1];
rz(-2.1173382) q[2];
sx q[2];
rz(-0.80596906) q[2];
sx q[2];
rz(0.96217996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2591178) q[1];
sx q[1];
rz(-1.7977409) q[1];
sx q[1];
rz(-2.4005753) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6131367) q[3];
sx q[3];
rz(-2.4848357) q[3];
sx q[3];
rz(0.17880759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3372779) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(0.24469963) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(-1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725175) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(-0.82292557) q[0];
rz(0.30934632) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(-1.8364505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5236429) q[0];
sx q[0];
rz(-2.4542913) q[0];
sx q[0];
rz(-2.6108517) q[0];
rz(-pi) q[1];
rz(-1.201198) q[2];
sx q[2];
rz(-0.9005138) q[2];
sx q[2];
rz(3.0073462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.23755632) q[1];
sx q[1];
rz(-0.95103969) q[1];
sx q[1];
rz(1.1035641) q[1];
x q[2];
rz(-0.023530258) q[3];
sx q[3];
rz(-1.9100034) q[3];
sx q[3];
rz(0.48285218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(-3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-1.7548783) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(1.6053258) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(2.4386491) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2962869) q[0];
sx q[0];
rz(-1.1830813) q[0];
sx q[0];
rz(-1.182343) q[0];
x q[1];
rz(1.6417575) q[2];
sx q[2];
rz(-0.5103726) q[2];
sx q[2];
rz(2.3662948) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1785537) q[1];
sx q[1];
rz(-0.56561618) q[1];
sx q[1];
rz(1.8815243) q[1];
rz(0.95548198) q[3];
sx q[3];
rz(-1.8468879) q[3];
sx q[3];
rz(2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90074173) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(-2.83589) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777578) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.9158069) q[0];
rz(-2.2380791) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43384957) q[0];
sx q[0];
rz(-1.9540457) q[0];
sx q[0];
rz(-1.3462523) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.956316) q[2];
sx q[2];
rz(-1.8744933) q[2];
sx q[2];
rz(-0.38244837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.66099) q[1];
sx q[1];
rz(-2.8669679) q[1];
sx q[1];
rz(1.3165228) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1115132) q[3];
sx q[3];
rz(-1.1772403) q[3];
sx q[3];
rz(2.4126088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(-1.1432077) q[2];
rz(-0.11463595) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464012) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(-0.62190965) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(0.86482277) q[2];
sx q[2];
rz(-1.0100126) q[2];
sx q[2];
rz(1.0457912) q[2];
rz(1.4680396) q[3];
sx q[3];
rz(-2.4874874) q[3];
sx q[3];
rz(-2.7174674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];