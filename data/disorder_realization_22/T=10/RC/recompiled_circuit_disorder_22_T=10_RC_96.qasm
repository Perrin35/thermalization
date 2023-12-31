OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(-2.1837283) q[0];
sx q[0];
rz(-0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(-2.1638343) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.190783) q[0];
sx q[0];
rz(-1.8291744) q[0];
sx q[0];
rz(-1.4825312) q[0];
rz(-0.30774967) q[2];
sx q[2];
rz(-2.0648742) q[2];
sx q[2];
rz(0.57927629) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8618776) q[1];
sx q[1];
rz(-1.4289083) q[1];
sx q[1];
rz(-0.14076294) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9631859) q[3];
sx q[3];
rz(-0.79531407) q[3];
sx q[3];
rz(-0.91245302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67291659) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(2.2093175) q[2];
rz(-2.9428234) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(2.1762302) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4345877) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(0.36112753) q[0];
rz(-1.4350285) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(2.3235869) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0571787) q[0];
sx q[0];
rz(-1.1857496) q[0];
sx q[0];
rz(-1.2731228) q[0];
x q[1];
rz(0.12273522) q[2];
sx q[2];
rz(-1.9397748) q[2];
sx q[2];
rz(0.69532794) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9243014) q[1];
sx q[1];
rz(-2.0089307) q[1];
sx q[1];
rz(0.89466722) q[1];
rz(1.9345476) q[3];
sx q[3];
rz(-1.517429) q[3];
sx q[3];
rz(-1.766891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25257418) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(1.7209631) q[2];
rz(-1.8255) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(-2.7533598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4017568) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(-0.87093583) q[0];
rz(2.8254106) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(2.8443764) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39697159) q[0];
sx q[0];
rz(-0.40980761) q[0];
sx q[0];
rz(2.9817392) q[0];
x q[1];
rz(-0.88163968) q[2];
sx q[2];
rz(-0.92703968) q[2];
sx q[2];
rz(-1.2780485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6857022) q[1];
sx q[1];
rz(-0.19280242) q[1];
sx q[1];
rz(1.8222068) q[1];
rz(0.8244332) q[3];
sx q[3];
rz(-1.0472877) q[3];
sx q[3];
rz(1.6214961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3729942) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(-0.42759839) q[2];
rz(1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(-2.6141613) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4363842) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(-2.3676681) q[0];
rz(0.71290839) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(-2.4598222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43142372) q[0];
sx q[0];
rz(-1.3623326) q[0];
sx q[0];
rz(-1.8375977) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21708023) q[2];
sx q[2];
rz(-2.1639369) q[2];
sx q[2];
rz(2.8341688) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7766429) q[1];
sx q[1];
rz(-2.4185191) q[1];
sx q[1];
rz(0.88100453) q[1];
x q[2];
rz(-1.2069615) q[3];
sx q[3];
rz(-1.2483276) q[3];
sx q[3];
rz(0.22842562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1327847) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(2.183389) q[2];
rz(-2.0751674) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(-2.7619894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.5565857) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(2.5812896) q[0];
rz(-2.141748) q[1];
sx q[1];
rz(-0.20345774) q[1];
sx q[1];
rz(1.6220185) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0409659) q[0];
sx q[0];
rz(-1.1275516) q[0];
sx q[0];
rz(-2.6020781) q[0];
x q[1];
rz(2.5298169) q[2];
sx q[2];
rz(-2.0423186) q[2];
sx q[2];
rz(-0.39786354) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2677742) q[1];
sx q[1];
rz(-0.78332892) q[1];
sx q[1];
rz(-0.28247139) q[1];
rz(-pi) q[2];
rz(-0.95789692) q[3];
sx q[3];
rz(-1.6727722) q[3];
sx q[3];
rz(2.0436055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4524298) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(-2.9186644) q[2];
rz(0.034742268) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(-0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3361622) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(-1.0700595) q[0];
rz(-1.7806212) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(1.2840575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43858466) q[0];
sx q[0];
rz(-2.2692338) q[0];
sx q[0];
rz(-0.75011487) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9610923) q[2];
sx q[2];
rz(-1.830606) q[2];
sx q[2];
rz(0.63894546) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9027054) q[1];
sx q[1];
rz(-2.0327912) q[1];
sx q[1];
rz(1.2483031) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0606023) q[3];
sx q[3];
rz(-0.59270699) q[3];
sx q[3];
rz(2.5268775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6665035) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(2.5040023) q[2];
rz(-2.3049138) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5752983) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(-0.56754011) q[0];
rz(-2.7138846) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(2.2033851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0465614) q[0];
sx q[0];
rz(-1.9782269) q[0];
sx q[0];
rz(-0.33336063) q[0];
x q[1];
rz(2.4110255) q[2];
sx q[2];
rz(-2.0148723) q[2];
sx q[2];
rz(-2.4161352) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6066305) q[1];
sx q[1];
rz(-2.3772847) q[1];
sx q[1];
rz(-1.0431837) q[1];
rz(-pi) q[2];
rz(1.2660965) q[3];
sx q[3];
rz(-0.88677553) q[3];
sx q[3];
rz(-1.6001493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2618835) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(-1.3860469) q[2];
rz(1.8188247) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(-0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26738527) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(1.7077131) q[0];
rz(-1.8677615) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(0.83126718) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76038218) q[0];
sx q[0];
rz(-1.4236649) q[0];
sx q[0];
rz(2.7358664) q[0];
rz(-1.3966884) q[2];
sx q[2];
rz(-2.499352) q[2];
sx q[2];
rz(-2.1070534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.23849328) q[1];
sx q[1];
rz(-1.0817173) q[1];
sx q[1];
rz(0.78519435) q[1];
rz(-pi) q[2];
rz(-0.23037489) q[3];
sx q[3];
rz(-1.3532234) q[3];
sx q[3];
rz(2.2090467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5635809) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-2.1772299) q[2];
rz(1.9780805) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37373856) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(-0.97737616) q[0];
rz(-1.7550229) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(-2.0358553) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0375263) q[0];
sx q[0];
rz(-2.6797047) q[0];
sx q[0];
rz(-0.26432963) q[0];
x q[1];
rz(3.0366304) q[2];
sx q[2];
rz(-2.2347921) q[2];
sx q[2];
rz(0.86519372) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3016262) q[1];
sx q[1];
rz(-1.3857538) q[1];
sx q[1];
rz(-1.806083) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68007277) q[3];
sx q[3];
rz(-2.701093) q[3];
sx q[3];
rz(2.3624453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(0.14979714) q[2];
rz(-1.3730565) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49366429) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(3.0143484) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-0.36247411) q[1];
sx q[1];
rz(2.9737934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0509335) q[0];
sx q[0];
rz(-1.5899842) q[0];
sx q[0];
rz(-3.1113383) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4160412) q[2];
sx q[2];
rz(-0.96826474) q[2];
sx q[2];
rz(-0.11520152) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.93712902) q[1];
sx q[1];
rz(-0.58848721) q[1];
sx q[1];
rz(-3.1292533) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4526669) q[3];
sx q[3];
rz(-1.6932634) q[3];
sx q[3];
rz(0.99692217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1148791) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(3.051565) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(0.95054039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5933843) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(0.38800115) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(2.3315196) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(-1.2896982) q[3];
sx q[3];
rz(-2.5847808) q[3];
sx q[3];
rz(1.9980711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
