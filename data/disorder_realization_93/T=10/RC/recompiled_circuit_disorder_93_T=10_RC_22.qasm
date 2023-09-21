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
rz(-1.7927875) q[1];
sx q[1];
rz(-0.92372149) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1169352) q[0];
sx q[0];
rz(-2.5664461) q[0];
sx q[0];
rz(-2.2206109) q[0];
x q[1];
rz(-1.107723) q[2];
sx q[2];
rz(-0.31422868) q[2];
sx q[2];
rz(0.83480922) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87286283) q[1];
sx q[1];
rz(-1.9951207) q[1];
sx q[1];
rz(0.55117589) q[1];
x q[2];
rz(-0.49433319) q[3];
sx q[3];
rz(-1.9088073) q[3];
sx q[3];
rz(2.5288343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.91360056) q[2];
sx q[2];
rz(-1.8929409) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.0733923) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(-1.1967999) q[0];
rz(-0.67990047) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.686036) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.992422) q[0];
sx q[0];
rz(-2.143321) q[0];
sx q[0];
rz(0.19897977) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37462072) q[2];
sx q[2];
rz(-1.6472367) q[2];
sx q[2];
rz(-0.74707109) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5160408) q[1];
sx q[1];
rz(-1.5572773) q[1];
sx q[1];
rz(1.0479755) q[1];
x q[2];
rz(-1.5442113) q[3];
sx q[3];
rz(-1.9823325) q[3];
sx q[3];
rz(-2.6708024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.7896174) q[2];
rz(-2.9591566) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(-0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75333726) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-2.341111) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(-1.172539) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3282933) q[0];
sx q[0];
rz(-2.1345703) q[0];
sx q[0];
rz(-1.9378807) q[0];
rz(-2.795479) q[2];
sx q[2];
rz(-2.3556404) q[2];
sx q[2];
rz(1.9217938) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.064627083) q[1];
sx q[1];
rz(-0.58276999) q[1];
sx q[1];
rz(-1.1175734) q[1];
x q[2];
rz(-1.3451194) q[3];
sx q[3];
rz(-2.884622) q[3];
sx q[3];
rz(0.32303177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0744434) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(2.2303936) q[2];
rz(-2.1905812) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38055414) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(-3.0134841) q[0];
rz(3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(-2.6180843) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9273705) q[0];
sx q[0];
rz(-2.4282051) q[0];
sx q[0];
rz(-0.58332304) q[0];
x q[1];
rz(0.33004327) q[2];
sx q[2];
rz(-1.651262) q[2];
sx q[2];
rz(-2.6835494) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3829141) q[1];
sx q[1];
rz(-0.74712336) q[1];
sx q[1];
rz(1.9488571) q[1];
rz(-0.13682271) q[3];
sx q[3];
rz(-2.1577583) q[3];
sx q[3];
rz(-0.50859355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(2.5775487) q[2];
rz(2.8530252) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(2.2619757) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(-0.99194828) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0569699) q[0];
sx q[0];
rz(-1.3966494) q[0];
sx q[0];
rz(1.4625545) q[0];
rz(-pi) q[1];
x q[1];
rz(3.132658) q[2];
sx q[2];
rz(-1.2679456) q[2];
sx q[2];
rz(2.6645899) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.76584133) q[1];
sx q[1];
rz(-1.5305133) q[1];
sx q[1];
rz(-2.8326616) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81327849) q[3];
sx q[3];
rz(-1.4762029) q[3];
sx q[3];
rz(-2.7300342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-0.27080718) q[2];
rz(0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(-0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145684) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.7927992) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(1.4250925) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3544918) q[0];
sx q[0];
rz(-2.0658501) q[0];
sx q[0];
rz(-1.3160734) q[0];
x q[1];
rz(-1.8259949) q[2];
sx q[2];
rz(-1.4824502) q[2];
sx q[2];
rz(-1.7905854) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52275601) q[1];
sx q[1];
rz(-2.6869876) q[1];
sx q[1];
rz(1.7230117) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2819727) q[3];
sx q[3];
rz(-2.4058127) q[3];
sx q[3];
rz(1.3501292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1340593) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(-0.18051906) q[3];
sx q[3];
rz(-1.6231977) q[3];
sx q[3];
rz(2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(1.58889) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(0.82180506) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9200631) q[0];
sx q[0];
rz(-0.95525817) q[0];
sx q[0];
rz(0.91381844) q[0];
rz(-pi) q[1];
rz(2.6452438) q[2];
sx q[2];
rz(-0.90663547) q[2];
sx q[2];
rz(1.4585444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8824749) q[1];
sx q[1];
rz(-1.3438517) q[1];
sx q[1];
rz(-0.74101733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.20007) q[3];
sx q[3];
rz(-1.0154187) q[3];
sx q[3];
rz(0.45645082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(-0.129536) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(-0.82292557) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(1.8364505) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1730723) q[0];
sx q[0];
rz(-2.1497796) q[0];
sx q[0];
rz(-1.1770244) q[0];
x q[1];
rz(-0.42758503) q[2];
sx q[2];
rz(-2.3901849) q[2];
sx q[2];
rz(-2.4497355) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23755632) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(-2.0380286) q[1];
x q[2];
rz(1.9100902) q[3];
sx q[3];
rz(-1.5929856) q[3];
sx q[3];
rz(2.0458178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0017073) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(1.6513599) q[2];
rz(-1.0772609) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867144) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(-0.3219147) q[0];
rz(1.6053258) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-0.70294356) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1212595) q[0];
sx q[0];
rz(-0.54176211) q[0];
sx q[0];
rz(0.74777491) q[0];
x q[1];
rz(-1.6417575) q[2];
sx q[2];
rz(-0.5103726) q[2];
sx q[2];
rz(0.77529782) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.96303899) q[1];
sx q[1];
rz(-2.5759765) q[1];
sx q[1];
rz(-1.2600684) q[1];
x q[2];
rz(0.95548198) q[3];
sx q[3];
rz(-1.8468879) q[3];
sx q[3];
rz(2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90074173) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(-1.3396324) q[2];
rz(0.30570269) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777578) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(2.6729565) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9194473) q[0];
sx q[0];
rz(-1.7788017) q[0];
sx q[0];
rz(0.39214765) q[0];
rz(-1.956316) q[2];
sx q[2];
rz(-1.8744933) q[2];
sx q[2];
rz(2.7591443) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.48060265) q[1];
sx q[1];
rz(-2.8669679) q[1];
sx q[1];
rz(-1.3165228) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0300794) q[3];
sx q[3];
rz(-1.1772403) q[3];
sx q[3];
rz(2.4126088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.1432077) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(-1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(2.519683) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(-0.86482277) q[2];
sx q[2];
rz(-2.13158) q[2];
sx q[2];
rz(-2.0958015) q[2];
rz(-1.4680396) q[3];
sx q[3];
rz(-0.65410528) q[3];
sx q[3];
rz(0.42412529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];