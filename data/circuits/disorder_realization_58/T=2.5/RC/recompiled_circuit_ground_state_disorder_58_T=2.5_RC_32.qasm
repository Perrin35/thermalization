OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.50103203) q[0];
sx q[0];
rz(-2.3409193) q[0];
sx q[0];
rz(-2.9963357) q[0];
rz(2.2486806) q[1];
sx q[1];
rz(-0.4141663) q[1];
sx q[1];
rz(2.034634) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021028886) q[0];
sx q[0];
rz(-1.0146838) q[0];
sx q[0];
rz(-0.70379852) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8584004) q[2];
sx q[2];
rz(-1.257466) q[2];
sx q[2];
rz(-1.4669795) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8095137) q[1];
sx q[1];
rz(-0.65843907) q[1];
sx q[1];
rz(-1.0926985) q[1];
x q[2];
rz(2.9248819) q[3];
sx q[3];
rz(-1.9418849) q[3];
sx q[3];
rz(-0.75066523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5408111) q[2];
sx q[2];
rz(-1.3593707) q[2];
sx q[2];
rz(0.092279807) q[2];
rz(1.4394834) q[3];
sx q[3];
rz(-1.1441792) q[3];
sx q[3];
rz(0.93851844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.0094902078) q[0];
sx q[0];
rz(-1.1630031) q[0];
sx q[0];
rz(0.60390419) q[0];
rz(1.6101135) q[1];
sx q[1];
rz(-1.8096626) q[1];
sx q[1];
rz(1.5276705) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2137215) q[0];
sx q[0];
rz(-1.2564714) q[0];
sx q[0];
rz(-0.95753352) q[0];
x q[1];
rz(-0.98601922) q[2];
sx q[2];
rz(-2.2246317) q[2];
sx q[2];
rz(2.3399835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2774573) q[1];
sx q[1];
rz(-1.2970595) q[1];
sx q[1];
rz(1.5587193) q[1];
x q[2];
rz(1.5701206) q[3];
sx q[3];
rz(-1.6689894) q[3];
sx q[3];
rz(-0.46969603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4332726) q[2];
sx q[2];
rz(-1.6927787) q[2];
sx q[2];
rz(-2.0937008) q[2];
rz(0.6692872) q[3];
sx q[3];
rz(-0.71988121) q[3];
sx q[3];
rz(-1.1650813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79353756) q[0];
sx q[0];
rz(-2.3568643) q[0];
sx q[0];
rz(-2.134557) q[0];
rz(-2.8496565) q[1];
sx q[1];
rz(-2.597229) q[1];
sx q[1];
rz(-0.67063355) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9898674) q[0];
sx q[0];
rz(-1.4833996) q[0];
sx q[0];
rz(0.94957994) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3449981) q[2];
sx q[2];
rz(-1.050794) q[2];
sx q[2];
rz(-2.7869625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2082767) q[1];
sx q[1];
rz(-1.2348935) q[1];
sx q[1];
rz(2.4927054) q[1];
rz(-pi) q[2];
rz(1.8996096) q[3];
sx q[3];
rz(-1.6805263) q[3];
sx q[3];
rz(-1.3789267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85155073) q[2];
sx q[2];
rz(-0.85323492) q[2];
sx q[2];
rz(0.89149371) q[2];
rz(-2.0391035) q[3];
sx q[3];
rz(-0.99415556) q[3];
sx q[3];
rz(-1.7360784) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82146984) q[0];
sx q[0];
rz(-1.7550884) q[0];
sx q[0];
rz(-1.0572877) q[0];
rz(3.140246) q[1];
sx q[1];
rz(-0.82095447) q[1];
sx q[1];
rz(2.3473306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47424537) q[0];
sx q[0];
rz(-2.4794331) q[0];
sx q[0];
rz(3.1116918) q[0];
rz(-pi) q[1];
rz(-0.35425953) q[2];
sx q[2];
rz(-1.2984167) q[2];
sx q[2];
rz(-2.543022) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92137488) q[1];
sx q[1];
rz(-1.1168606) q[1];
sx q[1];
rz(0.80516071) q[1];
rz(-pi) q[2];
rz(0.36093386) q[3];
sx q[3];
rz(-0.31692255) q[3];
sx q[3];
rz(1.8529056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1069676) q[2];
sx q[2];
rz(-0.62109533) q[2];
sx q[2];
rz(-1.0221488) q[2];
rz(1.243783) q[3];
sx q[3];
rz(-0.94977489) q[3];
sx q[3];
rz(-1.7826084) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46955243) q[0];
sx q[0];
rz(-1.7615027) q[0];
sx q[0];
rz(2.688038) q[0];
rz(1.0376616) q[1];
sx q[1];
rz(-1.1353759) q[1];
sx q[1];
rz(-0.79197788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8773408) q[0];
sx q[0];
rz(-1.2928354) q[0];
sx q[0];
rz(-3.0959913) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99292314) q[2];
sx q[2];
rz(-1.9067592) q[2];
sx q[2];
rz(-0.039388933) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9336613) q[1];
sx q[1];
rz(-1.9092036) q[1];
sx q[1];
rz(1.6945356) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7030992) q[3];
sx q[3];
rz(-3.0022394) q[3];
sx q[3];
rz(-2.9840368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3114634) q[2];
sx q[2];
rz(-0.66581231) q[2];
sx q[2];
rz(-2.8361481) q[2];
rz(-0.18181248) q[3];
sx q[3];
rz(-1.5287377) q[3];
sx q[3];
rz(-2.4055433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55564725) q[0];
sx q[0];
rz(-2.3190627) q[0];
sx q[0];
rz(3.0162051) q[0];
rz(-1.5665945) q[1];
sx q[1];
rz(-1.6927203) q[1];
sx q[1];
rz(0.032141846) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9956995) q[0];
sx q[0];
rz(-0.75135485) q[0];
sx q[0];
rz(-0.32352792) q[0];
rz(-pi) q[1];
rz(0.83937593) q[2];
sx q[2];
rz(-1.617086) q[2];
sx q[2];
rz(0.83781017) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74567262) q[1];
sx q[1];
rz(-0.88282871) q[1];
sx q[1];
rz(-0.93979551) q[1];
rz(-0.7456268) q[3];
sx q[3];
rz(-1.485482) q[3];
sx q[3];
rz(1.8861533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0142168) q[2];
sx q[2];
rz(-1.2140423) q[2];
sx q[2];
rz(1.1227013) q[2];
rz(-3.0681916) q[3];
sx q[3];
rz(-2.1843036) q[3];
sx q[3];
rz(-0.61298031) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70925322) q[0];
sx q[0];
rz(-1.4839577) q[0];
sx q[0];
rz(-2.6313229) q[0];
rz(0.36422745) q[1];
sx q[1];
rz(-0.42114708) q[1];
sx q[1];
rz(1.6417004) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9534665) q[0];
sx q[0];
rz(-1.1082442) q[0];
sx q[0];
rz(0.11818399) q[0];
rz(-pi) q[1];
rz(3.0232863) q[2];
sx q[2];
rz(-1.7276754) q[2];
sx q[2];
rz(1.9594693) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3245974) q[1];
sx q[1];
rz(-0.55068586) q[1];
sx q[1];
rz(0.43795069) q[1];
x q[2];
rz(-0.33320547) q[3];
sx q[3];
rz(-2.3764046) q[3];
sx q[3];
rz(0.122034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4437272) q[2];
sx q[2];
rz(-1.5865734) q[2];
sx q[2];
rz(-0.86722428) q[2];
rz(-1.5602268) q[3];
sx q[3];
rz(-2.9428704) q[3];
sx q[3];
rz(-1.8038512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4174058) q[0];
sx q[0];
rz(-3.0751808) q[0];
sx q[0];
rz(-0.41626406) q[0];
rz(-1.9789713) q[1];
sx q[1];
rz(-1.9527718) q[1];
sx q[1];
rz(-2.3847041) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.074919393) q[0];
sx q[0];
rz(-0.52353379) q[0];
sx q[0];
rz(0.96738167) q[0];
rz(-pi) q[1];
rz(0.9932809) q[2];
sx q[2];
rz(-1.2009635) q[2];
sx q[2];
rz(0.12060697) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1637266) q[1];
sx q[1];
rz(-1.9975348) q[1];
sx q[1];
rz(0.20789288) q[1];
x q[2];
rz(-2.5060095) q[3];
sx q[3];
rz(-1.9548237) q[3];
sx q[3];
rz(-1.4366729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4550712) q[2];
sx q[2];
rz(-1.7073809) q[2];
sx q[2];
rz(-1.3580458) q[2];
rz(-0.81651917) q[3];
sx q[3];
rz(-1.9101382) q[3];
sx q[3];
rz(2.9787298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002976) q[0];
sx q[0];
rz(-1.4485899) q[0];
sx q[0];
rz(-1.7342389) q[0];
rz(-2.3963212) q[1];
sx q[1];
rz(-0.73423568) q[1];
sx q[1];
rz(1.5596681) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2360412) q[0];
sx q[0];
rz(-2.2187382) q[0];
sx q[0];
rz(2.7929162) q[0];
rz(-pi) q[1];
rz(-2.4938131) q[2];
sx q[2];
rz(-1.5967007) q[2];
sx q[2];
rz(-1.6960953) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24589495) q[1];
sx q[1];
rz(-1.1283529) q[1];
sx q[1];
rz(0.14296431) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1464983) q[3];
sx q[3];
rz(-1.6213717) q[3];
sx q[3];
rz(-2.3114396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4173296) q[2];
sx q[2];
rz(-1.686692) q[2];
sx q[2];
rz(2.0261185) q[2];
rz(2.8880902) q[3];
sx q[3];
rz(-1.2616254) q[3];
sx q[3];
rz(-0.0089664627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.1215006) q[0];
sx q[0];
rz(-1.2637063) q[0];
sx q[0];
rz(0.76250917) q[0];
rz(-2.1992042) q[1];
sx q[1];
rz(-2.0768879) q[1];
sx q[1];
rz(1.9047838) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4651637) q[0];
sx q[0];
rz(-2.0662413) q[0];
sx q[0];
rz(1.6406607) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57970826) q[2];
sx q[2];
rz(-2.1695608) q[2];
sx q[2];
rz(2.6331462) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8071704) q[1];
sx q[1];
rz(-1.5958734) q[1];
sx q[1];
rz(-1.5958171) q[1];
rz(-pi) q[2];
rz(-1.4051626) q[3];
sx q[3];
rz(-0.40963033) q[3];
sx q[3];
rz(1.7543242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8952055) q[2];
sx q[2];
rz(-3.0463986) q[2];
sx q[2];
rz(0.0072366317) q[2];
rz(-0.17035189) q[3];
sx q[3];
rz(-1.4879613) q[3];
sx q[3];
rz(2.4836704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8728747) q[0];
sx q[0];
rz(-1.7418516) q[0];
sx q[0];
rz(-2.0280784) q[0];
rz(-3.0524104) q[1];
sx q[1];
rz(-1.5166278) q[1];
sx q[1];
rz(1.2022432) q[1];
rz(1.7066649) q[2];
sx q[2];
rz(-1.7262222) q[2];
sx q[2];
rz(-1.8651967) q[2];
rz(0.30059697) q[3];
sx q[3];
rz(-1.7491805) q[3];
sx q[3];
rz(0.94284369) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
