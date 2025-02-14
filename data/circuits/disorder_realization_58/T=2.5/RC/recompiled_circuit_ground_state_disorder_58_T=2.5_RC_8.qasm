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
rz(-0.89291209) q[1];
sx q[1];
rz(3.555759) q[1];
sx q[1];
rz(10.531737) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1482077) q[0];
sx q[0];
rz(-2.275106) q[0];
sx q[0];
rz(2.3763477) q[0];
rz(-0.85914454) q[2];
sx q[2];
rz(-0.41918735) q[2];
sx q[2];
rz(-0.91793467) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9913276) q[1];
sx q[1];
rz(-1.8561761) q[1];
sx q[1];
rz(2.1726514) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1916394) q[3];
sx q[3];
rz(-1.772545) q[3];
sx q[3];
rz(2.4011322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60078159) q[2];
sx q[2];
rz(-1.3593707) q[2];
sx q[2];
rz(-3.0493128) q[2];
rz(1.7021092) q[3];
sx q[3];
rz(-1.1441792) q[3];
sx q[3];
rz(2.2030742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0094902078) q[0];
sx q[0];
rz(-1.9785896) q[0];
sx q[0];
rz(0.60390419) q[0];
rz(-1.6101135) q[1];
sx q[1];
rz(-1.3319301) q[1];
sx q[1];
rz(1.5276705) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.912589) q[0];
sx q[0];
rz(-0.67978501) q[0];
sx q[0];
rz(-1.0566105) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74320044) q[2];
sx q[2];
rz(-2.024352) q[2];
sx q[2];
rz(-1.1519866) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4316667) q[1];
sx q[1];
rz(-1.559169) q[1];
sx q[1];
rz(2.8678368) q[1];
x q[2];
rz(-0.0068596938) q[3];
sx q[3];
rz(-0.09819542) q[3];
sx q[3];
rz(2.6650037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.70832002) q[2];
sx q[2];
rz(-1.4488139) q[2];
sx q[2];
rz(-2.0937008) q[2];
rz(-2.4723054) q[3];
sx q[3];
rz(-2.4217114) q[3];
sx q[3];
rz(-1.9765114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79353756) q[0];
sx q[0];
rz(-2.3568643) q[0];
sx q[0];
rz(-1.0070356) q[0];
rz(-2.8496565) q[1];
sx q[1];
rz(-0.54436362) q[1];
sx q[1];
rz(0.67063355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1517253) q[0];
sx q[0];
rz(-1.6581931) q[0];
sx q[0];
rz(-0.94957994) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67523662) q[2];
sx q[2];
rz(-0.91885447) q[2];
sx q[2];
rz(-0.76387954) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.047453316) q[1];
sx q[1];
rz(-0.71935868) q[1];
sx q[1];
rz(0.5238976) q[1];
x q[2];
rz(-1.241983) q[3];
sx q[3];
rz(-1.4610664) q[3];
sx q[3];
rz(1.3789267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85155073) q[2];
sx q[2];
rz(-2.2883577) q[2];
sx q[2];
rz(-0.89149371) q[2];
rz(2.0391035) q[3];
sx q[3];
rz(-0.99415556) q[3];
sx q[3];
rz(1.7360784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3201228) q[0];
sx q[0];
rz(-1.7550884) q[0];
sx q[0];
rz(1.0572877) q[0];
rz(-0.0013466324) q[1];
sx q[1];
rz(-0.82095447) q[1];
sx q[1];
rz(2.3473306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7052536) q[0];
sx q[0];
rz(-2.2326075) q[0];
sx q[0];
rz(-1.5940985) q[0];
x q[1];
rz(1.2813454) q[2];
sx q[2];
rz(-1.2301386) q[2];
sx q[2];
rz(-2.2685426) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0932622) q[1];
sx q[1];
rz(-2.2431313) q[1];
sx q[1];
rz(-0.59497084) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8438679) q[3];
sx q[3];
rz(-1.4605165) q[3];
sx q[3];
rz(3.0793234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0346251) q[2];
sx q[2];
rz(-2.5204973) q[2];
sx q[2];
rz(1.0221488) q[2];
rz(-1.8978097) q[3];
sx q[3];
rz(-0.94977489) q[3];
sx q[3];
rz(-1.7826084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46955243) q[0];
sx q[0];
rz(-1.38009) q[0];
sx q[0];
rz(-0.45355466) q[0];
rz(2.1039311) q[1];
sx q[1];
rz(-2.0062168) q[1];
sx q[1];
rz(-0.79197788) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4290473) q[0];
sx q[0];
rz(-0.28158108) q[0];
sx q[0];
rz(1.4123807) q[0];
rz(-2.7466082) q[2];
sx q[2];
rz(-2.1126267) q[2];
sx q[2];
rz(-1.8219558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8199948) q[1];
sx q[1];
rz(-1.6874847) q[1];
sx q[1];
rz(0.34081809) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4326454) q[3];
sx q[3];
rz(-1.589121) q[3];
sx q[3];
rz(1.5973171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3114634) q[2];
sx q[2];
rz(-0.66581231) q[2];
sx q[2];
rz(0.30544454) q[2];
rz(-0.18181248) q[3];
sx q[3];
rz(-1.5287377) q[3];
sx q[3];
rz(-2.4055433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(-1.5749982) q[1];
sx q[1];
rz(-1.4488723) q[1];
sx q[1];
rz(0.032141846) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18460546) q[0];
sx q[0];
rz(-1.3520387) q[0];
sx q[0];
rz(-0.72483988) q[0];
x q[1];
rz(-1.6400385) q[2];
sx q[2];
rz(-0.73261315) q[2];
sx q[2];
rz(-2.46012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6010203) q[1];
sx q[1];
rz(-2.2443612) q[1];
sx q[1];
rz(-2.5190398) q[1];
rz(-pi) q[2];
rz(2.3959659) q[3];
sx q[3];
rz(-1.6561106) q[3];
sx q[3];
rz(-1.8861533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0142168) q[2];
sx q[2];
rz(-1.9275503) q[2];
sx q[2];
rz(2.0188913) q[2];
rz(3.0681916) q[3];
sx q[3];
rz(-2.1843036) q[3];
sx q[3];
rz(0.61298031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4323394) q[0];
sx q[0];
rz(-1.4839577) q[0];
sx q[0];
rz(-0.51026979) q[0];
rz(-0.36422745) q[1];
sx q[1];
rz(-2.7204456) q[1];
sx q[1];
rz(1.6417004) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4356053) q[0];
sx q[0];
rz(-1.4650808) q[0];
sx q[0];
rz(2.0361498) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0232863) q[2];
sx q[2];
rz(-1.7276754) q[2];
sx q[2];
rz(1.1821234) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.31452306) q[1];
sx q[1];
rz(-2.0644958) q[1];
sx q[1];
rz(-1.3160454) q[1];
x q[2];
rz(-2.4046201) q[3];
sx q[3];
rz(-1.342257) q[3];
sx q[3];
rz(-1.6933954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.69786543) q[2];
sx q[2];
rz(-1.5550193) q[2];
sx q[2];
rz(-0.86722428) q[2];
rz(1.5813658) q[3];
sx q[3];
rz(-0.19872228) q[3];
sx q[3];
rz(-1.3377415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7241868) q[0];
sx q[0];
rz(-3.0751808) q[0];
sx q[0];
rz(2.7253286) q[0];
rz(-1.1626214) q[1];
sx q[1];
rz(-1.1888209) q[1];
sx q[1];
rz(-2.3847041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.074919393) q[0];
sx q[0];
rz(-2.6180589) q[0];
sx q[0];
rz(2.174211) q[0];
rz(-pi) q[1];
rz(-2.1882638) q[2];
sx q[2];
rz(-2.4673415) q[2];
sx q[2];
rz(1.9566388) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97786602) q[1];
sx q[1];
rz(-1.9975348) q[1];
sx q[1];
rz(2.9336998) q[1];
rz(-pi) q[2];
rz(-0.63558319) q[3];
sx q[3];
rz(-1.9548237) q[3];
sx q[3];
rz(1.4366729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4550712) q[2];
sx q[2];
rz(-1.7073809) q[2];
sx q[2];
rz(1.3580458) q[2];
rz(0.81651917) q[3];
sx q[3];
rz(-1.2314545) q[3];
sx q[3];
rz(-0.16286287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74129504) q[0];
sx q[0];
rz(-1.4485899) q[0];
sx q[0];
rz(1.7342389) q[0];
rz(2.3963212) q[1];
sx q[1];
rz(-0.73423568) q[1];
sx q[1];
rz(1.5819246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9055515) q[0];
sx q[0];
rz(-2.2187382) q[0];
sx q[0];
rz(-0.34867649) q[0];
rz(-0.042912622) q[2];
sx q[2];
rz(-2.4933698) q[2];
sx q[2];
rz(0.15951482) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.078439039) q[1];
sx q[1];
rz(-2.6780824) q[1];
sx q[1];
rz(-1.2786675) q[1];
x q[2];
rz(0.060272597) q[3];
sx q[3];
rz(-2.1456686) q[3];
sx q[3];
rz(0.70784345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72426307) q[2];
sx q[2];
rz(-1.686692) q[2];
sx q[2];
rz(-1.1154741) q[2];
rz(-2.8880902) q[3];
sx q[3];
rz(-1.8799672) q[3];
sx q[3];
rz(-0.0089664627) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1215006) q[0];
sx q[0];
rz(-1.8778863) q[0];
sx q[0];
rz(-0.76250917) q[0];
rz(-2.1992042) q[1];
sx q[1];
rz(-2.0768879) q[1];
sx q[1];
rz(1.9047838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4651637) q[0];
sx q[0];
rz(-1.0753514) q[0];
sx q[0];
rz(-1.5009319) q[0];
rz(2.2549596) q[2];
sx q[2];
rz(-2.0403452) q[2];
sx q[2];
rz(0.70882112) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2357465) q[1];
sx q[1];
rz(-1.5958092) q[1];
sx q[1];
rz(3.1165078) q[1];
rz(-1.1661766) q[3];
sx q[3];
rz(-1.6365094) q[3];
sx q[3];
rz(3.1102151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8952055) q[2];
sx q[2];
rz(-0.095194101) q[2];
sx q[2];
rz(-3.134356) q[2];
rz(-2.9712408) q[3];
sx q[3];
rz(-1.6536313) q[3];
sx q[3];
rz(2.4836704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2687179) q[0];
sx q[0];
rz(-1.7418516) q[0];
sx q[0];
rz(-2.0280784) q[0];
rz(-3.0524104) q[1];
sx q[1];
rz(-1.5166278) q[1];
sx q[1];
rz(1.2022432) q[1];
rz(-0.15684814) q[2];
sx q[2];
rz(-1.7050171) q[2];
sx q[2];
rz(-0.27324054) q[2];
rz(0.54696541) q[3];
sx q[3];
rz(-2.7934358) q[3];
sx q[3];
rz(-1.1478333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
