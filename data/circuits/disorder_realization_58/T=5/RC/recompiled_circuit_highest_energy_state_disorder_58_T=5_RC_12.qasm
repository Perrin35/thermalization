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
rz(0.27920029) q[0];
sx q[0];
rz(-1.2503799) q[0];
sx q[0];
rz(-0.00032902349) q[0];
rz(1.1286796) q[1];
sx q[1];
rz(5.1976701) q[1];
sx q[1];
rz(12.357098) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64855382) q[0];
sx q[0];
rz(-1.4161515) q[0];
sx q[0];
rz(-1.8039319) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26165719) q[2];
sx q[2];
rz(-2.5664866) q[2];
sx q[2];
rz(1.9272137) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0297552) q[1];
sx q[1];
rz(-1.2582114) q[1];
sx q[1];
rz(1.8769911) q[1];
rz(1.8162148) q[3];
sx q[3];
rz(-1.3384463) q[3];
sx q[3];
rz(1.3725252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2613232) q[2];
sx q[2];
rz(-1.3618733) q[2];
sx q[2];
rz(0.23639354) q[2];
rz(2.7586625) q[3];
sx q[3];
rz(-1.087944) q[3];
sx q[3];
rz(1.4517763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3799389) q[0];
sx q[0];
rz(-0.28872696) q[0];
sx q[0];
rz(0.54288236) q[0];
rz(-2.701243) q[1];
sx q[1];
rz(-1.1232702) q[1];
sx q[1];
rz(2.0139258) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5701582) q[0];
sx q[0];
rz(-1.2255403) q[0];
sx q[0];
rz(-2.4717114) q[0];
x q[1];
rz(-2.2067546) q[2];
sx q[2];
rz(-1.6976154) q[2];
sx q[2];
rz(0.93709125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4043413) q[1];
sx q[1];
rz(-1.4228428) q[1];
sx q[1];
rz(2.3322991) q[1];
rz(2.4324969) q[3];
sx q[3];
rz(-1.2768687) q[3];
sx q[3];
rz(2.9788168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27309624) q[2];
sx q[2];
rz(-0.79600969) q[2];
sx q[2];
rz(-1.8491171) q[2];
rz(0.83186197) q[3];
sx q[3];
rz(-1.6658733) q[3];
sx q[3];
rz(-2.7003435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6865987) q[0];
sx q[0];
rz(-0.43231493) q[0];
sx q[0];
rz(1.5469714) q[0];
rz(3.1349685) q[1];
sx q[1];
rz(-2.6057656) q[1];
sx q[1];
rz(2.5033902) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2062195) q[0];
sx q[0];
rz(-1.6068335) q[0];
sx q[0];
rz(-1.0956148) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41902991) q[2];
sx q[2];
rz(-0.81522504) q[2];
sx q[2];
rz(0.65302515) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91649721) q[1];
sx q[1];
rz(-1.2079122) q[1];
sx q[1];
rz(1.1919184) q[1];
x q[2];
rz(-0.21535899) q[3];
sx q[3];
rz(-1.0342068) q[3];
sx q[3];
rz(0.99765618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.745445) q[2];
sx q[2];
rz(-0.026693176) q[2];
sx q[2];
rz(-2.2849582) q[2];
rz(-0.73532295) q[3];
sx q[3];
rz(-1.4152799) q[3];
sx q[3];
rz(-2.0935238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51744866) q[0];
sx q[0];
rz(-3.090591) q[0];
sx q[0];
rz(-2.2080102) q[0];
rz(0.75992641) q[1];
sx q[1];
rz(-0.81587195) q[1];
sx q[1];
rz(1.5911721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226993) q[0];
sx q[0];
rz(-0.29445364) q[0];
sx q[0];
rz(0.30979777) q[0];
rz(0.43689219) q[2];
sx q[2];
rz(-1.5073188) q[2];
sx q[2];
rz(1.2530099) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5129169) q[1];
sx q[1];
rz(-0.73386359) q[1];
sx q[1];
rz(0.34725125) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7316667) q[3];
sx q[3];
rz(-1.817559) q[3];
sx q[3];
rz(2.3575229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7237225) q[2];
sx q[2];
rz(-1.3891862) q[2];
sx q[2];
rz(2.1088481) q[2];
rz(-0.51747733) q[3];
sx q[3];
rz(-1.4314707) q[3];
sx q[3];
rz(1.8364068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037769) q[0];
sx q[0];
rz(-0.10194898) q[0];
sx q[0];
rz(-1.8571437) q[0];
rz(2.2880554) q[1];
sx q[1];
rz(-1.4505016) q[1];
sx q[1];
rz(0.4761214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697107) q[0];
sx q[0];
rz(-2.5418542) q[0];
sx q[0];
rz(-2.1246596) q[0];
rz(0.65484859) q[2];
sx q[2];
rz(-2.691948) q[2];
sx q[2];
rz(-2.8596535) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.13083982) q[1];
sx q[1];
rz(-0.71341842) q[1];
sx q[1];
rz(-0.46229191) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6503324) q[3];
sx q[3];
rz(-0.74132468) q[3];
sx q[3];
rz(-1.1363514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4173172) q[2];
sx q[2];
rz(-0.83160916) q[2];
sx q[2];
rz(-1.5403436) q[2];
rz(0.075751461) q[3];
sx q[3];
rz(-1.3906393) q[3];
sx q[3];
rz(0.72995228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1002355) q[0];
sx q[0];
rz(-2.5900109) q[0];
sx q[0];
rz(-1.3496189) q[0];
rz(-1.2760108) q[1];
sx q[1];
rz(-2.3581235) q[1];
sx q[1];
rz(-1.0412201) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7079456) q[0];
sx q[0];
rz(-0.7793847) q[0];
sx q[0];
rz(0.31275605) q[0];
rz(1.5346933) q[2];
sx q[2];
rz(-2.3182333) q[2];
sx q[2];
rz(0.90561101) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7752616) q[1];
sx q[1];
rz(-1.1995254) q[1];
sx q[1];
rz(1.300578) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75821782) q[3];
sx q[3];
rz(-1.3617523) q[3];
sx q[3];
rz(2.2045639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87018806) q[2];
sx q[2];
rz(-1.1689593) q[2];
sx q[2];
rz(-0.99883396) q[2];
rz(-2.6073604) q[3];
sx q[3];
rz(-1.1778919) q[3];
sx q[3];
rz(-1.6629201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6967195) q[0];
sx q[0];
rz(-2.1598926) q[0];
sx q[0];
rz(1.556742) q[0];
rz(1.0818256) q[1];
sx q[1];
rz(-1.4993246) q[1];
sx q[1];
rz(-2.7781442) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7666192) q[0];
sx q[0];
rz(-0.21265499) q[0];
sx q[0];
rz(-0.029092832) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29691445) q[2];
sx q[2];
rz(-1.6726506) q[2];
sx q[2];
rz(0.43019003) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94564964) q[1];
sx q[1];
rz(-0.87313014) q[1];
sx q[1];
rz(-2.5673668) q[1];
x q[2];
rz(-2.9914175) q[3];
sx q[3];
rz(-1.2969103) q[3];
sx q[3];
rz(-1.3824579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58793679) q[2];
sx q[2];
rz(-2.9755972) q[2];
sx q[2];
rz(1.0774277) q[2];
rz(1.59498) q[3];
sx q[3];
rz(-1.3225821) q[3];
sx q[3];
rz(-0.67302978) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3308291) q[0];
sx q[0];
rz(-1.4973233) q[0];
sx q[0];
rz(-2.8610435) q[0];
rz(2.5827017) q[1];
sx q[1];
rz(-0.92550302) q[1];
sx q[1];
rz(1.9299318) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86211156) q[0];
sx q[0];
rz(-2.0923847) q[0];
sx q[0];
rz(-2.305899) q[0];
rz(-pi) q[1];
rz(1.3468387) q[2];
sx q[2];
rz(-1.0887165) q[2];
sx q[2];
rz(-2.659935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16226823) q[1];
sx q[1];
rz(-2.3613902) q[1];
sx q[1];
rz(-3.1092572) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1365355) q[3];
sx q[3];
rz(-1.3727875) q[3];
sx q[3];
rz(-0.34055948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9676548) q[2];
sx q[2];
rz(-1.329198) q[2];
sx q[2];
rz(-0.8052513) q[2];
rz(2.0179181) q[3];
sx q[3];
rz(-2.0011963) q[3];
sx q[3];
rz(2.8203188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5241908) q[0];
sx q[0];
rz(-0.62365714) q[0];
sx q[0];
rz(-2.6362841) q[0];
rz(-2.4303923) q[1];
sx q[1];
rz(-0.86831793) q[1];
sx q[1];
rz(2.0288846) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18234564) q[0];
sx q[0];
rz(-1.8646949) q[0];
sx q[0];
rz(-2.6698851) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0337721) q[2];
sx q[2];
rz(-1.6514213) q[2];
sx q[2];
rz(-0.18641678) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9660652) q[1];
sx q[1];
rz(-1.6972145) q[1];
sx q[1];
rz(2.011622) q[1];
rz(1.6203513) q[3];
sx q[3];
rz(-2.585404) q[3];
sx q[3];
rz(-2.9002528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.50735551) q[2];
sx q[2];
rz(-1.9017744) q[2];
sx q[2];
rz(-2.9061387) q[2];
rz(0.37774751) q[3];
sx q[3];
rz(-2.7538959) q[3];
sx q[3];
rz(-3.0915251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78558755) q[0];
sx q[0];
rz(-1.8536785) q[0];
sx q[0];
rz(-3.0872784) q[0];
rz(0.76447019) q[1];
sx q[1];
rz(-2.6722494) q[1];
sx q[1];
rz(-0.70942172) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.043886) q[0];
sx q[0];
rz(-1.4672252) q[0];
sx q[0];
rz(-1.3243789) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4357752) q[2];
sx q[2];
rz(-1.3121989) q[2];
sx q[2];
rz(1.0279581) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5125138) q[1];
sx q[1];
rz(-0.98328749) q[1];
sx q[1];
rz(1.4329989) q[1];
rz(0.34281369) q[3];
sx q[3];
rz(-1.6926994) q[3];
sx q[3];
rz(1.8203304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8556122) q[2];
sx q[2];
rz(-1.4748272) q[2];
sx q[2];
rz(-3.0484071) q[2];
rz(2.5423673) q[3];
sx q[3];
rz(-2.5803876) q[3];
sx q[3];
rz(-0.78527251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4335321) q[0];
sx q[0];
rz(-0.96994937) q[0];
sx q[0];
rz(2.1920828) q[0];
rz(-3.0769729) q[1];
sx q[1];
rz(-0.043793543) q[1];
sx q[1];
rz(0.66088062) q[1];
rz(2.1223161) q[2];
sx q[2];
rz(-1.5073771) q[2];
sx q[2];
rz(1.3310968) q[2];
rz(-0.23537469) q[3];
sx q[3];
rz(-0.67900758) q[3];
sx q[3];
rz(0.65386299) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
