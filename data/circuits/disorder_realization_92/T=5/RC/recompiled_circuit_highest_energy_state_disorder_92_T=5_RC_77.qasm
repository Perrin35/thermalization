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
rz(-1.5795213) q[0];
sx q[0];
rz(-3.002394) q[0];
sx q[0];
rz(-0.31779274) q[0];
rz(0.7572445) q[1];
sx q[1];
rz(-1.5905967) q[1];
sx q[1];
rz(1.6928033) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23148558) q[0];
sx q[0];
rz(-0.39662374) q[0];
sx q[0];
rz(-0.35901944) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7291387) q[2];
sx q[2];
rz(-1.9082365) q[2];
sx q[2];
rz(-0.98357633) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6283012) q[1];
sx q[1];
rz(-0.54524294) q[1];
sx q[1];
rz(2.1907275) q[1];
x q[2];
rz(1.9775544) q[3];
sx q[3];
rz(-1.4881322) q[3];
sx q[3];
rz(2.5661022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.455787) q[2];
sx q[2];
rz(-1.1263584) q[2];
sx q[2];
rz(2.8796999) q[2];
rz(2.453878) q[3];
sx q[3];
rz(-2.3759638) q[3];
sx q[3];
rz(-2.8253637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9460816) q[0];
sx q[0];
rz(-1.2844149) q[0];
sx q[0];
rz(-2.2971357) q[0];
rz(-2.5020245) q[1];
sx q[1];
rz(-0.25194672) q[1];
sx q[1];
rz(1.2676988) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7396916) q[0];
sx q[0];
rz(-0.88205719) q[0];
sx q[0];
rz(-2.9852377) q[0];
rz(-pi) q[1];
rz(2.9466509) q[2];
sx q[2];
rz(-2.4559074) q[2];
sx q[2];
rz(-0.99239381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87593791) q[1];
sx q[1];
rz(-2.3892683) q[1];
sx q[1];
rz(-0.0068412555) q[1];
rz(2.2926125) q[3];
sx q[3];
rz(-2.0343307) q[3];
sx q[3];
rz(3.0914538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0212448) q[2];
sx q[2];
rz(-2.6159888) q[2];
sx q[2];
rz(0.28908238) q[2];
rz(-0.31288475) q[3];
sx q[3];
rz(-2.2009234) q[3];
sx q[3];
rz(0.59278929) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631571) q[0];
sx q[0];
rz(-1.005123) q[0];
sx q[0];
rz(-2.4682755) q[0];
rz(-0.24766573) q[1];
sx q[1];
rz(-2.1238056) q[1];
sx q[1];
rz(0.30390513) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0756098) q[0];
sx q[0];
rz(-2.3572266) q[0];
sx q[0];
rz(0.83215587) q[0];
rz(-pi) q[1];
rz(-2.3788358) q[2];
sx q[2];
rz(-1.0029364) q[2];
sx q[2];
rz(-3.0936846) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0121234) q[1];
sx q[1];
rz(-0.15145603) q[1];
sx q[1];
rz(-1.1726332) q[1];
rz(-pi) q[2];
rz(-1.1506733) q[3];
sx q[3];
rz(-1.3748079) q[3];
sx q[3];
rz(-2.5677057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5725382) q[2];
sx q[2];
rz(-1.7819541) q[2];
sx q[2];
rz(3.1247826) q[2];
rz(2.861764) q[3];
sx q[3];
rz(-2.7357416) q[3];
sx q[3];
rz(0.60389891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0458321) q[0];
sx q[0];
rz(-2.026676) q[0];
sx q[0];
rz(1.8810077) q[0];
rz(2.6450805) q[1];
sx q[1];
rz(-1.1505726) q[1];
sx q[1];
rz(-1.6373434) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27816726) q[0];
sx q[0];
rz(-0.71411944) q[0];
sx q[0];
rz(-0.89366736) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38522933) q[2];
sx q[2];
rz(-0.57621114) q[2];
sx q[2];
rz(-0.50265898) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.28037) q[1];
sx q[1];
rz(-2.5739453) q[1];
sx q[1];
rz(1.2595792) q[1];
x q[2];
rz(-0.50187364) q[3];
sx q[3];
rz(-1.3089367) q[3];
sx q[3];
rz(-3.1002432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0513976) q[2];
sx q[2];
rz(-0.53956705) q[2];
sx q[2];
rz(-0.16057333) q[2];
rz(-1.4422902) q[3];
sx q[3];
rz(-1.620404) q[3];
sx q[3];
rz(3.1062533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2546688) q[0];
sx q[0];
rz(-1.897432) q[0];
sx q[0];
rz(2.1723893) q[0];
rz(1.2508378) q[1];
sx q[1];
rz(-0.67986095) q[1];
sx q[1];
rz(-2.9108237) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8849205) q[0];
sx q[0];
rz(-1.535059) q[0];
sx q[0];
rz(3.1362094) q[0];
rz(-pi) q[1];
rz(1.5577626) q[2];
sx q[2];
rz(-2.763204) q[2];
sx q[2];
rz(0.8410483) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.077037595) q[1];
sx q[1];
rz(-0.63743507) q[1];
sx q[1];
rz(-2.2861479) q[1];
rz(-pi) q[2];
rz(-2.4668815) q[3];
sx q[3];
rz(-2.5541948) q[3];
sx q[3];
rz(2.796668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.9427003) q[2];
sx q[2];
rz(-2.5653699) q[2];
sx q[2];
rz(-0.25299859) q[2];
rz(-1.7871855) q[3];
sx q[3];
rz(-1.1043045) q[3];
sx q[3];
rz(1.3033006) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.323728) q[0];
sx q[0];
rz(-2.4619894) q[0];
sx q[0];
rz(-3.0815109) q[0];
rz(-0.59728638) q[1];
sx q[1];
rz(-2.2164454) q[1];
sx q[1];
rz(-3.0937321) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40210783) q[0];
sx q[0];
rz(-3.0483709) q[0];
sx q[0];
rz(-1.7339105) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20345236) q[2];
sx q[2];
rz(-1.5849539) q[2];
sx q[2];
rz(-2.3490273) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4255449) q[1];
sx q[1];
rz(-2.7486781) q[1];
sx q[1];
rz(-2.1184738) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6512297) q[3];
sx q[3];
rz(-1.8917326) q[3];
sx q[3];
rz(1.3453573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73703274) q[2];
sx q[2];
rz(-0.76475778) q[2];
sx q[2];
rz(0.16727373) q[2];
rz(0.078992756) q[3];
sx q[3];
rz(-1.1827712) q[3];
sx q[3];
rz(2.3110068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9623742) q[0];
sx q[0];
rz(-3.0308864) q[0];
sx q[0];
rz(-3.0141444) q[0];
rz(-2.2178862) q[1];
sx q[1];
rz(-2.5741003) q[1];
sx q[1];
rz(0.73185241) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028721044) q[0];
sx q[0];
rz(-1.0339619) q[0];
sx q[0];
rz(1.7401933) q[0];
rz(-pi) q[1];
rz(-2.5212366) q[2];
sx q[2];
rz(-0.51232238) q[2];
sx q[2];
rz(0.27119766) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73790794) q[1];
sx q[1];
rz(-2.2353454) q[1];
sx q[1];
rz(-2.7079506) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4226888) q[3];
sx q[3];
rz(-1.2205924) q[3];
sx q[3];
rz(-1.4350791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7544516) q[2];
sx q[2];
rz(-2.0595422) q[2];
sx q[2];
rz(0.88073909) q[2];
rz(0.292101) q[3];
sx q[3];
rz(-0.62551522) q[3];
sx q[3];
rz(2.1706799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3683415) q[0];
sx q[0];
rz(-1.0057978) q[0];
sx q[0];
rz(-0.7811501) q[0];
rz(-1.7225522) q[1];
sx q[1];
rz(-2.1731989) q[1];
sx q[1];
rz(-2.9636545) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74465655) q[0];
sx q[0];
rz(-1.8783924) q[0];
sx q[0];
rz(1.164695) q[0];
x q[1];
rz(0.4660999) q[2];
sx q[2];
rz(-2.0498106) q[2];
sx q[2];
rz(-1.9164849) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7922349) q[1];
sx q[1];
rz(-0.48949896) q[1];
sx q[1];
rz(-1.6952747) q[1];
rz(-pi) q[2];
rz(2.8778361) q[3];
sx q[3];
rz(-1.7510771) q[3];
sx q[3];
rz(2.3041185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51131821) q[2];
sx q[2];
rz(-1.7061468) q[2];
sx q[2];
rz(1.1159631) q[2];
rz(-0.80439862) q[3];
sx q[3];
rz(-2.4837327) q[3];
sx q[3];
rz(2.4120954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(0.46111742) q[0];
sx q[0];
rz(-2.456433) q[0];
sx q[0];
rz(2.6035736) q[0];
rz(-2.8021725) q[1];
sx q[1];
rz(-1.598204) q[1];
sx q[1];
rz(3.0982049) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8953066) q[0];
sx q[0];
rz(-0.116838) q[0];
sx q[0];
rz(0.97703345) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0203982) q[2];
sx q[2];
rz(-1.699889) q[2];
sx q[2];
rz(-1.0132524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85203899) q[1];
sx q[1];
rz(-0.75467907) q[1];
sx q[1];
rz(0.11615495) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59915384) q[3];
sx q[3];
rz(-1.2470922) q[3];
sx q[3];
rz(0.97349461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1084719) q[2];
sx q[2];
rz(-2.7816009) q[2];
sx q[2];
rz(-1.6610422) q[2];
rz(-2.8451653) q[3];
sx q[3];
rz(-2.0015643) q[3];
sx q[3];
rz(-0.75291434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66210371) q[0];
sx q[0];
rz(-2.7487553) q[0];
sx q[0];
rz(-1.1641077) q[0];
rz(2.8083943) q[1];
sx q[1];
rz(-1.4771799) q[1];
sx q[1];
rz(-2.3936757) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40953982) q[0];
sx q[0];
rz(-2.4590069) q[0];
sx q[0];
rz(0.48614283) q[0];
x q[1];
rz(-1.2855154) q[2];
sx q[2];
rz(-0.65212265) q[2];
sx q[2];
rz(2.0959321) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1031706) q[1];
sx q[1];
rz(-1.5161247) q[1];
sx q[1];
rz(1.5780539) q[1];
rz(-2.2638948) q[3];
sx q[3];
rz(-1.5355645) q[3];
sx q[3];
rz(2.8352463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36772874) q[2];
sx q[2];
rz(-1.4089971) q[2];
sx q[2];
rz(1.9749953) q[2];
rz(-1.599132) q[3];
sx q[3];
rz(-1.6049933) q[3];
sx q[3];
rz(-2.0994073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.790697) q[0];
sx q[0];
rz(-2.6804374) q[0];
sx q[0];
rz(2.8275369) q[0];
rz(0.028388609) q[1];
sx q[1];
rz(-0.24463618) q[1];
sx q[1];
rz(-1.8520741) q[1];
rz(-0.51043541) q[2];
sx q[2];
rz(-2.107093) q[2];
sx q[2];
rz(0.34278266) q[2];
rz(-2.6131264) q[3];
sx q[3];
rz(-2.2615643) q[3];
sx q[3];
rz(-0.27558358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
