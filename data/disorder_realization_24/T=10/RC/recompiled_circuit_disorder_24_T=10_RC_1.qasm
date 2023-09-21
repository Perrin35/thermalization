OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5819117) q[0];
sx q[0];
rz(-2.0547325) q[0];
sx q[0];
rz(1.342919) q[0];
rz(7.8328447) q[1];
sx q[1];
rz(3.0631493) q[1];
sx q[1];
rz(6.9258239) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0138577) q[0];
sx q[0];
rz(-2.60175) q[0];
sx q[0];
rz(3.0873469) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8060993) q[2];
sx q[2];
rz(-2.6087084) q[2];
sx q[2];
rz(-2.3032308) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5583615) q[1];
sx q[1];
rz(-1.6915295) q[1];
sx q[1];
rz(-1.6912778) q[1];
x q[2];
rz(2.9896985) q[3];
sx q[3];
rz(-2.9048902) q[3];
sx q[3];
rz(-0.13329245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6671483) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(1.8623964) q[2];
rz(-0.71875087) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(2.8180502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779125) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(2.1610778) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(2.352879) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096743874) q[0];
sx q[0];
rz(-1.6941438) q[0];
sx q[0];
rz(0.11476536) q[0];
rz(1.6025701) q[2];
sx q[2];
rz(-1.3147768) q[2];
sx q[2];
rz(2.1928744) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2930254) q[1];
sx q[1];
rz(-1.2663406) q[1];
sx q[1];
rz(-0.57363631) q[1];
rz(-pi) q[2];
rz(-2.2867145) q[3];
sx q[3];
rz(-2.2023871) q[3];
sx q[3];
rz(0.43254334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7754037) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(0.186084) q[2];
rz(-0.6535334) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2114975) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(2.511456) q[0];
rz(3.0139626) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(-2.4198467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5436514) q[0];
sx q[0];
rz(-1.3184034) q[0];
sx q[0];
rz(-0.6681722) q[0];
rz(-2.6038405) q[2];
sx q[2];
rz(-2.3719412) q[2];
sx q[2];
rz(2.4485181) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4221146) q[1];
sx q[1];
rz(-1.1269224) q[1];
sx q[1];
rz(-1.5281954) q[1];
rz(2.0201163) q[3];
sx q[3];
rz(-2.46393) q[3];
sx q[3];
rz(0.47798702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3815986) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-0.90908137) q[2];
rz(2.6233853) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(2.853493) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58406126) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(-3.0990565) q[0];
rz(-2.361239) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(-0.76400486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1799058) q[0];
sx q[0];
rz(-2.4628203) q[0];
sx q[0];
rz(1.5907445) q[0];
rz(1.4444703) q[2];
sx q[2];
rz(-0.9300803) q[2];
sx q[2];
rz(1.7510406) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82356794) q[1];
sx q[1];
rz(-1.2122224) q[1];
sx q[1];
rz(-1.093868) q[1];
rz(-pi) q[2];
rz(-0.56960168) q[3];
sx q[3];
rz(-0.74865018) q[3];
sx q[3];
rz(2.1514055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(-2.6085473) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(-0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2247291) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(1.2438783) q[0];
rz(2.9149756) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(-2.7382543) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82542244) q[0];
sx q[0];
rz(-0.24992019) q[0];
sx q[0];
rz(2.4106246) q[0];
rz(-pi) q[1];
rz(1.9344159) q[2];
sx q[2];
rz(-2.223613) q[2];
sx q[2];
rz(-1.8967241) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12521872) q[1];
sx q[1];
rz(-2.5902936) q[1];
sx q[1];
rz(-0.2972879) q[1];
x q[2];
rz(-0.63286085) q[3];
sx q[3];
rz(-0.78213464) q[3];
sx q[3];
rz(-2.6485505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.32101813) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(-2.3590951) q[2];
rz(1.1123505) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(-2.6859786) q[0];
rz(-3.0420711) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(-0.0064370357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71953668) q[0];
sx q[0];
rz(-0.7932084) q[0];
sx q[0];
rz(0.35736812) q[0];
rz(-pi) q[1];
rz(-0.10613425) q[2];
sx q[2];
rz(-1.1960293) q[2];
sx q[2];
rz(-0.29667621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9812614) q[1];
sx q[1];
rz(-0.48109522) q[1];
sx q[1];
rz(1.0465924) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75565831) q[3];
sx q[3];
rz(-2.14058) q[3];
sx q[3];
rz(-0.33982402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.36879888) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(0.84645611) q[2];
rz(0.99772292) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(-1.6830106) q[3];
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
rz(-0.098175123) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(3.085882) q[0];
rz(-0.78272351) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(-1.1605211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46471483) q[0];
sx q[0];
rz(-2.5589295) q[0];
sx q[0];
rz(-0.90375264) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2968282) q[2];
sx q[2];
rz(-0.65462199) q[2];
sx q[2];
rz(-0.95915937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50642636) q[1];
sx q[1];
rz(-2.4392358) q[1];
sx q[1];
rz(-0.96794767) q[1];
rz(-pi) q[2];
rz(0.37961827) q[3];
sx q[3];
rz(-1.2109204) q[3];
sx q[3];
rz(-0.66756638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.065585) q[2];
sx q[2];
rz(-0.92131725) q[2];
sx q[2];
rz(-2.356142) q[2];
rz(-0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(-0.41539645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8287559) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(1.1428517) q[0];
rz(-1.8354592) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(0.41608861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7284262) q[0];
sx q[0];
rz(-0.51983716) q[0];
sx q[0];
rz(-1.1632989) q[0];
x q[1];
rz(-2.0390688) q[2];
sx q[2];
rz(-1.4434012) q[2];
sx q[2];
rz(2.1786736) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0033274) q[1];
sx q[1];
rz(-1.4107553) q[1];
sx q[1];
rz(-0.40469594) q[1];
rz(1.7851402) q[3];
sx q[3];
rz(-2.6594901) q[3];
sx q[3];
rz(-1.6899504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0351506) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(0.32315928) q[2];
rz(2.9390826) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(-0.57730738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768196) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(-1.1449822) q[0];
rz(1.9454983) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(-0.51913613) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7767169) q[0];
sx q[0];
rz(-1.4380699) q[0];
sx q[0];
rz(-1.9507292) q[0];
rz(1.3349418) q[2];
sx q[2];
rz(-1.354885) q[2];
sx q[2];
rz(0.73667919) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12431006) q[1];
sx q[1];
rz(-0.78622422) q[1];
sx q[1];
rz(-3.0148274) q[1];
rz(-pi) q[2];
rz(1.395123) q[3];
sx q[3];
rz(-1.8050977) q[3];
sx q[3];
rz(0.3570041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3141994) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(3.1006151) q[2];
rz(2.273902) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7460019) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(-1.6145153) q[0];
rz(-1.4279667) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(1.4987) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7259827) q[0];
sx q[0];
rz(-1.5506396) q[0];
sx q[0];
rz(-1.495342) q[0];
rz(1.8685568) q[2];
sx q[2];
rz(-2.9794663) q[2];
sx q[2];
rz(-2.8043384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1415256) q[1];
sx q[1];
rz(-2.2728517) q[1];
sx q[1];
rz(0.4483923) q[1];
rz(-pi) q[2];
rz(0.20262952) q[3];
sx q[3];
rz(-1.6912795) q[3];
sx q[3];
rz(2.1287624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3830118) q[2];
sx q[2];
rz(-2.0643533) q[2];
sx q[2];
rz(0.14653462) q[2];
rz(-2.3274029) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(-2.7248785) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2232589) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(1.8756443) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(0.11309359) q[2];
sx q[2];
rz(-1.332915) q[2];
sx q[2];
rz(-1.1168196) q[2];
rz(0.90445789) q[3];
sx q[3];
rz(-2.1790128) q[3];
sx q[3];
rz(2.0603767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
