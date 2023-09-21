OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(2.0547325) q[0];
sx q[0];
rz(7.6261043) q[0];
rz(1.5496594) q[1];
sx q[1];
rz(-0.078443371) q[1];
sx q[1];
rz(2.4989541) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064518236) q[0];
sx q[0];
rz(-2.1097578) q[0];
sx q[0];
rz(1.5383188) q[0];
rz(-1.7625916) q[2];
sx q[2];
rz(-1.0704874) q[2];
sx q[2];
rz(-1.2230011) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77036422) q[1];
sx q[1];
rz(-2.9712354) q[1];
sx q[1];
rz(2.360886) q[1];
rz(-pi) q[2];
rz(1.60728) q[3];
sx q[3];
rz(-1.3368703) q[3];
sx q[3];
rz(-2.8521188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.47444433) q[2];
sx q[2];
rz(-2.4953304) q[2];
sx q[2];
rz(-1.2791963) q[2];
rz(0.71875087) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(-2.8180502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779125) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(-0.98051488) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(-2.352879) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65608998) q[0];
sx q[0];
rz(-2.9733109) q[0];
sx q[0];
rz(-0.82505723) q[0];
rz(-pi) q[1];
rz(-2.8854495) q[2];
sx q[2];
rz(-1.6015341) q[2];
sx q[2];
rz(-2.5275633) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2930254) q[1];
sx q[1];
rz(-1.2663406) q[1];
sx q[1];
rz(-2.5679563) q[1];
rz(0.85487811) q[3];
sx q[3];
rz(-0.93920556) q[3];
sx q[3];
rz(2.7090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7754037) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(0.186084) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.7553522) q[3];
sx q[3];
rz(-3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2114975) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(2.511456) q[0];
rz(-3.0139626) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(-0.72174597) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5979413) q[0];
sx q[0];
rz(-1.3184034) q[0];
sx q[0];
rz(0.6681722) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6038405) q[2];
sx q[2];
rz(-2.3719412) q[2];
sx q[2];
rz(0.69307454) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6205412) q[1];
sx q[1];
rz(-0.44577814) q[1];
sx q[1];
rz(0.089322395) q[1];
x q[2];
rz(0.33629041) q[3];
sx q[3];
rz(-2.170917) q[3];
sx q[3];
rz(2.1093413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3815986) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-2.2325113) q[2];
rz(-2.6233853) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(2.853493) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58406126) q[0];
sx q[0];
rz(-0.73636213) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(0.78035367) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(-0.76400486) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9616868) q[0];
sx q[0];
rz(-0.67877239) q[0];
sx q[0];
rz(-1.5907445) q[0];
rz(1.6971223) q[2];
sx q[2];
rz(-0.9300803) q[2];
sx q[2];
rz(1.390552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82356794) q[1];
sx q[1];
rz(-1.9293702) q[1];
sx q[1];
rz(2.0477247) q[1];
rz(-0.66391151) q[3];
sx q[3];
rz(-1.9466562) q[3];
sx q[3];
rz(-1.0192878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(-0.53304535) q[2];
rz(-2.879203) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2247291) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(1.2438783) q[0];
rz(-0.22661701) q[1];
sx q[1];
rz(-2.367327) q[1];
sx q[1];
rz(-0.40333834) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6809083) q[0];
sx q[0];
rz(-1.4049238) q[0];
sx q[0];
rz(0.18780639) q[0];
rz(-1.2071768) q[2];
sx q[2];
rz(-2.223613) q[2];
sx q[2];
rz(-1.8967241) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4407318) q[1];
sx q[1];
rz(-1.7248389) q[1];
sx q[1];
rz(-2.6101019) q[1];
x q[2];
rz(-0.67540695) q[3];
sx q[3];
rz(-1.1408148) q[3];
sx q[3];
rz(-0.59795415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32101813) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(0.78249758) q[2];
rz(1.1123505) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(-2.6859786) q[0];
rz(-0.09952155) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(-0.0064370357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71953668) q[0];
sx q[0];
rz(-0.7932084) q[0];
sx q[0];
rz(0.35736812) q[0];
rz(-pi) q[1];
rz(1.3077277) q[2];
sx q[2];
rz(-2.7527713) q[2];
sx q[2];
rz(3.1281272) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2575063) q[1];
sx q[1];
rz(-1.8045366) q[1];
sx q[1];
rz(1.9952378) q[1];
x q[2];
rz(2.3859343) q[3];
sx q[3];
rz(-1.0010127) q[3];
sx q[3];
rz(-0.33982402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7727938) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(-2.2951365) q[2];
rz(0.99772292) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434175) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(3.085882) q[0];
rz(0.78272351) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.1605211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6877277) q[0];
sx q[0];
rz(-1.9181607) q[0];
sx q[0];
rz(1.0930644) q[0];
rz(-pi) q[1];
rz(0.84476446) q[2];
sx q[2];
rz(-2.4869707) q[2];
sx q[2];
rz(2.1824333) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9142368) q[1];
sx q[1];
rz(-2.1319234) q[1];
sx q[1];
rz(-0.44740541) q[1];
rz(-1.1858995) q[3];
sx q[3];
rz(-1.9249831) q[3];
sx q[3];
rz(-2.3779496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0760076) q[2];
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
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287559) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(1.1428517) q[0];
rz(1.3061334) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(-0.41608861) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6255105) q[0];
sx q[0];
rz(-1.372638) q[0];
sx q[0];
rz(-2.0546191) q[0];
rz(2.9990254) q[2];
sx q[2];
rz(-2.0349742) q[2];
sx q[2];
rz(0.67205059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5059698) q[1];
sx q[1];
rz(-1.1715679) q[1];
sx q[1];
rz(-1.3969621) q[1];
x q[2];
rz(1.7851402) q[3];
sx q[3];
rz(-2.6594901) q[3];
sx q[3];
rz(-1.6899504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1064421) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(0.32315928) q[2];
rz(2.9390826) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(0.57730738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768196) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(-1.1449822) q[0];
rz(-1.1960944) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(2.6224565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7767169) q[0];
sx q[0];
rz(-1.4380699) q[0];
sx q[0];
rz(-1.9507292) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3349418) q[2];
sx q[2];
rz(-1.354885) q[2];
sx q[2];
rz(0.73667919) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0172826) q[1];
sx q[1];
rz(-0.78622422) q[1];
sx q[1];
rz(3.0148274) q[1];
rz(0.63202745) q[3];
sx q[3];
rz(-2.8497189) q[3];
sx q[3];
rz(2.8458418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8273932) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(3.1006151) q[2];
rz(2.273902) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(-2.6303671) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39559078) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(-1.6145153) q[0];
rz(-1.4279667) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(-1.4987) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.98793) q[0];
sx q[0];
rz(-1.4953574) q[0];
sx q[0];
rz(0.020214202) q[0];
rz(-pi) q[1];
x q[1];
rz(0.047949009) q[2];
sx q[2];
rz(-1.7257294) q[2];
sx q[2];
rz(2.5028554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1415256) q[1];
sx q[1];
rz(-2.2728517) q[1];
sx q[1];
rz(2.6932004) q[1];
rz(1.4478217) q[3];
sx q[3];
rz(-1.7719367) q[3];
sx q[3];
rz(-2.5589383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(-2.995058) q[2];
rz(-2.3274029) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(-2.7248785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.2232589) q[0];
sx q[0];
rz(-1.2467361) q[0];
sx q[0];
rz(0.99714739) q[0];
rz(1.2659484) q[1];
sx q[1];
rz(-2.1389778) q[1];
sx q[1];
rz(-1.9139342) q[1];
rz(-1.3314432) q[2];
sx q[2];
rz(-1.6806921) q[2];
sx q[2];
rz(0.42721911) q[2];
rz(0.72600611) q[3];
sx q[3];
rz(-2.2719759) q[3];
sx q[3];
rz(-2.0235973) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
