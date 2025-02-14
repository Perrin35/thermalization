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
rz(-0.64033163) q[0];
sx q[0];
rz(-0.92136541) q[0];
sx q[0];
rz(1.6093572) q[0];
rz(0.20707239) q[1];
sx q[1];
rz(4.3011811) q[1];
sx q[1];
rz(11.094697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4764413) q[0];
sx q[0];
rz(-1.8855321) q[0];
sx q[0];
rz(1.2405618) q[0];
rz(0.81469131) q[2];
sx q[2];
rz(-2.1213795) q[2];
sx q[2];
rz(-1.4514635) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5122657) q[1];
sx q[1];
rz(-1.4719324) q[1];
sx q[1];
rz(-1.8105641) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59572409) q[3];
sx q[3];
rz(-1.4200153) q[3];
sx q[3];
rz(2.8861256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49429911) q[2];
sx q[2];
rz(-0.99227253) q[2];
sx q[2];
rz(-0.54599071) q[2];
rz(-2.2030988) q[3];
sx q[3];
rz(-2.6280845) q[3];
sx q[3];
rz(-0.45309666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0168125) q[0];
sx q[0];
rz(-1.1174959) q[0];
sx q[0];
rz(-0.6413396) q[0];
rz(-1.9524139) q[1];
sx q[1];
rz(-2.7180505) q[1];
sx q[1];
rz(1.10434) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3743312) q[0];
sx q[0];
rz(-2.4924954) q[0];
sx q[0];
rz(0.24477203) q[0];
rz(2.10675) q[2];
sx q[2];
rz(-1.286473) q[2];
sx q[2];
rz(0.51229561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84766372) q[1];
sx q[1];
rz(-2.1184455) q[1];
sx q[1];
rz(0.54137648) q[1];
x q[2];
rz(0.087314815) q[3];
sx q[3];
rz(-1.5305291) q[3];
sx q[3];
rz(-2.1873459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39893338) q[2];
sx q[2];
rz(-1.089596) q[2];
sx q[2];
rz(0.99417865) q[2];
rz(-1.5587156) q[3];
sx q[3];
rz(-0.67928687) q[3];
sx q[3];
rz(2.5905632) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0068552103) q[0];
sx q[0];
rz(-1.0900494) q[0];
sx q[0];
rz(-0.63049522) q[0];
rz(-2.9948803) q[1];
sx q[1];
rz(-1.2598597) q[1];
sx q[1];
rz(2.2781118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7566301) q[0];
sx q[0];
rz(-2.4061086) q[0];
sx q[0];
rz(2.0277268) q[0];
rz(-pi) q[1];
rz(2.9296257) q[2];
sx q[2];
rz(-1.8463928) q[2];
sx q[2];
rz(2.8805062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5456704) q[1];
sx q[1];
rz(-2.0059667) q[1];
sx q[1];
rz(-0.11637139) q[1];
x q[2];
rz(1.5210739) q[3];
sx q[3];
rz(-1.8843302) q[3];
sx q[3];
rz(0.16253838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33519393) q[2];
sx q[2];
rz(-1.4525745) q[2];
sx q[2];
rz(2.2550968) q[2];
rz(0.42996201) q[3];
sx q[3];
rz(-0.49476799) q[3];
sx q[3];
rz(-0.89216843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6524413) q[0];
sx q[0];
rz(-2.8471071) q[0];
sx q[0];
rz(2.6926706) q[0];
rz(-2.4849675) q[1];
sx q[1];
rz(-1.4337599) q[1];
sx q[1];
rz(0.33033672) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9245488) q[0];
sx q[0];
rz(-1.7946436) q[0];
sx q[0];
rz(-1.2428817) q[0];
x q[1];
rz(1.1524989) q[2];
sx q[2];
rz(-1.7798063) q[2];
sx q[2];
rz(-1.542926) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1691087) q[1];
sx q[1];
rz(-2.0044998) q[1];
sx q[1];
rz(0.32167158) q[1];
rz(-pi) q[2];
rz(-0.0090325677) q[3];
sx q[3];
rz(-1.4497644) q[3];
sx q[3];
rz(1.5435404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0535447) q[2];
sx q[2];
rz(-1.6036754) q[2];
sx q[2];
rz(-2.9223082) q[2];
rz(-2.3413279) q[3];
sx q[3];
rz(-2.181668) q[3];
sx q[3];
rz(-0.76103359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6680172) q[0];
sx q[0];
rz(-2.0125084) q[0];
sx q[0];
rz(-1.4833204) q[0];
rz(-2.5738916) q[1];
sx q[1];
rz(-1.9300108) q[1];
sx q[1];
rz(-1.6426881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8627825) q[0];
sx q[0];
rz(-1.7176976) q[0];
sx q[0];
rz(2.6082619) q[0];
rz(1.5597867) q[2];
sx q[2];
rz(-1.6686474) q[2];
sx q[2];
rz(2.7400573) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77563846) q[1];
sx q[1];
rz(-1.3729551) q[1];
sx q[1];
rz(2.5256846) q[1];
rz(1.0077644) q[3];
sx q[3];
rz(-0.22528409) q[3];
sx q[3];
rz(-0.14984261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6296926) q[2];
sx q[2];
rz(-0.7581768) q[2];
sx q[2];
rz(2.8153815) q[2];
rz(-0.1263667) q[3];
sx q[3];
rz(-2.5765403) q[3];
sx q[3];
rz(-0.27157426) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3417974) q[0];
sx q[0];
rz(-1.3038776) q[0];
sx q[0];
rz(-0.52249348) q[0];
rz(0.74784589) q[1];
sx q[1];
rz(-1.6035085) q[1];
sx q[1];
rz(-1.1572256) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25230468) q[0];
sx q[0];
rz(-1.2784908) q[0];
sx q[0];
rz(-2.0415123) q[0];
rz(-2.0339478) q[2];
sx q[2];
rz(-1.2105296) q[2];
sx q[2];
rz(-0.86371326) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95871201) q[1];
sx q[1];
rz(-0.17016115) q[1];
sx q[1];
rz(0.99693701) q[1];
x q[2];
rz(-0.22208235) q[3];
sx q[3];
rz(-0.33470585) q[3];
sx q[3];
rz(3.0396653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.010765643) q[2];
sx q[2];
rz(-0.5126493) q[2];
sx q[2];
rz(1.958468) q[2];
rz(-3.0986687) q[3];
sx q[3];
rz(-1.5019417) q[3];
sx q[3];
rz(2.8973798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6965028) q[0];
sx q[0];
rz(-0.88560167) q[0];
sx q[0];
rz(2.4687299) q[0];
rz(-0.54975763) q[1];
sx q[1];
rz(-1.2943228) q[1];
sx q[1];
rz(-1.4901644) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8880418) q[0];
sx q[0];
rz(-1.7490938) q[0];
sx q[0];
rz(0.56714296) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0668903) q[2];
sx q[2];
rz(-2.0286244) q[2];
sx q[2];
rz(-1.7448448) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71201176) q[1];
sx q[1];
rz(-0.44613555) q[1];
sx q[1];
rz(0.22408102) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5225386) q[3];
sx q[3];
rz(-1.537206) q[3];
sx q[3];
rz(0.03625333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0979536) q[2];
sx q[2];
rz(-2.592228) q[2];
sx q[2];
rz(-1.9926386) q[2];
rz(-2.7502381) q[3];
sx q[3];
rz(-1.9214182) q[3];
sx q[3];
rz(-2.2933188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5452165) q[0];
sx q[0];
rz(-1.2059809) q[0];
sx q[0];
rz(2.4429831) q[0];
rz(0.85083234) q[1];
sx q[1];
rz(-0.60092503) q[1];
sx q[1];
rz(2.1934379) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20336172) q[0];
sx q[0];
rz(-1.4885509) q[0];
sx q[0];
rz(-2.4191678) q[0];
rz(-pi) q[1];
rz(-1.7265375) q[2];
sx q[2];
rz(-0.33340764) q[2];
sx q[2];
rz(1.6108212) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8906456) q[1];
sx q[1];
rz(-1.0140103) q[1];
sx q[1];
rz(-2.6614038) q[1];
rz(-0.83017577) q[3];
sx q[3];
rz(-2.1599033) q[3];
sx q[3];
rz(-1.0419221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.66905388) q[2];
sx q[2];
rz(-1.4488181) q[2];
sx q[2];
rz(-1.6834458) q[2];
rz(1.091188) q[3];
sx q[3];
rz(-0.89717054) q[3];
sx q[3];
rz(0.18479656) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882196) q[0];
sx q[0];
rz(-0.18167697) q[0];
sx q[0];
rz(-1.3004119) q[0];
rz(-0.45937195) q[1];
sx q[1];
rz(-1.6875024) q[1];
sx q[1];
rz(-0.5836817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65410173) q[0];
sx q[0];
rz(-2.0493453) q[0];
sx q[0];
rz(0.74639456) q[0];
x q[1];
rz(1.5351717) q[2];
sx q[2];
rz(-1.6093264) q[2];
sx q[2];
rz(-0.57254475) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4897399) q[1];
sx q[1];
rz(-2.0803343) q[1];
sx q[1];
rz(-2.7368828) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0784057) q[3];
sx q[3];
rz(-0.48590966) q[3];
sx q[3];
rz(2.8775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61648458) q[2];
sx q[2];
rz(-0.27518347) q[2];
sx q[2];
rz(0.52213651) q[2];
rz(0.76256847) q[3];
sx q[3];
rz(-1.6430166) q[3];
sx q[3];
rz(2.7975119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6793215) q[0];
sx q[0];
rz(-1.9336047) q[0];
sx q[0];
rz(-2.5479877) q[0];
rz(2.4484334) q[1];
sx q[1];
rz(-0.79265541) q[1];
sx q[1];
rz(2.3902182) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63735139) q[0];
sx q[0];
rz(-0.48287409) q[0];
sx q[0];
rz(0.58321799) q[0];
x q[1];
rz(-0.85752731) q[2];
sx q[2];
rz(-0.81999841) q[2];
sx q[2];
rz(1.6423026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9819465) q[1];
sx q[1];
rz(-1.6345191) q[1];
sx q[1];
rz(-0.66566408) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4937711) q[3];
sx q[3];
rz(-1.4024156) q[3];
sx q[3];
rz(1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0005325) q[2];
sx q[2];
rz(-2.4324721) q[2];
sx q[2];
rz(-0.023690311) q[2];
rz(2.0247816) q[3];
sx q[3];
rz(-1.6938554) q[3];
sx q[3];
rz(-1.1344596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5524207) q[0];
sx q[0];
rz(-1.9235274) q[0];
sx q[0];
rz(-2.0148475) q[0];
rz(1.3202271) q[1];
sx q[1];
rz(-1.8658493) q[1];
sx q[1];
rz(1.0125926) q[1];
rz(-1.3158952) q[2];
sx q[2];
rz(-1.0409689) q[2];
sx q[2];
rz(2.5036176) q[2];
rz(1.8018525) q[3];
sx q[3];
rz(-1.7920296) q[3];
sx q[3];
rz(-1.5248004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
