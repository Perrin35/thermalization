OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0260789) q[0];
sx q[0];
rz(-1.6576515) q[0];
sx q[0];
rz(-2.8154362) q[0];
rz(-1.1905319) q[1];
sx q[1];
rz(-1.3500554) q[1];
sx q[1];
rz(-1.5989369) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9261949) q[0];
sx q[0];
rz(-0.30456671) q[0];
sx q[0];
rz(0.79415168) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36122303) q[2];
sx q[2];
rz(-0.63399705) q[2];
sx q[2];
rz(1.0499954) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.029763873) q[1];
sx q[1];
rz(-1.3297237) q[1];
sx q[1];
rz(1.3328711) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78674973) q[3];
sx q[3];
rz(-1.5610715) q[3];
sx q[3];
rz(-1.7973289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2797543) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(-0.88511434) q[2];
rz(-0.72201133) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(-0.0074145934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136114) q[0];
sx q[0];
rz(-0.95887029) q[0];
sx q[0];
rz(2.0425178) q[0];
rz(-0.66501578) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(2.2639993) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646334) q[0];
sx q[0];
rz(-1.3898802) q[0];
sx q[0];
rz(2.4448256) q[0];
x q[1];
rz(0.36466937) q[2];
sx q[2];
rz(-1.696535) q[2];
sx q[2];
rz(-1.0276577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5223915) q[1];
sx q[1];
rz(-3.0295105) q[1];
sx q[1];
rz(1.8735318) q[1];
rz(1.9545994) q[3];
sx q[3];
rz(-1.2591397) q[3];
sx q[3];
rz(2.4812738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7426804) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(2.0111283) q[2];
rz(-1.2997262) q[3];
sx q[3];
rz(-1.2239417) q[3];
sx q[3];
rz(-1.6931504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.995342) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(-0.28999844) q[0];
rz(-0.6668123) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(0.07382948) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683257) q[0];
sx q[0];
rz(-3.0330015) q[0];
sx q[0];
rz(-2.5383839) q[0];
rz(2.8507112) q[2];
sx q[2];
rz(-1.8996432) q[2];
sx q[2];
rz(-1.3128624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5362629) q[1];
sx q[1];
rz(-1.7483835) q[1];
sx q[1];
rz(-3.0847343) q[1];
rz(-pi) q[2];
rz(2.79014) q[3];
sx q[3];
rz(-1.6822527) q[3];
sx q[3];
rz(0.66672882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4905711) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(-0.64669615) q[2];
rz(-2.0329287) q[3];
sx q[3];
rz(-2.3587148) q[3];
sx q[3];
rz(2.0126608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17764238) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(1.9529163) q[0];
rz(2.1229318) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(1.4368988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3387412) q[0];
sx q[0];
rz(-0.84206284) q[0];
sx q[0];
rz(1.4472423) q[0];
x q[1];
rz(2.3827299) q[2];
sx q[2];
rz(-2.7042411) q[2];
sx q[2];
rz(0.26668374) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23197933) q[1];
sx q[1];
rz(-2.5464499) q[1];
sx q[1];
rz(0.91722782) q[1];
x q[2];
rz(-0.8021637) q[3];
sx q[3];
rz(-0.18076104) q[3];
sx q[3];
rz(2.4651431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58549762) q[2];
sx q[2];
rz(-1.4336339) q[2];
sx q[2];
rz(-1.1222703) q[2];
rz(1.026011) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(2.1508353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68200237) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(-0.37297747) q[0];
rz(2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(-1.3164828) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3217357) q[0];
sx q[0];
rz(-0.61367354) q[0];
sx q[0];
rz(0.86921285) q[0];
x q[1];
rz(1.3189795) q[2];
sx q[2];
rz(-1.3164815) q[2];
sx q[2];
rz(0.11676678) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.398048) q[1];
sx q[1];
rz(-0.053357031) q[1];
sx q[1];
rz(1.3392901) q[1];
rz(-pi) q[2];
rz(-1.0373711) q[3];
sx q[3];
rz(-1.08764) q[3];
sx q[3];
rz(-2.0869568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0405154) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(-0.80580795) q[2];
rz(-0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(-1.0114975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7508115) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(-0.090963013) q[0];
rz(-2.2816351) q[1];
sx q[1];
rz(-2.0188315) q[1];
sx q[1];
rz(-1.8213173) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1307756) q[0];
sx q[0];
rz(-0.59033075) q[0];
sx q[0];
rz(-0.70461313) q[0];
rz(-pi) q[1];
rz(-1.7485511) q[2];
sx q[2];
rz(-2.1338935) q[2];
sx q[2];
rz(0.69001889) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7399866) q[1];
sx q[1];
rz(-1.1005797) q[1];
sx q[1];
rz(-2.2120038) q[1];
rz(-pi) q[2];
rz(-0.91464197) q[3];
sx q[3];
rz(-1.3387036) q[3];
sx q[3];
rz(2.7616012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0423353) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(-0.60097224) q[2];
rz(-2.6565334) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(-1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27959529) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(-2.5860508) q[0];
rz(0.034596054) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(1.3909891) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23404183) q[0];
sx q[0];
rz(-0.64767917) q[0];
sx q[0];
rz(-2.5729022) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5552942) q[2];
sx q[2];
rz(-2.7099897) q[2];
sx q[2];
rz(1.7238215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6185551) q[1];
sx q[1];
rz(-2.4589775) q[1];
sx q[1];
rz(-1.7726266) q[1];
rz(1.6452541) q[3];
sx q[3];
rz(-2.3450608) q[3];
sx q[3];
rz(-1.2021241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3283219) q[2];
sx q[2];
rz(-0.44181028) q[2];
sx q[2];
rz(-0.39548809) q[2];
rz(1.8528806) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46273461) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(-1.6495552) q[0];
rz(2.18816) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(1.4377726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8673082) q[0];
sx q[0];
rz(-1.058393) q[0];
sx q[0];
rz(-0.209765) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67280025) q[2];
sx q[2];
rz(-0.2163419) q[2];
sx q[2];
rz(1.0018444) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2807323) q[1];
sx q[1];
rz(-0.8287462) q[1];
sx q[1];
rz(2.4651205) q[1];
rz(-2.0182761) q[3];
sx q[3];
rz(-1.5527225) q[3];
sx q[3];
rz(1.3656838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4961204) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(2.9628741) q[2];
rz(0.86137613) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9361967) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(1.0937011) q[0];
rz(0.73668346) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(-1.0587943) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5537542) q[0];
sx q[0];
rz(-2.5231579) q[0];
sx q[0];
rz(1.0879602) q[0];
rz(-pi) q[1];
rz(2.043622) q[2];
sx q[2];
rz(-1.9377922) q[2];
sx q[2];
rz(-2.8851913) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.79566369) q[1];
sx q[1];
rz(-1.1946031) q[1];
sx q[1];
rz(0.13161195) q[1];
rz(-pi) q[2];
rz(-1.9731673) q[3];
sx q[3];
rz(-2.0109004) q[3];
sx q[3];
rz(-0.518706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(-1.4808562) q[2];
rz(2.7311834) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(-0.15795344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720298) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(2.8503382) q[0];
rz(2.5323396) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(1.7094918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0390022) q[0];
sx q[0];
rz(-1.0800835) q[0];
sx q[0];
rz(-1.9309994) q[0];
rz(-1.3347689) q[2];
sx q[2];
rz(-1.126295) q[2];
sx q[2];
rz(2.1399463) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6546302) q[1];
sx q[1];
rz(-1.527206) q[1];
sx q[1];
rz(-0.55200465) q[1];
x q[2];
rz(1.1389144) q[3];
sx q[3];
rz(-1.4915885) q[3];
sx q[3];
rz(-1.6457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46135205) q[2];
sx q[2];
rz(-2.0246918) q[2];
sx q[2];
rz(-2.0001901) q[2];
rz(1.6067778) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(-2.6721568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948982) q[0];
sx q[0];
rz(-1.6225157) q[0];
sx q[0];
rz(-1.7058104) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(-0.73451191) q[2];
sx q[2];
rz(-0.30167689) q[2];
sx q[2];
rz(-2.6405356) q[2];
rz(-0.41714824) q[3];
sx q[3];
rz(-0.45776164) q[3];
sx q[3];
rz(-0.2840857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
