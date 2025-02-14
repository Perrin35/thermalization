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
rz(2.501261) q[0];
sx q[0];
rz(-2.2202272) q[0];
sx q[0];
rz(-1.6093572) q[0];
rz(-2.9345203) q[1];
sx q[1];
rz(-1.1595885) q[1];
sx q[1];
rz(-1.6699189) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5014316) q[0];
sx q[0];
rz(-2.6894173) q[0];
sx q[0];
rz(0.78340952) q[0];
rz(-pi) q[1];
rz(-2.4406836) q[2];
sx q[2];
rz(-0.94628382) q[2];
sx q[2];
rz(0.33901843) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.629327) q[1];
sx q[1];
rz(-1.6696603) q[1];
sx q[1];
rz(-1.3310286) q[1];
rz(-pi) q[2];
rz(-0.59572409) q[3];
sx q[3];
rz(-1.7215773) q[3];
sx q[3];
rz(0.25546701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49429911) q[2];
sx q[2];
rz(-0.99227253) q[2];
sx q[2];
rz(0.54599071) q[2];
rz(-0.93849385) q[3];
sx q[3];
rz(-2.6280845) q[3];
sx q[3];
rz(-2.688496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-0.42354217) q[1];
sx q[1];
rz(2.0372527) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1416407) q[0];
sx q[0];
rz(-1.4237836) q[0];
sx q[0];
rz(-2.5070058) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0510089) q[2];
sx q[2];
rz(-2.5414645) q[2];
sx q[2];
rz(-2.5243197) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2939289) q[1];
sx q[1];
rz(-2.1184455) q[1];
sx q[1];
rz(-0.54137648) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0542778) q[3];
sx q[3];
rz(-1.5305291) q[3];
sx q[3];
rz(-0.95424679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7426593) q[2];
sx q[2];
rz(-1.089596) q[2];
sx q[2];
rz(-2.147414) q[2];
rz(-1.582877) q[3];
sx q[3];
rz(-2.4623058) q[3];
sx q[3];
rz(2.5905632) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1347374) q[0];
sx q[0];
rz(-2.0515433) q[0];
sx q[0];
rz(0.63049522) q[0];
rz(0.14671239) q[1];
sx q[1];
rz(-1.881733) q[1];
sx q[1];
rz(0.86348081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16374198) q[0];
sx q[0];
rz(-1.8713179) q[0];
sx q[0];
rz(-0.88877166) q[0];
rz(-pi) q[1];
rz(2.9296257) q[2];
sx q[2];
rz(-1.8463928) q[2];
sx q[2];
rz(-0.26108643) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.024115) q[1];
sx q[1];
rz(-1.4653135) q[1];
sx q[1];
rz(2.0085667) q[1];
x q[2];
rz(0.15210882) q[3];
sx q[3];
rz(-2.8242691) q[3];
sx q[3];
rz(-3.139024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33519393) q[2];
sx q[2];
rz(-1.4525745) q[2];
sx q[2];
rz(-0.88649583) q[2];
rz(-2.7116306) q[3];
sx q[3];
rz(-2.6468247) q[3];
sx q[3];
rz(-2.2494242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48915136) q[0];
sx q[0];
rz(-2.8471071) q[0];
sx q[0];
rz(-2.6926706) q[0];
rz(-2.4849675) q[1];
sx q[1];
rz(-1.4337599) q[1];
sx q[1];
rz(0.33033672) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2170439) q[0];
sx q[0];
rz(-1.346949) q[0];
sx q[0];
rz(-1.2428817) q[0];
x q[1];
rz(-1.0895715) q[2];
sx q[2];
rz(-0.46483609) q[2];
sx q[2];
rz(-2.6770153) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8395665) q[1];
sx q[1];
rz(-2.6077787) q[1];
sx q[1];
rz(0.97179504) q[1];
rz(-3.1325601) q[3];
sx q[3];
rz(-1.6918283) q[3];
sx q[3];
rz(1.5435404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.088048) q[2];
sx q[2];
rz(-1.5379173) q[2];
sx q[2];
rz(-0.21928445) q[2];
rz(-2.3413279) q[3];
sx q[3];
rz(-0.95992464) q[3];
sx q[3];
rz(0.76103359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6680172) q[0];
sx q[0];
rz(-1.1290843) q[0];
sx q[0];
rz(-1.6582723) q[0];
rz(2.5738916) q[1];
sx q[1];
rz(-1.9300108) q[1];
sx q[1];
rz(1.6426881) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53500861) q[0];
sx q[0];
rz(-2.5902915) q[0];
sx q[0];
rz(-0.28321858) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11168555) q[2];
sx q[2];
rz(-3.0431261) q[2];
sx q[2];
rz(2.6278327) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.65693457) q[1];
sx q[1];
rz(-0.96862786) q[1];
sx q[1];
rz(1.8116211) q[1];
rz(-pi) q[2];
rz(-1.3793702) q[3];
sx q[3];
rz(-1.4512806) q[3];
sx q[3];
rz(-1.9724979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6296926) q[2];
sx q[2];
rz(-0.7581768) q[2];
sx q[2];
rz(-2.8153815) q[2];
rz(0.1263667) q[3];
sx q[3];
rz(-0.56505239) q[3];
sx q[3];
rz(-0.27157426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7997953) q[0];
sx q[0];
rz(-1.3038776) q[0];
sx q[0];
rz(-2.6190992) q[0];
rz(-2.3937468) q[1];
sx q[1];
rz(-1.6035085) q[1];
sx q[1];
rz(-1.1572256) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25230468) q[0];
sx q[0];
rz(-1.8631019) q[0];
sx q[0];
rz(-2.0415123) q[0];
rz(-2.2713178) q[2];
sx q[2];
rz(-0.57862568) q[2];
sx q[2];
rz(-3.0491997) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5392235) q[1];
sx q[1];
rz(-1.7134949) q[1];
sx q[1];
rz(-0.093009503) q[1];
rz(-pi) q[2];
rz(0.22208235) q[3];
sx q[3];
rz(-2.8068868) q[3];
sx q[3];
rz(-0.10192733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.130827) q[2];
sx q[2];
rz(-2.6289434) q[2];
sx q[2];
rz(-1.1831247) q[2];
rz(0.042923953) q[3];
sx q[3];
rz(-1.5019417) q[3];
sx q[3];
rz(-0.24421282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44508988) q[0];
sx q[0];
rz(-0.88560167) q[0];
sx q[0];
rz(2.4687299) q[0];
rz(0.54975763) q[1];
sx q[1];
rz(-1.8472698) q[1];
sx q[1];
rz(-1.4901644) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8880418) q[0];
sx q[0];
rz(-1.7490938) q[0];
sx q[0];
rz(0.56714296) q[0];
rz(0.76810683) q[2];
sx q[2];
rz(-0.66181493) q[2];
sx q[2];
rz(2.2827374) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0615599) q[1];
sx q[1];
rz(-1.6668238) q[1];
sx q[1];
rz(-2.705234) q[1];
rz(-2.1792322) q[3];
sx q[3];
rz(-3.0828028) q[3];
sx q[3];
rz(2.1421681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0979536) q[2];
sx q[2];
rz(-2.592228) q[2];
sx q[2];
rz(-1.9926386) q[2];
rz(2.7502381) q[3];
sx q[3];
rz(-1.2201744) q[3];
sx q[3];
rz(-2.2933188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5963762) q[0];
sx q[0];
rz(-1.9356118) q[0];
sx q[0];
rz(-2.4429831) q[0];
rz(0.85083234) q[1];
sx q[1];
rz(-2.5406676) q[1];
sx q[1];
rz(0.94815475) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2744926) q[0];
sx q[0];
rz(-2.4153408) q[0];
sx q[0];
rz(-0.12402835) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7265375) q[2];
sx q[2];
rz(-2.808185) q[2];
sx q[2];
rz(-1.6108212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25094704) q[1];
sx q[1];
rz(-2.1275824) q[1];
sx q[1];
rz(-2.6614038) q[1];
x q[2];
rz(0.73581477) q[3];
sx q[3];
rz(-2.1664005) q[3];
sx q[3];
rz(-0.99887139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4725388) q[2];
sx q[2];
rz(-1.6927745) q[2];
sx q[2];
rz(1.4581468) q[2];
rz(2.0504047) q[3];
sx q[3];
rz(-0.89717054) q[3];
sx q[3];
rz(2.9567961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75337306) q[0];
sx q[0];
rz(-0.18167697) q[0];
sx q[0];
rz(-1.8411807) q[0];
rz(-0.45937195) q[1];
sx q[1];
rz(-1.6875024) q[1];
sx q[1];
rz(2.557911) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65410173) q[0];
sx q[0];
rz(-1.0922474) q[0];
sx q[0];
rz(-2.3951981) q[0];
rz(-pi) q[1];
rz(0.038554474) q[2];
sx q[2];
rz(-1.5351982) q[2];
sx q[2];
rz(2.1447139) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.7129659) q[1];
sx q[1];
rz(-1.9216864) q[1];
sx q[1];
rz(-2.116973) q[1];
rz(-pi) q[2];
rz(-2.8902938) q[3];
sx q[3];
rz(-1.150395) q[3];
sx q[3];
rz(0.29743089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5251081) q[2];
sx q[2];
rz(-2.8664092) q[2];
sx q[2];
rz(2.6194561) q[2];
rz(0.76256847) q[3];
sx q[3];
rz(-1.6430166) q[3];
sx q[3];
rz(-0.34408072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6793215) q[0];
sx q[0];
rz(-1.9336047) q[0];
sx q[0];
rz(2.5479877) q[0];
rz(2.4484334) q[1];
sx q[1];
rz(-0.79265541) q[1];
sx q[1];
rz(2.3902182) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8639785) q[0];
sx q[0];
rz(-1.1728012) q[0];
sx q[0];
rz(-1.2897183) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2840653) q[2];
sx q[2];
rz(-0.81999841) q[2];
sx q[2];
rz(-1.4992901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6804816) q[1];
sx q[1];
rz(-2.2348677) q[1];
sx q[1];
rz(1.6517488) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7165224) q[3];
sx q[3];
rz(-0.18501013) q[3];
sx q[3];
rz(-1.4671221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0005325) q[2];
sx q[2];
rz(-2.4324721) q[2];
sx q[2];
rz(-3.1179023) q[2];
rz(1.116811) q[3];
sx q[3];
rz(-1.4477372) q[3];
sx q[3];
rz(2.007133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5524207) q[0];
sx q[0];
rz(-1.2180653) q[0];
sx q[0];
rz(1.1267452) q[0];
rz(-1.8213656) q[1];
sx q[1];
rz(-1.8658493) q[1];
sx q[1];
rz(1.0125926) q[1];
rz(0.54426469) q[2];
sx q[2];
rz(-1.7901292) q[2];
sx q[2];
rz(-2.3397056) q[2];
rz(-2.9145225) q[3];
sx q[3];
rz(-1.3454701) q[3];
sx q[3];
rz(-3.0440192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
