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
rz(0.47705874) q[0];
sx q[0];
rz(-0.30204371) q[0];
sx q[0];
rz(1.8855236) q[0];
rz(-2.7657901) q[1];
sx q[1];
rz(-1.8355651) q[1];
sx q[1];
rz(-2.0130872) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80917847) q[0];
sx q[0];
rz(-2.1398394) q[0];
sx q[0];
rz(-0.33359417) q[0];
rz(-pi) q[1];
rz(-1.8774421) q[2];
sx q[2];
rz(-2.1696343) q[2];
sx q[2];
rz(1.6055589) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0839786) q[1];
sx q[1];
rz(-0.65257665) q[1];
sx q[1];
rz(-2.233391) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5079191) q[3];
sx q[3];
rz(-1.9320935) q[3];
sx q[3];
rz(-0.47508815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6065373) q[2];
sx q[2];
rz(-2.8474319) q[2];
sx q[2];
rz(-1.199022) q[2];
rz(-2.8626275) q[3];
sx q[3];
rz(-1.6821034) q[3];
sx q[3];
rz(0.038486686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26175427) q[0];
sx q[0];
rz(-2.6848875) q[0];
sx q[0];
rz(2.7754011) q[0];
rz(-1.2884033) q[1];
sx q[1];
rz(-0.90694702) q[1];
sx q[1];
rz(0.2624951) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27124559) q[0];
sx q[0];
rz(-2.7923005) q[0];
sx q[0];
rz(2.5906934) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8196202) q[2];
sx q[2];
rz(-1.5142528) q[2];
sx q[2];
rz(2.4849934) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9001635) q[1];
sx q[1];
rz(-1.2821322) q[1];
sx q[1];
rz(-2.7735387) q[1];
rz(-pi) q[2];
rz(-1.8257687) q[3];
sx q[3];
rz(-1.5951968) q[3];
sx q[3];
rz(-2.5687664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.75258201) q[2];
sx q[2];
rz(-1.2685308) q[2];
sx q[2];
rz(0.1788204) q[2];
rz(-3.1255152) q[3];
sx q[3];
rz(-2.6060846) q[3];
sx q[3];
rz(2.2107562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43612424) q[0];
sx q[0];
rz(-3.0146764) q[0];
sx q[0];
rz(-0.57417589) q[0];
rz(2.9899959) q[1];
sx q[1];
rz(-0.44436374) q[1];
sx q[1];
rz(-1.9134329) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41802619) q[0];
sx q[0];
rz(-2.0834588) q[0];
sx q[0];
rz(2.961757) q[0];
rz(-pi) q[1];
rz(-1.4692233) q[2];
sx q[2];
rz(-0.59528643) q[2];
sx q[2];
rz(-0.02073076) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0361402) q[1];
sx q[1];
rz(-2.5220118) q[1];
sx q[1];
rz(-1.416227) q[1];
x q[2];
rz(-2.0545928) q[3];
sx q[3];
rz(-1.3143149) q[3];
sx q[3];
rz(-2.2352122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0000618) q[2];
sx q[2];
rz(-2.5059293) q[2];
sx q[2];
rz(-2.2504811) q[2];
rz(-0.38854232) q[3];
sx q[3];
rz(-1.5107692) q[3];
sx q[3];
rz(-0.33963206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6252839) q[0];
sx q[0];
rz(-3.0966274) q[0];
sx q[0];
rz(2.6594824) q[0];
rz(-0.66857839) q[1];
sx q[1];
rz(-1.6079638) q[1];
sx q[1];
rz(-0.63749981) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26698819) q[0];
sx q[0];
rz(-2.1529196) q[0];
sx q[0];
rz(-0.28018392) q[0];
rz(-0.11580616) q[2];
sx q[2];
rz(-2.2409332) q[2];
sx q[2];
rz(-2.0320803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1177534) q[1];
sx q[1];
rz(-1.6477403) q[1];
sx q[1];
rz(-1.5368098) q[1];
rz(-pi) q[2];
rz(-0.25084514) q[3];
sx q[3];
rz(-1.224784) q[3];
sx q[3];
rz(3.0437247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.94347) q[2];
sx q[2];
rz(-2.1169457) q[2];
sx q[2];
rz(3.0812954) q[2];
rz(-0.30334011) q[3];
sx q[3];
rz(-3.0637432) q[3];
sx q[3];
rz(-2.8764603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(2.9918793) q[0];
sx q[0];
rz(-0.35357058) q[0];
sx q[0];
rz(0.41325945) q[0];
rz(2.7464271) q[1];
sx q[1];
rz(-1.7839849) q[1];
sx q[1];
rz(-1.0023592) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7089139) q[0];
sx q[0];
rz(-2.1122547) q[0];
sx q[0];
rz(-1.0546636) q[0];
x q[1];
rz(-1.7384612) q[2];
sx q[2];
rz(-0.21279003) q[2];
sx q[2];
rz(-2.8820702) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.614568) q[1];
sx q[1];
rz(-2.2348512) q[1];
sx q[1];
rz(-2.8675957) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2484364) q[3];
sx q[3];
rz(-1.3731706) q[3];
sx q[3];
rz(1.2569497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5137382) q[2];
sx q[2];
rz(-2.0444874) q[2];
sx q[2];
rz(1.6045137) q[2];
rz(1.9419935) q[3];
sx q[3];
rz(-2.0204085) q[3];
sx q[3];
rz(-2.3698923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62726778) q[0];
sx q[0];
rz(-1.5031313) q[0];
sx q[0];
rz(-2.7728873) q[0];
rz(0.74327028) q[1];
sx q[1];
rz(-2.5786046) q[1];
sx q[1];
rz(2.7875913) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6021332) q[0];
sx q[0];
rz(-1.6230427) q[0];
sx q[0];
rz(0.9297855) q[0];
x q[1];
rz(-2.6145489) q[2];
sx q[2];
rz(-1.2979002) q[2];
sx q[2];
rz(2.5986586) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56997715) q[1];
sx q[1];
rz(-1.0034402) q[1];
sx q[1];
rz(-2.9302271) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5052222) q[3];
sx q[3];
rz(-1.2546347) q[3];
sx q[3];
rz(-0.13383257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.02504286) q[2];
sx q[2];
rz(-1.9321238) q[2];
sx q[2];
rz(0.0054736007) q[2];
rz(-0.51760751) q[3];
sx q[3];
rz(-2.7766683) q[3];
sx q[3];
rz(-3.1143809) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7356877) q[0];
sx q[0];
rz(-0.55957496) q[0];
sx q[0];
rz(2.9285808) q[0];
rz(-1.6464015) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(-0.84943736) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0809459) q[0];
sx q[0];
rz(-0.95260006) q[0];
sx q[0];
rz(-1.5175876) q[0];
x q[1];
rz(-3.1108625) q[2];
sx q[2];
rz(-1.5275914) q[2];
sx q[2];
rz(-1.2655366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.60472) q[1];
sx q[1];
rz(-1.1013796) q[1];
sx q[1];
rz(-1.0632656) q[1];
x q[2];
rz(1.0424741) q[3];
sx q[3];
rz(-1.436621) q[3];
sx q[3];
rz(3.1229916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9269632) q[2];
sx q[2];
rz(-1.9081076) q[2];
sx q[2];
rz(-2.5978973) q[2];
rz(2.8804273) q[3];
sx q[3];
rz(-1.5615014) q[3];
sx q[3];
rz(-2.9847667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1866622) q[0];
sx q[0];
rz(-2.08827) q[0];
sx q[0];
rz(2.9539811) q[0];
rz(-2.018351) q[1];
sx q[1];
rz(-0.76485991) q[1];
sx q[1];
rz(-0.16256464) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99062067) q[0];
sx q[0];
rz(-0.91998581) q[0];
sx q[0];
rz(-1.9810173) q[0];
x q[1];
rz(1.8674054) q[2];
sx q[2];
rz(-2.5979418) q[2];
sx q[2];
rz(-2.6988876) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9342207) q[1];
sx q[1];
rz(-1.6160864) q[1];
sx q[1];
rz(0.95691935) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1670893) q[3];
sx q[3];
rz(-1.664229) q[3];
sx q[3];
rz(3.1395885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8673914) q[2];
sx q[2];
rz(-0.31495366) q[2];
sx q[2];
rz(-0.20142041) q[2];
rz(2.6650186) q[3];
sx q[3];
rz(-1.7628935) q[3];
sx q[3];
rz(0.99982888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.9607214) q[0];
sx q[0];
rz(-2.9044386) q[0];
sx q[0];
rz(-0.57688212) q[0];
rz(-2.8374529) q[1];
sx q[1];
rz(-1.1241333) q[1];
sx q[1];
rz(-2.0374128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7542412) q[0];
sx q[0];
rz(-1.1925444) q[0];
sx q[0];
rz(3.0463112) q[0];
x q[1];
rz(-0.56892345) q[2];
sx q[2];
rz(-1.4683791) q[2];
sx q[2];
rz(-1.2889287) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4889646) q[1];
sx q[1];
rz(-1.9581989) q[1];
sx q[1];
rz(2.1968958) q[1];
x q[2];
rz(-1.5315834) q[3];
sx q[3];
rz(-2.750142) q[3];
sx q[3];
rz(1.9658364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.96678174) q[2];
sx q[2];
rz(-0.2468144) q[2];
sx q[2];
rz(0.68320572) q[2];
rz(-1.5934058) q[3];
sx q[3];
rz(-2.4631409) q[3];
sx q[3];
rz(-0.21827503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.9504647) q[0];
sx q[0];
rz(-2.2948965) q[0];
sx q[0];
rz(0.51093131) q[0];
rz(1.1448942) q[1];
sx q[1];
rz(-1.1868008) q[1];
sx q[1];
rz(-1.6330382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3630229) q[0];
sx q[0];
rz(-1.8196882) q[0];
sx q[0];
rz(0.56015862) q[0];
rz(-pi) q[1];
rz(-2.0548923) q[2];
sx q[2];
rz(-1.8082613) q[2];
sx q[2];
rz(-2.9959842) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3737693) q[1];
sx q[1];
rz(-2.0160185) q[1];
sx q[1];
rz(-1.7675464) q[1];
rz(-pi) q[2];
rz(-2.8415751) q[3];
sx q[3];
rz(-1.3277866) q[3];
sx q[3];
rz(-0.21081616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.42023429) q[2];
sx q[2];
rz(-2.5098462) q[2];
sx q[2];
rz(1.7101804) q[2];
rz(-0.015627705) q[3];
sx q[3];
rz(-1.2651919) q[3];
sx q[3];
rz(-2.6354852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90841993) q[0];
sx q[0];
rz(-1.7411727) q[0];
sx q[0];
rz(2.2199051) q[0];
rz(-1.635101) q[1];
sx q[1];
rz(-1.2163305) q[1];
sx q[1];
rz(-0.46179927) q[1];
rz(2.2057381) q[2];
sx q[2];
rz(-2.6781185) q[2];
sx q[2];
rz(0.92739633) q[2];
rz(-1.9962068) q[3];
sx q[3];
rz(-0.085193188) q[3];
sx q[3];
rz(-1.9151158) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
