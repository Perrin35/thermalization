OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.86201) q[0];
sx q[0];
rz(-0.53505889) q[0];
sx q[0];
rz(0.23393272) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(-1.7276126) q[1];
sx q[1];
rz(1.6230621) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46759638) q[0];
sx q[0];
rz(-2.5270871) q[0];
sx q[0];
rz(-0.55802457) q[0];
rz(-pi) q[1];
rz(0.21295453) q[2];
sx q[2];
rz(-1.9969654) q[2];
sx q[2];
rz(-2.1338685) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3668574) q[1];
sx q[1];
rz(-2.1090057) q[1];
sx q[1];
rz(0.66267207) q[1];
rz(-pi) q[2];
x q[2];
rz(2.150334) q[3];
sx q[3];
rz(-2.391444) q[3];
sx q[3];
rz(-1.0446435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5412377) q[2];
sx q[2];
rz(-1.4404094) q[2];
sx q[2];
rz(-0.8055299) q[2];
rz(2.7496998) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(-2.384757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5517752) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(2.7031194) q[0];
rz(-1.27502) q[1];
sx q[1];
rz(-1.4507111) q[1];
sx q[1];
rz(-3.0335887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94919449) q[0];
sx q[0];
rz(-0.11365232) q[0];
sx q[0];
rz(1.7920919) q[0];
rz(2.503452) q[2];
sx q[2];
rz(-0.9305939) q[2];
sx q[2];
rz(-1.3253044) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0846363) q[1];
sx q[1];
rz(-2.8743187) q[1];
sx q[1];
rz(1.2483424) q[1];
rz(2.0612848) q[3];
sx q[3];
rz(-1.401859) q[3];
sx q[3];
rz(2.6712457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7763623) q[2];
sx q[2];
rz(-2.6680816) q[2];
sx q[2];
rz(3.0253809) q[2];
rz(-2.727437) q[3];
sx q[3];
rz(-2.0678346) q[3];
sx q[3];
rz(-0.80348408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47495833) q[0];
sx q[0];
rz(-1.3131498) q[0];
sx q[0];
rz(-2.3057002) q[0];
rz(1.5180786) q[1];
sx q[1];
rz(-1.514879) q[1];
sx q[1];
rz(1.9035043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9758751) q[0];
sx q[0];
rz(-1.200186) q[0];
sx q[0];
rz(2.1134645) q[0];
rz(-pi) q[1];
x q[1];
rz(2.798271) q[2];
sx q[2];
rz(-1.3385217) q[2];
sx q[2];
rz(0.49988036) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1819396) q[1];
sx q[1];
rz(-1.6283855) q[1];
sx q[1];
rz(-1.4532386) q[1];
x q[2];
rz(-0.28210552) q[3];
sx q[3];
rz(-1.8477401) q[3];
sx q[3];
rz(1.1462584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8572924) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(0.71883744) q[2];
rz(-2.5607064) q[3];
sx q[3];
rz(-1.418117) q[3];
sx q[3];
rz(-2.412793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.523664) q[0];
sx q[0];
rz(-1.5476462) q[0];
sx q[0];
rz(1.1267598) q[0];
rz(1.2976546) q[1];
sx q[1];
rz(-1.490386) q[1];
sx q[1];
rz(1.0430956) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8422416) q[0];
sx q[0];
rz(-1.8411921) q[0];
sx q[0];
rz(0.031456703) q[0];
rz(1.1642493) q[2];
sx q[2];
rz(-2.5669328) q[2];
sx q[2];
rz(-0.80917796) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.265643) q[1];
sx q[1];
rz(-1.416872) q[1];
sx q[1];
rz(-1.2624082) q[1];
rz(-pi) q[2];
rz(3.0303589) q[3];
sx q[3];
rz(-2.4832442) q[3];
sx q[3];
rz(-1.8452132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8635233) q[2];
sx q[2];
rz(-1.3845283) q[2];
sx q[2];
rz(2.1045904) q[2];
rz(2.1841124) q[3];
sx q[3];
rz(-1.0629531) q[3];
sx q[3];
rz(-2.3069042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1400414) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(3.0539883) q[0];
rz(-2.7032848) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(1.3137438) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8439595) q[0];
sx q[0];
rz(-0.54415138) q[0];
sx q[0];
rz(-1.1060957) q[0];
rz(-pi) q[1];
rz(1.9319973) q[2];
sx q[2];
rz(-3.0362077) q[2];
sx q[2];
rz(-3.1230645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.524595) q[1];
sx q[1];
rz(-1.4819078) q[1];
sx q[1];
rz(1.7434381) q[1];
rz(-0.99788061) q[3];
sx q[3];
rz(-2.8393203) q[3];
sx q[3];
rz(0.59461601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4801415) q[2];
sx q[2];
rz(-1.9037312) q[2];
sx q[2];
rz(2.5782149) q[2];
rz(-1.9208469) q[3];
sx q[3];
rz(-1.3910339) q[3];
sx q[3];
rz(-2.5105072) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0193943) q[0];
sx q[0];
rz(-2.4833184) q[0];
sx q[0];
rz(3.1296545) q[0];
rz(-2.7519233) q[1];
sx q[1];
rz(-2.4395112) q[1];
sx q[1];
rz(-0.24519244) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.747904) q[0];
sx q[0];
rz(-0.98391279) q[0];
sx q[0];
rz(-0.52397195) q[0];
x q[1];
rz(2.8831364) q[2];
sx q[2];
rz(-1.3535796) q[2];
sx q[2];
rz(2.1803792) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6840382) q[1];
sx q[1];
rz(-2.0054617) q[1];
sx q[1];
rz(-0.60575374) q[1];
x q[2];
rz(2.526029) q[3];
sx q[3];
rz(-2.6608753) q[3];
sx q[3];
rz(-1.2445039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7616854) q[2];
sx q[2];
rz(-1.8319538) q[2];
sx q[2];
rz(1.7725819) q[2];
rz(-0.53064972) q[3];
sx q[3];
rz(-2.4574418) q[3];
sx q[3];
rz(-1.0859038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1928007) q[0];
sx q[0];
rz(-1.3230319) q[0];
sx q[0];
rz(-0.2970933) q[0];
rz(-1.6311215) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(2.7105601) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6490899) q[0];
sx q[0];
rz(-1.8592872) q[0];
sx q[0];
rz(1.3516851) q[0];
x q[1];
rz(0.23081339) q[2];
sx q[2];
rz(-1.5195519) q[2];
sx q[2];
rz(1.9153999) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7707211) q[1];
sx q[1];
rz(-1.8951192) q[1];
sx q[1];
rz(2.1065708) q[1];
x q[2];
rz(2.0991574) q[3];
sx q[3];
rz(-1.8260806) q[3];
sx q[3];
rz(-1.135889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83583528) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(2.8744899) q[2];
rz(-3.109572) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(3.035868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0778462) q[0];
sx q[0];
rz(-1.2910605) q[0];
sx q[0];
rz(2.5015976) q[0];
rz(0.85482875) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(1.7291501) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4025637) q[0];
sx q[0];
rz(-1.1629346) q[0];
sx q[0];
rz(-0.5824851) q[0];
rz(-1.6663297) q[2];
sx q[2];
rz(-0.49411202) q[2];
sx q[2];
rz(1.6803368) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.58929491) q[1];
sx q[1];
rz(-1.6550242) q[1];
sx q[1];
rz(-0.34119795) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9962003) q[3];
sx q[3];
rz(-1.2266876) q[3];
sx q[3];
rz(0.95208012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5414446) q[2];
sx q[2];
rz(-1.7057799) q[2];
sx q[2];
rz(2.7195462) q[2];
rz(2.6680434) q[3];
sx q[3];
rz(-2.519042) q[3];
sx q[3];
rz(-2.9464909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7710829) q[0];
sx q[0];
rz(-0.93733731) q[0];
sx q[0];
rz(-1.7899845) q[0];
rz(-0.31556684) q[1];
sx q[1];
rz(-1.4548929) q[1];
sx q[1];
rz(1.1531166) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095799965) q[0];
sx q[0];
rz(-1.6216767) q[0];
sx q[0];
rz(-1.6232383) q[0];
rz(-pi) q[1];
rz(-2.6062327) q[2];
sx q[2];
rz(-1.6125814) q[2];
sx q[2];
rz(-0.9524065) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1156731) q[1];
sx q[1];
rz(-0.5957091) q[1];
sx q[1];
rz(0.22819789) q[1];
rz(-2.7615776) q[3];
sx q[3];
rz(-1.1855864) q[3];
sx q[3];
rz(0.64982254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6844668) q[2];
sx q[2];
rz(-1.212333) q[2];
sx q[2];
rz(-0.37377629) q[2];
rz(1.2260381) q[3];
sx q[3];
rz(-0.91807476) q[3];
sx q[3];
rz(-1.3396243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1821197) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(1.3207588) q[0];
rz(-0.93718115) q[1];
sx q[1];
rz(-1.8056185) q[1];
sx q[1];
rz(-0.46930596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.824911) q[0];
sx q[0];
rz(-1.4535722) q[0];
sx q[0];
rz(0.076923142) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9515418) q[2];
sx q[2];
rz(-2.1441438) q[2];
sx q[2];
rz(-2.9721476) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9633858) q[1];
sx q[1];
rz(-1.4904516) q[1];
sx q[1];
rz(-2.9909913) q[1];
rz(-2.8402249) q[3];
sx q[3];
rz(-1.3346938) q[3];
sx q[3];
rz(2.2599041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70539537) q[2];
sx q[2];
rz(-0.86279482) q[2];
sx q[2];
rz(1.2104642) q[2];
rz(-0.56810275) q[3];
sx q[3];
rz(-1.3848687) q[3];
sx q[3];
rz(0.97584045) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9217011) q[0];
sx q[0];
rz(-2.2911063) q[0];
sx q[0];
rz(-1.6794857) q[0];
rz(-0.89314356) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(0.22994269) q[2];
sx q[2];
rz(-1.8529006) q[2];
sx q[2];
rz(0.012250031) q[2];
rz(0.94598573) q[3];
sx q[3];
rz(-2.47249) q[3];
sx q[3];
rz(-2.5853018) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
