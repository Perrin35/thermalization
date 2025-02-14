OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.42216766) q[0];
sx q[0];
rz(-2.4165805) q[0];
sx q[0];
rz(-0.83591914) q[0];
rz(-2.3140276) q[1];
sx q[1];
rz(-0.27048549) q[1];
sx q[1];
rz(-2.7343813) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93210685) q[0];
sx q[0];
rz(-1.3582241) q[0];
sx q[0];
rz(1.2680156) q[0];
x q[1];
rz(1.9303481) q[2];
sx q[2];
rz(-1.3463194) q[2];
sx q[2];
rz(1.6752072) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29556698) q[1];
sx q[1];
rz(-2.7252619) q[1];
sx q[1];
rz(-2.8041072) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9933957) q[3];
sx q[3];
rz(-0.90809408) q[3];
sx q[3];
rz(-2.4277594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3621346) q[2];
sx q[2];
rz(-2.802765) q[2];
sx q[2];
rz(-2.1753878) q[2];
rz(-2.2400014) q[3];
sx q[3];
rz(-0.85351557) q[3];
sx q[3];
rz(-1.0454267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40039429) q[0];
sx q[0];
rz(-0.61120954) q[0];
sx q[0];
rz(0.10261593) q[0];
rz(-1.1688894) q[1];
sx q[1];
rz(-0.97016197) q[1];
sx q[1];
rz(0.51365596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4054905) q[0];
sx q[0];
rz(-1.0879192) q[0];
sx q[0];
rz(3.0391284) q[0];
rz(-pi) q[1];
rz(-1.7134708) q[2];
sx q[2];
rz(-2.1032984) q[2];
sx q[2];
rz(-0.14217686) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8508965) q[1];
sx q[1];
rz(-2.4833931) q[1];
sx q[1];
rz(1.0653516) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30269514) q[3];
sx q[3];
rz(-1.9871357) q[3];
sx q[3];
rz(0.2442322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1932842) q[2];
sx q[2];
rz(-0.19364348) q[2];
sx q[2];
rz(0.2085169) q[2];
rz(0.68566132) q[3];
sx q[3];
rz(-2.0833569) q[3];
sx q[3];
rz(2.9764777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.044428069) q[0];
sx q[0];
rz(-1.8655638) q[0];
sx q[0];
rz(1.9097419) q[0];
rz(0.34230289) q[1];
sx q[1];
rz(-1.1016223) q[1];
sx q[1];
rz(-1.5922348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78944713) q[0];
sx q[0];
rz(-1.2593126) q[0];
sx q[0];
rz(-1.3793263) q[0];
rz(-pi) q[1];
rz(0.64022417) q[2];
sx q[2];
rz(-2.9415543) q[2];
sx q[2];
rz(-2.6286516) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3433088) q[1];
sx q[1];
rz(-0.73241703) q[1];
sx q[1];
rz(1.364752) q[1];
rz(2.1560043) q[3];
sx q[3];
rz(-1.1051559) q[3];
sx q[3];
rz(1.5805494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0937097) q[2];
sx q[2];
rz(-1.1992531) q[2];
sx q[2];
rz(-1.017978) q[2];
rz(-1.7300946) q[3];
sx q[3];
rz(-0.74677765) q[3];
sx q[3];
rz(1.5311034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8133076) q[0];
sx q[0];
rz(-1.3628553) q[0];
sx q[0];
rz(0.12411975) q[0];
rz(0.92503754) q[1];
sx q[1];
rz(-1.3003474) q[1];
sx q[1];
rz(-0.73840028) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9772241) q[0];
sx q[0];
rz(-1.508655) q[0];
sx q[0];
rz(-2.3232099) q[0];
rz(2.5321711) q[2];
sx q[2];
rz(-1.8928693) q[2];
sx q[2];
rz(0.24064482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99216813) q[1];
sx q[1];
rz(-2.1271884) q[1];
sx q[1];
rz(1.8144649) q[1];
rz(1.4703512) q[3];
sx q[3];
rz(-2.1444706) q[3];
sx q[3];
rz(2.3015882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0783405) q[2];
sx q[2];
rz(-1.1880778) q[2];
sx q[2];
rz(2.038302) q[2];
rz(-2.8830146) q[3];
sx q[3];
rz(-1.0753814) q[3];
sx q[3];
rz(-2.8598089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2461808) q[0];
sx q[0];
rz(-0.87765944) q[0];
sx q[0];
rz(1.3377162) q[0];
rz(0.61996639) q[1];
sx q[1];
rz(-1.4620616) q[1];
sx q[1];
rz(-1.5043129) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97457394) q[0];
sx q[0];
rz(-1.7593547) q[0];
sx q[0];
rz(1.2861286) q[0];
rz(-pi) q[1];
rz(-0.70233924) q[2];
sx q[2];
rz(-1.5389256) q[2];
sx q[2];
rz(-1.7447217) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.729709) q[1];
sx q[1];
rz(-2.4296654) q[1];
sx q[1];
rz(0.61213778) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4535459) q[3];
sx q[3];
rz(-0.73620376) q[3];
sx q[3];
rz(0.29090009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6164246) q[2];
sx q[2];
rz(-1.1526356) q[2];
sx q[2];
rz(-2.0174513) q[2];
rz(-0.21640402) q[3];
sx q[3];
rz(-2.3085322) q[3];
sx q[3];
rz(2.0719349) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6935317) q[0];
sx q[0];
rz(-0.8140642) q[0];
sx q[0];
rz(2.6779209) q[0];
rz(-3.0060153) q[1];
sx q[1];
rz(-2.1453073) q[1];
sx q[1];
rz(-2.7574976) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4736872) q[0];
sx q[0];
rz(-1.1039005) q[0];
sx q[0];
rz(-1.6993194) q[0];
rz(-pi) q[1];
rz(-0.76445396) q[2];
sx q[2];
rz(-1.4234418) q[2];
sx q[2];
rz(2.5558228) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2781394) q[1];
sx q[1];
rz(-1.7378686) q[1];
sx q[1];
rz(-1.9394622) q[1];
x q[2];
rz(1.6946402) q[3];
sx q[3];
rz(-2.5899124) q[3];
sx q[3];
rz(1.1744896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7163081) q[2];
sx q[2];
rz(-2.7649438) q[2];
sx q[2];
rz(2.6728805) q[2];
rz(-2.5675755) q[3];
sx q[3];
rz(-1.6711957) q[3];
sx q[3];
rz(2.4911657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18844093) q[0];
sx q[0];
rz(-1.3583536) q[0];
sx q[0];
rz(2.4515732) q[0];
rz(0.72235876) q[1];
sx q[1];
rz(-1.1411618) q[1];
sx q[1];
rz(-2.3648327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9463766) q[0];
sx q[0];
rz(-1.4958943) q[0];
sx q[0];
rz(-2.9086065) q[0];
x q[1];
rz(2.5525307) q[2];
sx q[2];
rz(-1.8839594) q[2];
sx q[2];
rz(2.2472266) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0563363) q[1];
sx q[1];
rz(-1.83167) q[1];
sx q[1];
rz(0.44020997) q[1];
rz(-1.5473821) q[3];
sx q[3];
rz(-0.6802313) q[3];
sx q[3];
rz(1.7643203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40778273) q[2];
sx q[2];
rz(-1.5747384) q[2];
sx q[2];
rz(3.0877647) q[2];
rz(0.6221866) q[3];
sx q[3];
rz(-0.48431188) q[3];
sx q[3];
rz(-0.079306451) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4902041) q[0];
sx q[0];
rz(-1.1505928) q[0];
sx q[0];
rz(-1.401061) q[0];
rz(2.1531064) q[1];
sx q[1];
rz(-1.2683615) q[1];
sx q[1];
rz(-2.7338457) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9765104) q[0];
sx q[0];
rz(-0.4645949) q[0];
sx q[0];
rz(2.6376758) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2456535) q[2];
sx q[2];
rz(-1.4714339) q[2];
sx q[2];
rz(-0.69420136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0235879) q[1];
sx q[1];
rz(-2.0366001) q[1];
sx q[1];
rz(-0.76264571) q[1];
rz(-pi) q[2];
rz(1.4325822) q[3];
sx q[3];
rz(-1.301389) q[3];
sx q[3];
rz(2.9215653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.72655788) q[2];
sx q[2];
rz(-1.1316391) q[2];
sx q[2];
rz(2.4416907) q[2];
rz(-2.0578201) q[3];
sx q[3];
rz(-2.2743069) q[3];
sx q[3];
rz(2.9054902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.20339762) q[0];
sx q[0];
rz(-1.6827826) q[0];
sx q[0];
rz(-0.38988018) q[0];
rz(-1.1979206) q[1];
sx q[1];
rz(-0.094466297) q[1];
sx q[1];
rz(0.87497154) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.724204) q[0];
sx q[0];
rz(-1.5483138) q[0];
sx q[0];
rz(2.1940986) q[0];
rz(1.0880555) q[2];
sx q[2];
rz(-2.3953279) q[2];
sx q[2];
rz(-1.7604699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7851049) q[1];
sx q[1];
rz(-0.56863368) q[1];
sx q[1];
rz(-2.4604753) q[1];
rz(-2.7194831) q[3];
sx q[3];
rz(-1.0448467) q[3];
sx q[3];
rz(1.6933481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8209057) q[2];
sx q[2];
rz(-0.69183886) q[2];
sx q[2];
rz(2.4917277) q[2];
rz(-2.0248905) q[3];
sx q[3];
rz(-1.2973123) q[3];
sx q[3];
rz(0.19702774) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7277302) q[0];
sx q[0];
rz(-3.08857) q[0];
sx q[0];
rz(-2.9470288) q[0];
rz(0.77231705) q[1];
sx q[1];
rz(-1.1338502) q[1];
sx q[1];
rz(0.17759855) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5935737) q[0];
sx q[0];
rz(-1.5691891) q[0];
sx q[0];
rz(3.1409825) q[0];
rz(-2.3760711) q[2];
sx q[2];
rz(-2.3280848) q[2];
sx q[2];
rz(-0.080527079) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2706621) q[1];
sx q[1];
rz(-2.1829603) q[1];
sx q[1];
rz(0.98302676) q[1];
rz(-pi) q[2];
rz(2.1924413) q[3];
sx q[3];
rz(-2.2165934) q[3];
sx q[3];
rz(-0.27816712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3557768) q[2];
sx q[2];
rz(-2.5444701) q[2];
sx q[2];
rz(2.7115278) q[2];
rz(2.0307342) q[3];
sx q[3];
rz(-1.6836124) q[3];
sx q[3];
rz(-2.7250169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72660245) q[0];
sx q[0];
rz(-1.5885329) q[0];
sx q[0];
rz(2.9581099) q[0];
rz(1.968374) q[1];
sx q[1];
rz(-1.219974) q[1];
sx q[1];
rz(-3.0164607) q[1];
rz(-2.0424615) q[2];
sx q[2];
rz(-2.5644414) q[2];
sx q[2];
rz(1.9664398) q[2];
rz(0.54456768) q[3];
sx q[3];
rz(-0.93251317) q[3];
sx q[3];
rz(2.9208356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
