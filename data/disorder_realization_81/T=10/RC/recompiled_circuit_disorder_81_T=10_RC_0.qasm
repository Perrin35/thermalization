OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(-1.7629495) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(4.6255914) q[1];
sx q[1];
rz(8.9738823) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.742813) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(0.16422693) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0300006) q[2];
sx q[2];
rz(-1.0942232) q[2];
sx q[2];
rz(0.0026207844) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3528459) q[1];
sx q[1];
rz(-0.55410085) q[1];
sx q[1];
rz(-0.43177859) q[1];
x q[2];
rz(1.6879184) q[3];
sx q[3];
rz(-2.4544567) q[3];
sx q[3];
rz(0.62108921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(0.84428865) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-2.8785008) q[0];
rz(0.94353765) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.9553604) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3620136) q[0];
sx q[0];
rz(-1.5148666) q[0];
sx q[0];
rz(1.5520142) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8430014) q[2];
sx q[2];
rz(-0.208374) q[2];
sx q[2];
rz(1.6076455) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7344208) q[1];
sx q[1];
rz(-2.2312575) q[1];
sx q[1];
rz(0.15863005) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4466189) q[3];
sx q[3];
rz(-1.4222099) q[3];
sx q[3];
rz(-0.23651628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0120323) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(-2.7705079) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(2.8306567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.7611258) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(-2.3348715) q[0];
rz(0.21356788) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7263111) q[0];
sx q[0];
rz(-1.6831241) q[0];
sx q[0];
rz(2.2583654) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21913146) q[2];
sx q[2];
rz(-1.0872772) q[2];
sx q[2];
rz(-2.6205274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.080938235) q[1];
sx q[1];
rz(-1.2186236) q[1];
sx q[1];
rz(2.0209795) q[1];
rz(1.3385653) q[3];
sx q[3];
rz(-2.665822) q[3];
sx q[3];
rz(0.25482086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.31072581) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(2.2107928) q[2];
rz(0.15549774) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(0.91127515) q[0];
rz(-0.43831929) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(1.8211676) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60162773) q[0];
sx q[0];
rz(-2.7318582) q[0];
sx q[0];
rz(-1.5967303) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7772061) q[2];
sx q[2];
rz(-2.0136535) q[2];
sx q[2];
rz(-2.4368311) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.98852324) q[1];
sx q[1];
rz(-1.3077902) q[1];
sx q[1];
rz(-0.70446976) q[1];
x q[2];
rz(0.96418013) q[3];
sx q[3];
rz(-2.2607431) q[3];
sx q[3];
rz(2.6026158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(-2.7992115) q[2];
rz(-2.9648182) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(-2.2763021) q[0];
rz(1.226549) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(-1.3006166) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1039935) q[0];
sx q[0];
rz(-0.35498699) q[0];
sx q[0];
rz(-0.07261891) q[0];
x q[1];
rz(-0.36655764) q[2];
sx q[2];
rz(-1.6379116) q[2];
sx q[2];
rz(-1.3973438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.865766) q[1];
sx q[1];
rz(-1.8991125) q[1];
sx q[1];
rz(2.5715716) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6404068) q[3];
sx q[3];
rz(-0.29705829) q[3];
sx q[3];
rz(-2.3230769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(2.664393) q[2];
rz(0.19208433) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(-2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79648298) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(3.1298424) q[0];
rz(0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(1.5531497) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71205157) q[0];
sx q[0];
rz(-2.6459604) q[0];
sx q[0];
rz(0.80722157) q[0];
rz(-1.6278218) q[2];
sx q[2];
rz(-2.7751623) q[2];
sx q[2];
rz(2.1232405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.8223871) q[1];
sx q[1];
rz(-2.2777595) q[1];
sx q[1];
rz(-1.6824526) q[1];
rz(2.2736069) q[3];
sx q[3];
rz(-1.6276974) q[3];
sx q[3];
rz(1.9753319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6340296) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(-1.2825512) q[2];
rz(1.7717308) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(-1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(-2.989785) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(2.1645434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86977406) q[0];
sx q[0];
rz(-0.5973814) q[0];
sx q[0];
rz(-0.13418829) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1408429) q[2];
sx q[2];
rz(-0.13204083) q[2];
sx q[2];
rz(-0.11242871) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9338785) q[1];
sx q[1];
rz(-1.2320476) q[1];
sx q[1];
rz(-3.0599041) q[1];
x q[2];
rz(-0.3427152) q[3];
sx q[3];
rz(-1.8748771) q[3];
sx q[3];
rz(-2.8976687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3892422) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(2.1288669) q[2];
rz(-1.1879454) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(-2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1241207) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(2.0195122) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(1.2493856) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.218924) q[0];
sx q[0];
rz(-1.7814753) q[0];
sx q[0];
rz(1.6181437) q[0];
rz(0.84455372) q[2];
sx q[2];
rz(-2.484998) q[2];
sx q[2];
rz(-3.055228) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3932712) q[1];
sx q[1];
rz(-1.5652579) q[1];
sx q[1];
rz(-3.1194411) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2035355) q[3];
sx q[3];
rz(-1.8055827) q[3];
sx q[3];
rz(0.42890047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32968783) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(1.9630986) q[2];
rz(-1.684749) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6417398) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(1.8485803) q[0];
rz(1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(2.5440149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9874728) q[0];
sx q[0];
rz(-0.14053908) q[0];
sx q[0];
rz(-1.8494291) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3915646) q[2];
sx q[2];
rz(-2.4366597) q[2];
sx q[2];
rz(0.36545576) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9012007) q[1];
sx q[1];
rz(-1.8933834) q[1];
sx q[1];
rz(2.8278973) q[1];
x q[2];
rz(-1.4937917) q[3];
sx q[3];
rz(-2.1553851) q[3];
sx q[3];
rz(-1.0494174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9188345) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(-1.9082327) q[2];
rz(-0.90138609) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(-1.6433158) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047886588) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(3.1179324) q[0];
rz(-0.95611447) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-2.4694209) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4155054) q[0];
sx q[0];
rz(-3.130075) q[0];
sx q[0];
rz(-1.0938209) q[0];
x q[1];
rz(1.5800843) q[2];
sx q[2];
rz(-2.0822968) q[2];
sx q[2];
rz(1.7363422) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2939261) q[1];
sx q[1];
rz(-0.77984174) q[1];
sx q[1];
rz(-2.6508209) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8176114) q[3];
sx q[3];
rz(-2.6880662) q[3];
sx q[3];
rz(-1.6739664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6293634) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(-2.771634) q[2];
rz(1.5036748) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(-1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.5621915) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(-2.4178986) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(-1.9696708) q[2];
sx q[2];
rz(-1.4234067) q[2];
sx q[2];
rz(1.1248551) q[2];
rz(-1.1643812) q[3];
sx q[3];
rz(-1.8462528) q[3];
sx q[3];
rz(0.13312199) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
