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
rz(0.089224815) q[0];
sx q[0];
rz(-0.74639809) q[0];
sx q[0];
rz(-2.449692) q[0];
rz(0.46299419) q[1];
sx q[1];
rz(-1.6166592) q[1];
sx q[1];
rz(-1.0301627) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44987401) q[0];
sx q[0];
rz(-2.1263459) q[0];
sx q[0];
rz(1.8434974) q[0];
rz(-pi) q[1];
rz(-1.1724418) q[2];
sx q[2];
rz(-1.7400558) q[2];
sx q[2];
rz(-1.5072849) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6840238) q[1];
sx q[1];
rz(-1.6872703) q[1];
sx q[1];
rz(-1.2435786) q[1];
rz(-pi) q[2];
rz(2.452677) q[3];
sx q[3];
rz(-1.6647881) q[3];
sx q[3];
rz(1.1979562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87972155) q[2];
sx q[2];
rz(-0.77360952) q[2];
sx q[2];
rz(-2.775178) q[2];
rz(1.0877747) q[3];
sx q[3];
rz(-2.352114) q[3];
sx q[3];
rz(-2.4521949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3159897) q[0];
sx q[0];
rz(-0.082254224) q[0];
sx q[0];
rz(0.91973037) q[0];
rz(2.3986744) q[1];
sx q[1];
rz(-1.8835604) q[1];
sx q[1];
rz(-2.0222576) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9539328) q[0];
sx q[0];
rz(-1.7155678) q[0];
sx q[0];
rz(-0.06186084) q[0];
rz(0.7190312) q[2];
sx q[2];
rz(-1.708548) q[2];
sx q[2];
rz(-3.0437623) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85038919) q[1];
sx q[1];
rz(-0.77050401) q[1];
sx q[1];
rz(2.797956) q[1];
rz(0.51471424) q[3];
sx q[3];
rz(-1.052894) q[3];
sx q[3];
rz(1.6165773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4982345) q[2];
sx q[2];
rz(-1.1946119) q[2];
sx q[2];
rz(-3.0894792) q[2];
rz(2.0335967) q[3];
sx q[3];
rz(-0.85927695) q[3];
sx q[3];
rz(3.0767379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68010083) q[0];
sx q[0];
rz(-1.9100186) q[0];
sx q[0];
rz(1.8517866) q[0];
rz(0.75311226) q[1];
sx q[1];
rz(-1.4537922) q[1];
sx q[1];
rz(-2.6444816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6996993) q[0];
sx q[0];
rz(-1.9632247) q[0];
sx q[0];
rz(-2.5470069) q[0];
rz(-pi) q[1];
rz(-2.3894823) q[2];
sx q[2];
rz(-2.3083938) q[2];
sx q[2];
rz(1.1583661) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5376007) q[1];
sx q[1];
rz(-1.6087172) q[1];
sx q[1];
rz(-2.7673519) q[1];
rz(2.4889491) q[3];
sx q[3];
rz(-1.6473928) q[3];
sx q[3];
rz(0.21993876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6734753) q[2];
sx q[2];
rz(-2.868729) q[2];
sx q[2];
rz(0.72431481) q[2];
rz(-2.916548) q[3];
sx q[3];
rz(-1.8658172) q[3];
sx q[3];
rz(-0.79506522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.986079) q[0];
sx q[0];
rz(-2.2561095) q[0];
sx q[0];
rz(-2.8610863) q[0];
rz(1.4836503) q[1];
sx q[1];
rz(-1.2018485) q[1];
sx q[1];
rz(-0.4712421) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0126033) q[0];
sx q[0];
rz(-1.5480729) q[0];
sx q[0];
rz(-1.2034334) q[0];
x q[1];
rz(1.4578826) q[2];
sx q[2];
rz(-1.6364003) q[2];
sx q[2];
rz(-2.0173617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1859101) q[1];
sx q[1];
rz(-1.5985399) q[1];
sx q[1];
rz(1.8415336) q[1];
rz(-pi) q[2];
rz(2.1625278) q[3];
sx q[3];
rz(-0.85460288) q[3];
sx q[3];
rz(-1.5997052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19807854) q[2];
sx q[2];
rz(-0.55800617) q[2];
sx q[2];
rz(2.1264709) q[2];
rz(-2.4070168) q[3];
sx q[3];
rz(-1.7769122) q[3];
sx q[3];
rz(2.6034897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839142) q[0];
sx q[0];
rz(-1.9240802) q[0];
sx q[0];
rz(-0.96018803) q[0];
rz(-0.045711191) q[1];
sx q[1];
rz(-0.8784596) q[1];
sx q[1];
rz(0.033230573) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93532256) q[0];
sx q[0];
rz(-1.066813) q[0];
sx q[0];
rz(-2.5978243) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.303238) q[2];
sx q[2];
rz(-1.55624) q[2];
sx q[2];
rz(0.47786682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14151621) q[1];
sx q[1];
rz(-1.7011397) q[1];
sx q[1];
rz(-1.9059859) q[1];
rz(-pi) q[2];
rz(1.1636249) q[3];
sx q[3];
rz(-1.9800751) q[3];
sx q[3];
rz(0.57151505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9783322) q[2];
sx q[2];
rz(-2.1425118) q[2];
sx q[2];
rz(2.8283289) q[2];
rz(-0.35798171) q[3];
sx q[3];
rz(-2.101818) q[3];
sx q[3];
rz(-3.0430005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1545496) q[0];
sx q[0];
rz(-1.6313666) q[0];
sx q[0];
rz(-0.48598591) q[0];
rz(-1.9588574) q[1];
sx q[1];
rz(-2.4458838) q[1];
sx q[1];
rz(-0.33014306) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81583038) q[0];
sx q[0];
rz(-1.6005524) q[0];
sx q[0];
rz(1.8682006) q[0];
x q[1];
rz(-2.5146944) q[2];
sx q[2];
rz(-2.2316602) q[2];
sx q[2];
rz(-3.0783184) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3482326) q[1];
sx q[1];
rz(-2.4168682) q[1];
sx q[1];
rz(-0.7591922) q[1];
rz(2.7948912) q[3];
sx q[3];
rz(-1.7809521) q[3];
sx q[3];
rz(-0.29276785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5064454) q[2];
sx q[2];
rz(-2.5250285) q[2];
sx q[2];
rz(0.35663024) q[2];
rz(-0.59456524) q[3];
sx q[3];
rz(-1.4831355) q[3];
sx q[3];
rz(2.4026292) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67721382) q[0];
sx q[0];
rz(-0.76736275) q[0];
sx q[0];
rz(-0.58694029) q[0];
rz(-0.37239536) q[1];
sx q[1];
rz(-1.7814691) q[1];
sx q[1];
rz(2.8940103) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0599759) q[0];
sx q[0];
rz(-2.6756786) q[0];
sx q[0];
rz(-1.640725) q[0];
rz(1.3273986) q[2];
sx q[2];
rz(-2.1888424) q[2];
sx q[2];
rz(-2.7400561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3629439) q[1];
sx q[1];
rz(-1.2374068) q[1];
sx q[1];
rz(0.89493805) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60815292) q[3];
sx q[3];
rz(-2.7233015) q[3];
sx q[3];
rz(-2.2564512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26329142) q[2];
sx q[2];
rz(-1.2211439) q[2];
sx q[2];
rz(-0.2090052) q[2];
rz(-0.59669295) q[3];
sx q[3];
rz(-2.7946819) q[3];
sx q[3];
rz(1.6056812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40808943) q[0];
sx q[0];
rz(-0.96036378) q[0];
sx q[0];
rz(-2.7485513) q[0];
rz(0.6706925) q[1];
sx q[1];
rz(-1.9978943) q[1];
sx q[1];
rz(1.447698) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2358952) q[0];
sx q[0];
rz(-1.6950771) q[0];
sx q[0];
rz(1.0727364) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1850519) q[2];
sx q[2];
rz(-1.8280085) q[2];
sx q[2];
rz(-1.9458015) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.61673855) q[1];
sx q[1];
rz(-0.752512) q[1];
sx q[1];
rz(2.3630211) q[1];
rz(-pi) q[2];
rz(2.6781769) q[3];
sx q[3];
rz(-0.8532195) q[3];
sx q[3];
rz(-1.0427102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.712901) q[2];
sx q[2];
rz(-0.8608326) q[2];
sx q[2];
rz(-0.11670308) q[2];
rz(-2.1278925) q[3];
sx q[3];
rz(-1.2319535) q[3];
sx q[3];
rz(-1.3056171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75334615) q[0];
sx q[0];
rz(-1.4660864) q[0];
sx q[0];
rz(1.7062794) q[0];
rz(0.995579) q[1];
sx q[1];
rz(-1.3187871) q[1];
sx q[1];
rz(-0.54042712) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89622241) q[0];
sx q[0];
rz(-2.5858552) q[0];
sx q[0];
rz(-0.51387365) q[0];
x q[1];
rz(-0.86825722) q[2];
sx q[2];
rz(-0.80108713) q[2];
sx q[2];
rz(-0.67745393) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0173954) q[1];
sx q[1];
rz(-0.81995839) q[1];
sx q[1];
rz(-0.52738054) q[1];
rz(2.5164817) q[3];
sx q[3];
rz(-1.907699) q[3];
sx q[3];
rz(0.57527104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2702177) q[2];
sx q[2];
rz(-1.3048708) q[2];
sx q[2];
rz(0.77667856) q[2];
rz(-3.1089879) q[3];
sx q[3];
rz(-1.6083345) q[3];
sx q[3];
rz(-1.5672654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33083522) q[0];
sx q[0];
rz(-2.7726655) q[0];
sx q[0];
rz(2.708013) q[0];
rz(-0.66917229) q[1];
sx q[1];
rz(-1.1829665) q[1];
sx q[1];
rz(-0.42632595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0678372) q[0];
sx q[0];
rz(-1.3569419) q[0];
sx q[0];
rz(-2.9392089) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9258921) q[2];
sx q[2];
rz(-0.85596687) q[2];
sx q[2];
rz(-1.0155884) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59407104) q[1];
sx q[1];
rz(-1.921355) q[1];
sx q[1];
rz(-0.5188397) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0845523) q[3];
sx q[3];
rz(-2.1003688) q[3];
sx q[3];
rz(-2.2467006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7508037) q[2];
sx q[2];
rz(-2.8533253) q[2];
sx q[2];
rz(2.0354347) q[2];
rz(-0.11354167) q[3];
sx q[3];
rz(-0.93324408) q[3];
sx q[3];
rz(-0.043702628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9857367) q[0];
sx q[0];
rz(-0.899213) q[0];
sx q[0];
rz(-0.56094299) q[0];
rz(-2.7583495) q[1];
sx q[1];
rz(-1.1226729) q[1];
sx q[1];
rz(2.9164006) q[1];
rz(0.55632527) q[2];
sx q[2];
rz(-1.3090324) q[2];
sx q[2];
rz(-0.39568452) q[2];
rz(-0.81859891) q[3];
sx q[3];
rz(-1.9055453) q[3];
sx q[3];
rz(0.26858134) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
