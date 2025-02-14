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
rz(-2.9774732) q[0];
sx q[0];
rz(-1.0241221) q[0];
sx q[0];
rz(2.6574988) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(-1.5634544) q[1];
sx q[1];
rz(-0.054952316) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30357519) q[0];
sx q[0];
rz(-0.9009046) q[0];
sx q[0];
rz(-0.039867) q[0];
rz(-pi) q[1];
rz(-1.7968494) q[2];
sx q[2];
rz(-0.25616562) q[2];
sx q[2];
rz(2.0992172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.41356) q[1];
sx q[1];
rz(-1.7031341) q[1];
sx q[1];
rz(1.4721666) q[1];
x q[2];
rz(1.1160158) q[3];
sx q[3];
rz(-0.6729799) q[3];
sx q[3];
rz(-2.2448886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.310828) q[2];
sx q[2];
rz(-0.97042933) q[2];
sx q[2];
rz(-0.27466276) q[2];
rz(0.87485391) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(2.8680475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65983588) q[0];
sx q[0];
rz(-2.4517224) q[0];
sx q[0];
rz(-2.7428395) q[0];
rz(-2.8975471) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(2.0358548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4717446) q[0];
sx q[0];
rz(-2.6378184) q[0];
sx q[0];
rz(-2.6543447) q[0];
rz(2.2695838) q[2];
sx q[2];
rz(-2.3559743) q[2];
sx q[2];
rz(2.9532022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9404906) q[1];
sx q[1];
rz(-0.0041714287) q[1];
sx q[1];
rz(-0.68912403) q[1];
rz(-pi) q[2];
rz(1.9883184) q[3];
sx q[3];
rz(-2.8510749) q[3];
sx q[3];
rz(-2.6444496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70046052) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(2.0665118) q[2];
rz(-2.4454146) q[3];
sx q[3];
rz(-1.2735561) q[3];
sx q[3];
rz(-0.16304326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27580801) q[0];
sx q[0];
rz(-1.2514665) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(-1.7614583) q[1];
sx q[1];
rz(-0.41788995) q[1];
sx q[1];
rz(-2.5221141) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.769128) q[0];
sx q[0];
rz(-1.541168) q[0];
sx q[0];
rz(1.6086786) q[0];
rz(-pi) q[1];
rz(1.2080433) q[2];
sx q[2];
rz(-0.13158509) q[2];
sx q[2];
rz(-0.059748273) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.463355) q[1];
sx q[1];
rz(-1.2690407) q[1];
sx q[1];
rz(0.22308087) q[1];
rz(-pi) q[2];
rz(-0.9914753) q[3];
sx q[3];
rz(-2.1001171) q[3];
sx q[3];
rz(-2.4059699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2565903) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(0.090506434) q[2];
rz(-2.6835486) q[3];
sx q[3];
rz(-2.6543255) q[3];
sx q[3];
rz(-1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3312382) q[0];
sx q[0];
rz(-0.88307035) q[0];
sx q[0];
rz(0.93165398) q[0];
rz(0.16904198) q[1];
sx q[1];
rz(-2.5773498) q[1];
sx q[1];
rz(2.1021252) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19379481) q[0];
sx q[0];
rz(-0.76776869) q[0];
sx q[0];
rz(-0.34996535) q[0];
x q[1];
rz(2.9157964) q[2];
sx q[2];
rz(-1.3232627) q[2];
sx q[2];
rz(-3.053726) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.89371496) q[1];
sx q[1];
rz(-1.7778686) q[1];
sx q[1];
rz(3.0805001) q[1];
x q[2];
rz(1.3576034) q[3];
sx q[3];
rz(-1.8490095) q[3];
sx q[3];
rz(1.9813862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.52183759) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(0.49631611) q[2];
rz(0.79658341) q[3];
sx q[3];
rz(-0.28592548) q[3];
sx q[3];
rz(-2.0485785) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26688823) q[0];
sx q[0];
rz(-2.5570091) q[0];
sx q[0];
rz(2.1697178) q[0];
rz(-0.12174363) q[1];
sx q[1];
rz(-1.6106482) q[1];
sx q[1];
rz(1.8121388) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055453528) q[0];
sx q[0];
rz(-2.8290966) q[0];
sx q[0];
rz(-1.3700831) q[0];
rz(1.0779641) q[2];
sx q[2];
rz(-1.9415641) q[2];
sx q[2];
rz(-2.855122) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54534528) q[1];
sx q[1];
rz(-1.9831428) q[1];
sx q[1];
rz(-0.72280563) q[1];
x q[2];
rz(1.6869808) q[3];
sx q[3];
rz(-0.079616485) q[3];
sx q[3];
rz(-1.5594359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10151265) q[2];
sx q[2];
rz(-1.9209346) q[2];
sx q[2];
rz(2.9902048) q[2];
rz(-2.3769412) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5966655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0695892) q[0];
sx q[0];
rz(-1.66865) q[0];
sx q[0];
rz(1.0714916) q[0];
rz(-0.53120652) q[1];
sx q[1];
rz(-1.4899645) q[1];
sx q[1];
rz(0.62613097) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1086639) q[0];
sx q[0];
rz(-1.5877921) q[0];
sx q[0];
rz(-3.1368318) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7886969) q[2];
sx q[2];
rz(-2.5145686) q[2];
sx q[2];
rz(0.74909808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.31589139) q[1];
sx q[1];
rz(-2.4448312) q[1];
sx q[1];
rz(0.48721643) q[1];
rz(-pi) q[2];
rz(0.49895309) q[3];
sx q[3];
rz(-1.1756983) q[3];
sx q[3];
rz(0.56624352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7834187) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(-2.1273071) q[2];
rz(1.2554393) q[3];
sx q[3];
rz(-1.3557326) q[3];
sx q[3];
rz(2.0881418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805873) q[0];
sx q[0];
rz(-0.87562457) q[0];
sx q[0];
rz(-0.92051202) q[0];
rz(3.0275184) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(-1.1309518) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3944954) q[0];
sx q[0];
rz(-2.1685026) q[0];
sx q[0];
rz(1.5639485) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3384096) q[2];
sx q[2];
rz(-2.4121662) q[2];
sx q[2];
rz(0.27952295) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1753833) q[1];
sx q[1];
rz(-1.7031324) q[1];
sx q[1];
rz(-2.2369409) q[1];
rz(-pi) q[2];
rz(-1.6507963) q[3];
sx q[3];
rz(-1.0558075) q[3];
sx q[3];
rz(1.1752216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7941234) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(-0.40880173) q[2];
rz(0.18320228) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(2.0028152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.113753) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(-0.29079944) q[0];
rz(0.94789061) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(-1.1999757) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5644218) q[0];
sx q[0];
rz(-2.0977712) q[0];
sx q[0];
rz(2.0354664) q[0];
rz(-pi) q[1];
rz(0.84917111) q[2];
sx q[2];
rz(-2.2361922) q[2];
sx q[2];
rz(0.1553436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0713745) q[1];
sx q[1];
rz(-1.6171675) q[1];
sx q[1];
rz(1.2769615) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7539976) q[3];
sx q[3];
rz(-1.6699413) q[3];
sx q[3];
rz(2.9354036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7802508) q[2];
sx q[2];
rz(-1.4486518) q[2];
sx q[2];
rz(-1.7928436) q[2];
rz(-1.5032984) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(0.1851113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1143484) q[0];
sx q[0];
rz(-0.37186563) q[0];
sx q[0];
rz(1.0741023) q[0];
rz(2.7117924) q[1];
sx q[1];
rz(-2.0975515) q[1];
sx q[1];
rz(2.6780186) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.539042) q[0];
sx q[0];
rz(-2.9304404) q[0];
sx q[0];
rz(1.1013012) q[0];
rz(0.65166574) q[2];
sx q[2];
rz(-2.8164688) q[2];
sx q[2];
rz(2.5420957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1164828) q[1];
sx q[1];
rz(-1.9466725) q[1];
sx q[1];
rz(2.2459338) q[1];
rz(-pi) q[2];
x q[2];
rz(1.660668) q[3];
sx q[3];
rz(-2.1276132) q[3];
sx q[3];
rz(1.1695605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8672436) q[2];
sx q[2];
rz(-2.5090019) q[2];
sx q[2];
rz(-0.80844936) q[2];
rz(0.48458734) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(-0.43420473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3633858) q[0];
sx q[0];
rz(-0.45553842) q[0];
sx q[0];
rz(-2.0090012) q[0];
rz(0.038381902) q[1];
sx q[1];
rz(-1.3382341) q[1];
sx q[1];
rz(2.8448232) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1365741) q[0];
sx q[0];
rz(-1.8474192) q[0];
sx q[0];
rz(1.0199976) q[0];
rz(0.44906868) q[2];
sx q[2];
rz(-1.9757604) q[2];
sx q[2];
rz(-1.8835889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3546875) q[1];
sx q[1];
rz(-2.9070435) q[1];
sx q[1];
rz(-0.3292747) q[1];
x q[2];
rz(-0.88480437) q[3];
sx q[3];
rz(-1.2594885) q[3];
sx q[3];
rz(-1.5112359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.177629) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(-2.5753944) q[2];
rz(-2.0171793) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2531256) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(-0.26168564) q[1];
sx q[1];
rz(-1.9410004) q[1];
sx q[1];
rz(0.95071361) q[1];
rz(3.059584) q[2];
sx q[2];
rz(-1.5951984) q[2];
sx q[2];
rz(-2.3328608) q[2];
rz(0.86033173) q[3];
sx q[3];
rz(-2.0900149) q[3];
sx q[3];
rz(2.7505977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
