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
rz(1.053831) q[0];
sx q[0];
rz(3.6917917) q[0];
sx q[0];
rz(10.79296) q[0];
rz(-0.47222459) q[1];
sx q[1];
rz(-3.0250186) q[1];
sx q[1];
rz(-0.20224686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41101) q[0];
sx q[0];
rz(-0.29313335) q[0];
sx q[0];
rz(-0.75384753) q[0];
rz(-pi) q[1];
rz(-2.6361385) q[2];
sx q[2];
rz(-1.1040282) q[2];
sx q[2];
rz(2.9708178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7649012) q[1];
sx q[1];
rz(-0.67923949) q[1];
sx q[1];
rz(-2.0308073) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1038001) q[3];
sx q[3];
rz(-1.1049912) q[3];
sx q[3];
rz(-0.58096262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54174417) q[2];
sx q[2];
rz(-0.84459633) q[2];
sx q[2];
rz(-2.2005626) q[2];
rz(2.8372676) q[3];
sx q[3];
rz(-0.86186886) q[3];
sx q[3];
rz(0.26722515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32352725) q[0];
sx q[0];
rz(-0.37779385) q[0];
sx q[0];
rz(-0.7793119) q[0];
rz(-1.7794973) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(2.0077226) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2529566) q[0];
sx q[0];
rz(-1.1373113) q[0];
sx q[0];
rz(0.10826464) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.629584) q[2];
sx q[2];
rz(-1.4530621) q[2];
sx q[2];
rz(1.1716154) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.459455) q[1];
sx q[1];
rz(-1.3054779) q[1];
sx q[1];
rz(3.0299761) q[1];
x q[2];
rz(1.5725582) q[3];
sx q[3];
rz(-2.632318) q[3];
sx q[3];
rz(0.60520303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4506932) q[2];
sx q[2];
rz(-2.1096114) q[2];
sx q[2];
rz(-2.9449985) q[2];
rz(0.96902668) q[3];
sx q[3];
rz(-1.6220379) q[3];
sx q[3];
rz(1.2381747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75152385) q[0];
sx q[0];
rz(-2.5876434) q[0];
sx q[0];
rz(2.4098136) q[0];
rz(-2.7732908) q[1];
sx q[1];
rz(-0.75804561) q[1];
sx q[1];
rz(-0.36645737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5598504) q[0];
sx q[0];
rz(-1.3936211) q[0];
sx q[0];
rz(-0.26818256) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24213893) q[2];
sx q[2];
rz(-1.79091) q[2];
sx q[2];
rz(1.5907703) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.428047) q[1];
sx q[1];
rz(-2.4427892) q[1];
sx q[1];
rz(-2.2627463) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59753527) q[3];
sx q[3];
rz(-0.96216494) q[3];
sx q[3];
rz(-0.017692117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1825819) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(0.89243531) q[2];
rz(2.6871032) q[3];
sx q[3];
rz(-2.317704) q[3];
sx q[3];
rz(-0.23196001) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7682122) q[0];
sx q[0];
rz(-0.1425655) q[0];
sx q[0];
rz(-0.40661231) q[0];
rz(-2.3136102) q[1];
sx q[1];
rz(-1.7384638) q[1];
sx q[1];
rz(2.3060395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260561) q[0];
sx q[0];
rz(-2.9565507) q[0];
sx q[0];
rz(-1.3547784) q[0];
x q[1];
rz(-1.2242975) q[2];
sx q[2];
rz(-1.76148) q[2];
sx q[2];
rz(2.4257223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9753127) q[1];
sx q[1];
rz(-1.1381835) q[1];
sx q[1];
rz(1.7325425) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9560247) q[3];
sx q[3];
rz(-1.480417) q[3];
sx q[3];
rz(-0.65672311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39997175) q[2];
sx q[2];
rz(-0.30210945) q[2];
sx q[2];
rz(-0.92140222) q[2];
rz(-1.7737927) q[3];
sx q[3];
rz(-2.1801345) q[3];
sx q[3];
rz(-3.0000946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88857404) q[0];
sx q[0];
rz(-1.0971917) q[0];
sx q[0];
rz(-0.15884037) q[0];
rz(1.5171492) q[1];
sx q[1];
rz(-0.30888638) q[1];
sx q[1];
rz(-1.812017) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3811262) q[0];
sx q[0];
rz(-1.5601349) q[0];
sx q[0];
rz(-3.1233379) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9407523) q[2];
sx q[2];
rz(-2.5479377) q[2];
sx q[2];
rz(-1.6933189) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.161474) q[1];
sx q[1];
rz(-1.1921645) q[1];
sx q[1];
rz(2.9254854) q[1];
rz(-pi) q[2];
rz(-1.8908126) q[3];
sx q[3];
rz(-0.59139267) q[3];
sx q[3];
rz(2.0689912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5627731) q[2];
sx q[2];
rz(-0.69391888) q[2];
sx q[2];
rz(2.5895183) q[2];
rz(-0.79365802) q[3];
sx q[3];
rz(-2.5439883) q[3];
sx q[3];
rz(2.1551267) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8296705) q[0];
sx q[0];
rz(-0.86450082) q[0];
sx q[0];
rz(-0.70575869) q[0];
rz(-3.0041079) q[1];
sx q[1];
rz(-0.64422137) q[1];
sx q[1];
rz(-0.71298832) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5891083) q[0];
sx q[0];
rz(-0.45408598) q[0];
sx q[0];
rz(0.54291351) q[0];
rz(-pi) q[1];
rz(2.627264) q[2];
sx q[2];
rz(-1.9514958) q[2];
sx q[2];
rz(-1.2027539) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2299774) q[1];
sx q[1];
rz(-1.4739081) q[1];
sx q[1];
rz(-1.1712892) q[1];
rz(1.2778132) q[3];
sx q[3];
rz(-1.6821386) q[3];
sx q[3];
rz(-2.7822943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9299499) q[2];
sx q[2];
rz(-0.86923081) q[2];
sx q[2];
rz(0.28406528) q[2];
rz(0.12187135) q[3];
sx q[3];
rz(-2.8594696) q[3];
sx q[3];
rz(-1.4819283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.876494) q[0];
sx q[0];
rz(-2.5229186) q[0];
sx q[0];
rz(0.70736831) q[0];
rz(2.7012198) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(-1.3622989) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.52584) q[0];
sx q[0];
rz(-2.7670015) q[0];
sx q[0];
rz(2.2994141) q[0];
rz(-pi) q[1];
rz(2.0411891) q[2];
sx q[2];
rz(-0.78436034) q[2];
sx q[2];
rz(-2.228594) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55864298) q[1];
sx q[1];
rz(-2.042965) q[1];
sx q[1];
rz(0.98319816) q[1];
rz(-0.23162095) q[3];
sx q[3];
rz(-0.4670139) q[3];
sx q[3];
rz(-1.2488332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1661487) q[2];
sx q[2];
rz(-0.41836172) q[2];
sx q[2];
rz(0.56900209) q[2];
rz(2.1042018) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(2.67498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.55050945) q[0];
sx q[0];
rz(-0.35174462) q[0];
sx q[0];
rz(1.2048703) q[0];
rz(-2.6288746) q[1];
sx q[1];
rz(-2.3725489) q[1];
sx q[1];
rz(0.18096322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0807225) q[0];
sx q[0];
rz(-1.1747169) q[0];
sx q[0];
rz(-1.164308) q[0];
rz(2.7920807) q[2];
sx q[2];
rz(-1.4415359) q[2];
sx q[2];
rz(1.845128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1026336) q[1];
sx q[1];
rz(-2.5397083) q[1];
sx q[1];
rz(-1.9485056) q[1];
x q[2];
rz(3.0309903) q[3];
sx q[3];
rz(-1.6440653) q[3];
sx q[3];
rz(1.4721118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5391431) q[2];
sx q[2];
rz(-2.4931543) q[2];
sx q[2];
rz(-3.0446206) q[2];
rz(0.55475956) q[3];
sx q[3];
rz(-1.3223038) q[3];
sx q[3];
rz(0.35840148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48034126) q[0];
sx q[0];
rz(-2.3619409) q[0];
sx q[0];
rz(-0.10777792) q[0];
rz(-0.55094552) q[1];
sx q[1];
rz(-0.44083732) q[1];
sx q[1];
rz(-1.5239747) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0227528) q[0];
sx q[0];
rz(-1.700239) q[0];
sx q[0];
rz(1.3368827) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0650159) q[2];
sx q[2];
rz(-0.98978251) q[2];
sx q[2];
rz(2.6083824) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2366166) q[1];
sx q[1];
rz(-1.1091241) q[1];
sx q[1];
rz(3.0932337) q[1];
rz(-1.518256) q[3];
sx q[3];
rz(-2.0618304) q[3];
sx q[3];
rz(2.4201916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9109351) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(-2.0111734) q[2];
rz(2.9039827) q[3];
sx q[3];
rz(-2.7311324) q[3];
sx q[3];
rz(-0.75806481) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1199353) q[0];
sx q[0];
rz(-0.1788685) q[0];
sx q[0];
rz(-0.94104952) q[0];
rz(-0.36048105) q[1];
sx q[1];
rz(-1.4070114) q[1];
sx q[1];
rz(-1.9182659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4775765) q[0];
sx q[0];
rz(-1.7682791) q[0];
sx q[0];
rz(2.5504179) q[0];
x q[1];
rz(0.98127301) q[2];
sx q[2];
rz(-2.2726577) q[2];
sx q[2];
rz(-1.6035994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66050038) q[1];
sx q[1];
rz(-1.6410488) q[1];
sx q[1];
rz(0.41892799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9452865) q[3];
sx q[3];
rz(-2.098791) q[3];
sx q[3];
rz(-1.4083901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2337522) q[2];
sx q[2];
rz(-1.3774104) q[2];
sx q[2];
rz(-0.09093786) q[2];
rz(2.9371373) q[3];
sx q[3];
rz(-0.8046059) q[3];
sx q[3];
rz(-1.0819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7627056) q[0];
sx q[0];
rz(-1.611041) q[0];
sx q[0];
rz(1.8086717) q[0];
rz(-1.7145722) q[1];
sx q[1];
rz(-0.49976977) q[1];
sx q[1];
rz(1.5907092) q[1];
rz(-2.0841523) q[2];
sx q[2];
rz(-0.79887894) q[2];
sx q[2];
rz(0.32499921) q[2];
rz(-2.9279999) q[3];
sx q[3];
rz(-1.5145258) q[3];
sx q[3];
rz(0.35888844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
