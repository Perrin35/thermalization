OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9659757) q[0];
sx q[0];
rz(-1.5953925) q[0];
sx q[0];
rz(1.6311128) q[0];
rz(-2.053396) q[1];
sx q[1];
rz(4.4463867) q[1];
sx q[1];
rz(10.211853) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2505456) q[0];
sx q[0];
rz(-0.68792397) q[0];
sx q[0];
rz(2.964818) q[0];
x q[1];
rz(-1.6252615) q[2];
sx q[2];
rz(-1.174841) q[2];
sx q[2];
rz(-1.2690282) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8382227) q[1];
sx q[1];
rz(-1.7821687) q[1];
sx q[1];
rz(1.6081393) q[1];
rz(2.3863411) q[3];
sx q[3];
rz(-0.77338615) q[3];
sx q[3];
rz(2.7735047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8481019) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(2.9254986) q[2];
rz(-0.50348336) q[3];
sx q[3];
rz(-1.2516021) q[3];
sx q[3];
rz(-0.066472806) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9965432) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(0.37522069) q[0];
rz(-2.5253865) q[1];
sx q[1];
rz(-1.0169949) q[1];
sx q[1];
rz(0.18888758) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8804042) q[0];
sx q[0];
rz(-1.1587901) q[0];
sx q[0];
rz(2.3212577) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99449956) q[2];
sx q[2];
rz(-0.81953555) q[2];
sx q[2];
rz(-2.799965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5423373) q[1];
sx q[1];
rz(-2.2399678) q[1];
sx q[1];
rz(-2.8821844) q[1];
rz(-1.9697207) q[3];
sx q[3];
rz(-2.0634983) q[3];
sx q[3];
rz(-1.0176095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6836493) q[2];
sx q[2];
rz(-1.8234437) q[2];
sx q[2];
rz(-1.1434435) q[2];
rz(-2.5055366) q[3];
sx q[3];
rz(-0.053074107) q[3];
sx q[3];
rz(-1.599357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77883333) q[0];
sx q[0];
rz(-0.23985671) q[0];
sx q[0];
rz(-1.1676316) q[0];
rz(3.0150343) q[1];
sx q[1];
rz(-1.5153706) q[1];
sx q[1];
rz(-0.57356858) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9605947) q[0];
sx q[0];
rz(-0.74686909) q[0];
sx q[0];
rz(-1.7043714) q[0];
rz(-pi) q[1];
rz(-0.36678183) q[2];
sx q[2];
rz(-1.3202808) q[2];
sx q[2];
rz(-0.17630746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4497609) q[1];
sx q[1];
rz(-1.1402601) q[1];
sx q[1];
rz(-0.052816347) q[1];
rz(-1.1352121) q[3];
sx q[3];
rz(-0.95802973) q[3];
sx q[3];
rz(1.719914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.48587376) q[2];
sx q[2];
rz(-1.3081009) q[2];
sx q[2];
rz(-1.5029079) q[2];
rz(-1.2949299) q[3];
sx q[3];
rz(-0.27532268) q[3];
sx q[3];
rz(2.3143229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5591705) q[0];
sx q[0];
rz(-2.9952413) q[0];
sx q[0];
rz(-2.225112) q[0];
rz(-0.16042635) q[1];
sx q[1];
rz(-1.8564936) q[1];
sx q[1];
rz(-0.697335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.682919) q[0];
sx q[0];
rz(-1.3240755) q[0];
sx q[0];
rz(-1.9015067) q[0];
rz(2.7576912) q[2];
sx q[2];
rz(-2.4426201) q[2];
sx q[2];
rz(2.0052103) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5494255) q[1];
sx q[1];
rz(-2.2695955) q[1];
sx q[1];
rz(-1.1392085) q[1];
x q[2];
rz(-3.0352108) q[3];
sx q[3];
rz(-0.13088317) q[3];
sx q[3];
rz(-1.037055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9852898) q[2];
sx q[2];
rz(-2.6142945) q[2];
sx q[2];
rz(-2.018833) q[2];
rz(-1.5491693) q[3];
sx q[3];
rz(-1.6808108) q[3];
sx q[3];
rz(0.39337015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3093981) q[0];
sx q[0];
rz(-0.10646146) q[0];
sx q[0];
rz(-1.1291809) q[0];
rz(-0.87942266) q[1];
sx q[1];
rz(-1.8303454) q[1];
sx q[1];
rz(2.5130491) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.382269) q[0];
sx q[0];
rz(-0.17062561) q[0];
sx q[0];
rz(-2.3294446) q[0];
rz(-pi) q[1];
rz(-2.4833335) q[2];
sx q[2];
rz(-1.1193917) q[2];
sx q[2];
rz(-2.4413204) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4820833) q[1];
sx q[1];
rz(-2.1225559) q[1];
sx q[1];
rz(-3.0943046) q[1];
rz(-pi) q[2];
rz(-1.4416297) q[3];
sx q[3];
rz(-2.3979275) q[3];
sx q[3];
rz(-0.63275933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92155543) q[2];
sx q[2];
rz(-0.76801378) q[2];
sx q[2];
rz(-1.7680232) q[2];
rz(-2.8288362) q[3];
sx q[3];
rz(-2.0650568) q[3];
sx q[3];
rz(-3.055174) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81109989) q[0];
sx q[0];
rz(-0.70527768) q[0];
sx q[0];
rz(-0.16264859) q[0];
rz(-1.9074408) q[1];
sx q[1];
rz(-1.1470497) q[1];
sx q[1];
rz(-0.92076921) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61495691) q[0];
sx q[0];
rz(-1.4702073) q[0];
sx q[0];
rz(-1.7593764) q[0];
rz(-2.2170794) q[2];
sx q[2];
rz(-2.0343896) q[2];
sx q[2];
rz(1.4626423) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7850254) q[1];
sx q[1];
rz(-1.5058702) q[1];
sx q[1];
rz(0.75218997) q[1];
rz(2.9716039) q[3];
sx q[3];
rz(-2.3521479) q[3];
sx q[3];
rz(1.1657749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84772253) q[2];
sx q[2];
rz(-0.72124481) q[2];
sx q[2];
rz(3.0976307) q[2];
rz(1.4219159) q[3];
sx q[3];
rz(-0.43216643) q[3];
sx q[3];
rz(0.034984263) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1304355) q[0];
sx q[0];
rz(-2.3373117) q[0];
sx q[0];
rz(0.093611896) q[0];
rz(2.1162107) q[1];
sx q[1];
rz(-1.5954115) q[1];
sx q[1];
rz(-2.3536918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2219395) q[0];
sx q[0];
rz(-2.1062231) q[0];
sx q[0];
rz(-0.8485283) q[0];
rz(2.5727243) q[2];
sx q[2];
rz(-1.7143357) q[2];
sx q[2];
rz(1.6678866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3026784) q[1];
sx q[1];
rz(-0.60015772) q[1];
sx q[1];
rz(0.072235302) q[1];
x q[2];
rz(-0.1206081) q[3];
sx q[3];
rz(-1.4292307) q[3];
sx q[3];
rz(2.0742576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0660144) q[2];
sx q[2];
rz(-2.5983073) q[2];
sx q[2];
rz(-2.1755966) q[2];
rz(-2.994359) q[3];
sx q[3];
rz(-1.4598264) q[3];
sx q[3];
rz(0.60432965) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886993) q[0];
sx q[0];
rz(-1.7558782) q[0];
sx q[0];
rz(-1.8889486) q[0];
rz(-3.0839651) q[1];
sx q[1];
rz(-1.7689972) q[1];
sx q[1];
rz(-0.39841121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1712043) q[0];
sx q[0];
rz(-2.1792767) q[0];
sx q[0];
rz(-1.358991) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4008825) q[2];
sx q[2];
rz(-1.0771829) q[2];
sx q[2];
rz(-0.60442558) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9908473) q[1];
sx q[1];
rz(-1.7707033) q[1];
sx q[1];
rz(0.0041109494) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54315217) q[3];
sx q[3];
rz(-1.6622322) q[3];
sx q[3];
rz(-2.9989797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5742699) q[2];
sx q[2];
rz(-1.8696573) q[2];
sx q[2];
rz(2.2197913) q[2];
rz(0.70720339) q[3];
sx q[3];
rz(-1.0450398) q[3];
sx q[3];
rz(2.8309477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1929753) q[0];
sx q[0];
rz(-3.0204168) q[0];
sx q[0];
rz(2.2344672) q[0];
rz(2.6616197) q[1];
sx q[1];
rz(-2.0860806) q[1];
sx q[1];
rz(2.5352246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7146695) q[0];
sx q[0];
rz(-1.6712609) q[0];
sx q[0];
rz(0.31025525) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1290413) q[2];
sx q[2];
rz(-1.1079567) q[2];
sx q[2];
rz(-1.3430581) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1923101) q[1];
sx q[1];
rz(-2.3282166) q[1];
sx q[1];
rz(-1.1622914) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6799029) q[3];
sx q[3];
rz(-0.91489313) q[3];
sx q[3];
rz(-2.4885268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.85764641) q[2];
sx q[2];
rz(-1.5479156) q[2];
sx q[2];
rz(2.471981) q[2];
rz(-0.56704632) q[3];
sx q[3];
rz(-1.0457057) q[3];
sx q[3];
rz(1.4001575) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6794716) q[0];
sx q[0];
rz(-1.3422048) q[0];
sx q[0];
rz(-2.2456428) q[0];
rz(-0.040291928) q[1];
sx q[1];
rz(-0.94172421) q[1];
sx q[1];
rz(1.5641854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1791081) q[0];
sx q[0];
rz(-0.024027457) q[0];
sx q[0];
rz(1.7683309) q[0];
rz(-pi) q[1];
rz(1.1383923) q[2];
sx q[2];
rz(-0.93195019) q[2];
sx q[2];
rz(2.2496719) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60301723) q[1];
sx q[1];
rz(-1.9780206) q[1];
sx q[1];
rz(-1.0124221) q[1];
x q[2];
rz(1.0684929) q[3];
sx q[3];
rz(-2.4728259) q[3];
sx q[3];
rz(2.7947102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8335874) q[2];
sx q[2];
rz(-0.31978017) q[2];
sx q[2];
rz(1.8316815) q[2];
rz(-1.9561907) q[3];
sx q[3];
rz(-2.2535321) q[3];
sx q[3];
rz(1.5230702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1401405) q[0];
sx q[0];
rz(-1.8403213) q[0];
sx q[0];
rz(2.1885827) q[0];
rz(-0.078527191) q[1];
sx q[1];
rz(-2.0546866) q[1];
sx q[1];
rz(2.3436117) q[1];
rz(-0.25987249) q[2];
sx q[2];
rz(-0.91283471) q[2];
sx q[2];
rz(1.2531493) q[2];
rz(-0.6497186) q[3];
sx q[3];
rz(-0.25593723) q[3];
sx q[3];
rz(-1.2091936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
