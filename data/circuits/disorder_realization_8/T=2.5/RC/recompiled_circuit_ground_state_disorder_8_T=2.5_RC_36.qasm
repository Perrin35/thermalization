OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4158004) q[0];
sx q[0];
rz(-1.9771201) q[0];
sx q[0];
rz(-1.8902984) q[0];
rz(2.8407821) q[1];
sx q[1];
rz(-1.7164813) q[1];
sx q[1];
rz(-3.0256174) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2271954) q[0];
sx q[0];
rz(-1.6956788) q[0];
sx q[0];
rz(-1.8950805) q[0];
x q[1];
rz(-2.0550344) q[2];
sx q[2];
rz(-1.3645483) q[2];
sx q[2];
rz(-2.1833486) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8311316) q[1];
sx q[1];
rz(-0.77378213) q[1];
sx q[1];
rz(1.4813862) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7099781) q[3];
sx q[3];
rz(-1.0324036) q[3];
sx q[3];
rz(2.0288426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3731602) q[2];
sx q[2];
rz(-3.0765522) q[2];
sx q[2];
rz(1.6981286) q[2];
rz(2.4453435) q[3];
sx q[3];
rz(-1.6257476) q[3];
sx q[3];
rz(-2.7878917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.26831216) q[0];
sx q[0];
rz(-0.074957632) q[0];
sx q[0];
rz(0.36175501) q[0];
rz(-1.7573645) q[1];
sx q[1];
rz(-2.7499966) q[1];
sx q[1];
rz(-1.0272383) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8105668) q[0];
sx q[0];
rz(-0.89052708) q[0];
sx q[0];
rz(0.23305889) q[0];
rz(-pi) q[1];
rz(-1.7387668) q[2];
sx q[2];
rz(-2.3609475) q[2];
sx q[2];
rz(1.7602242) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9578182) q[1];
sx q[1];
rz(-1.8190067) q[1];
sx q[1];
rz(-0.91285984) q[1];
x q[2];
rz(-1.6546685) q[3];
sx q[3];
rz(-1.7954651) q[3];
sx q[3];
rz(-0.39473907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0072173) q[2];
sx q[2];
rz(-0.5642887) q[2];
sx q[2];
rz(-1.0312274) q[2];
rz(-0.68650308) q[3];
sx q[3];
rz(-1.406823) q[3];
sx q[3];
rz(-2.8197207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.6512076) q[0];
sx q[0];
rz(-1.729916) q[0];
sx q[0];
rz(1.7899293) q[0];
rz(-1.5216113) q[1];
sx q[1];
rz(-0.23623513) q[1];
sx q[1];
rz(1.1865541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51031798) q[0];
sx q[0];
rz(-2.5532852) q[0];
sx q[0];
rz(-0.73008213) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6461238) q[2];
sx q[2];
rz(-1.1161303) q[2];
sx q[2];
rz(-2.1476434) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12315519) q[1];
sx q[1];
rz(-1.5644524) q[1];
sx q[1];
rz(3.1325587) q[1];
x q[2];
rz(-2.260798) q[3];
sx q[3];
rz(-1.8024416) q[3];
sx q[3];
rz(-0.85154545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97570193) q[2];
sx q[2];
rz(-1.2624792) q[2];
sx q[2];
rz(2.8603437) q[2];
rz(2.7665372) q[3];
sx q[3];
rz(-0.91383362) q[3];
sx q[3];
rz(2.8121172) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81958812) q[0];
sx q[0];
rz(-0.81398886) q[0];
sx q[0];
rz(-2.8170012) q[0];
rz(-1.548467) q[1];
sx q[1];
rz(-0.32061583) q[1];
sx q[1];
rz(0.89210192) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2587713) q[0];
sx q[0];
rz(-0.020892512) q[0];
sx q[0];
rz(2.5740795) q[0];
x q[1];
rz(-1.6692144) q[2];
sx q[2];
rz(-0.8618671) q[2];
sx q[2];
rz(2.6977328) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39246757) q[1];
sx q[1];
rz(-1.8294847) q[1];
sx q[1];
rz(1.2767747) q[1];
rz(-3.1079477) q[3];
sx q[3];
rz(-0.97760289) q[3];
sx q[3];
rz(-0.9481687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0804245) q[2];
sx q[2];
rz(-1.755038) q[2];
sx q[2];
rz(0.57677734) q[2];
rz(-0.12442496) q[3];
sx q[3];
rz(-2.3915229) q[3];
sx q[3];
rz(-0.09593825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4633923) q[0];
sx q[0];
rz(-0.41825774) q[0];
sx q[0];
rz(-0.28045714) q[0];
rz(-0.29817835) q[1];
sx q[1];
rz(-0.069998048) q[1];
sx q[1];
rz(2.9895474) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9579118) q[0];
sx q[0];
rz(-0.059558161) q[0];
sx q[0];
rz(-1.0899215) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0934257) q[2];
sx q[2];
rz(-1.4844456) q[2];
sx q[2];
rz(0.62355838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3144532) q[1];
sx q[1];
rz(-1.1734073) q[1];
sx q[1];
rz(2.023815) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5923813) q[3];
sx q[3];
rz(-0.38281554) q[3];
sx q[3];
rz(2.5408753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8449479) q[2];
sx q[2];
rz(-1.351202) q[2];
sx q[2];
rz(-3.1015934) q[2];
rz(0.1376888) q[3];
sx q[3];
rz(-2.6929893) q[3];
sx q[3];
rz(0.78534809) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8765151) q[0];
sx q[0];
rz(-3.0410933) q[0];
sx q[0];
rz(2.5608089) q[0];
rz(0.68897092) q[1];
sx q[1];
rz(-0.13801485) q[1];
sx q[1];
rz(1.2977915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.101872) q[0];
sx q[0];
rz(-0.61148171) q[0];
sx q[0];
rz(0.78973706) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6766729) q[2];
sx q[2];
rz(-0.54621241) q[2];
sx q[2];
rz(-0.27413163) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4124085) q[1];
sx q[1];
rz(-1.7910335) q[1];
sx q[1];
rz(-2.7413648) q[1];
rz(-pi) q[2];
rz(-2.0125603) q[3];
sx q[3];
rz(-1.0241177) q[3];
sx q[3];
rz(-2.2571079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7350404) q[2];
sx q[2];
rz(-1.8031969) q[2];
sx q[2];
rz(-1.5945386) q[2];
rz(0.43153396) q[3];
sx q[3];
rz(-1.1003234) q[3];
sx q[3];
rz(0.6507473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4887017) q[0];
sx q[0];
rz(-0.015691375) q[0];
sx q[0];
rz(2.5456862) q[0];
rz(-1.7911004) q[1];
sx q[1];
rz(-0.48521438) q[1];
sx q[1];
rz(0.57292604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8329215) q[0];
sx q[0];
rz(-3.1317319) q[0];
sx q[0];
rz(0.29247756) q[0];
x q[1];
rz(2.7155034) q[2];
sx q[2];
rz(-1.3182148) q[2];
sx q[2];
rz(1.9468228) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.10526722) q[1];
sx q[1];
rz(-2.46439) q[1];
sx q[1];
rz(-0.52097229) q[1];
x q[2];
rz(1.9455276) q[3];
sx q[3];
rz(-0.84991377) q[3];
sx q[3];
rz(1.8955613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4493745) q[2];
sx q[2];
rz(-2.2276679) q[2];
sx q[2];
rz(-2.2268028) q[2];
rz(1.6080914) q[3];
sx q[3];
rz(-1.1398818) q[3];
sx q[3];
rz(-2.1515414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6221301) q[0];
sx q[0];
rz(-1.7161481) q[0];
sx q[0];
rz(-3.07716) q[0];
rz(-2.875115) q[1];
sx q[1];
rz(-0.12424145) q[1];
sx q[1];
rz(1.2640094) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5910931) q[0];
sx q[0];
rz(-1.2471218) q[0];
sx q[0];
rz(1.8514015) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1830215) q[2];
sx q[2];
rz(-2.2335548) q[2];
sx q[2];
rz(0.52802625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.79962423) q[1];
sx q[1];
rz(-1.2576629) q[1];
sx q[1];
rz(1.3596441) q[1];
x q[2];
rz(-1.0059935) q[3];
sx q[3];
rz(-2.0801615) q[3];
sx q[3];
rz(-2.6107174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6182575) q[2];
sx q[2];
rz(-1.9674415) q[2];
sx q[2];
rz(2.0972283) q[2];
rz(2.4339645) q[3];
sx q[3];
rz(-2.7510567) q[3];
sx q[3];
rz(2.0948758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3656965) q[0];
sx q[0];
rz(-1.5702268) q[0];
sx q[0];
rz(-0.43029341) q[0];
rz(2.4039092) q[1];
sx q[1];
rz(-3.0569515) q[1];
sx q[1];
rz(2.784909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9707613) q[0];
sx q[0];
rz(-0.30730844) q[0];
sx q[0];
rz(2.8224562) q[0];
rz(-pi) q[1];
rz(-1.6820774) q[2];
sx q[2];
rz(-0.26646895) q[2];
sx q[2];
rz(2.1499718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1445279) q[1];
sx q[1];
rz(-2.4127919) q[1];
sx q[1];
rz(2.4841734) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5857626) q[3];
sx q[3];
rz(-0.14066108) q[3];
sx q[3];
rz(-2.6859796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.77013993) q[2];
sx q[2];
rz(-0.74718237) q[2];
sx q[2];
rz(-2.746197) q[2];
rz(1.4622408) q[3];
sx q[3];
rz(-1.3926927) q[3];
sx q[3];
rz(3.1126378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28703654) q[0];
sx q[0];
rz(-2.365878) q[0];
sx q[0];
rz(-2.4391644) q[0];
rz(2.9632945) q[1];
sx q[1];
rz(-0.65419227) q[1];
sx q[1];
rz(-1.6184695) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010945646) q[0];
sx q[0];
rz(-0.48640448) q[0];
sx q[0];
rz(-2.0964699) q[0];
x q[1];
rz(1.6687001) q[2];
sx q[2];
rz(-1.5563097) q[2];
sx q[2];
rz(1.9928428) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6206086) q[1];
sx q[1];
rz(-1.0080326) q[1];
sx q[1];
rz(-1.4042793) q[1];
rz(0.50242063) q[3];
sx q[3];
rz(-1.2278886) q[3];
sx q[3];
rz(-0.08758647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6452597) q[2];
sx q[2];
rz(-2.1767904) q[2];
sx q[2];
rz(2.1920835) q[2];
rz(2.005714) q[3];
sx q[3];
rz(-0.94616008) q[3];
sx q[3];
rz(0.56376636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62354325) q[0];
sx q[0];
rz(-1.8907056) q[0];
sx q[0];
rz(-0.84778669) q[0];
rz(-1.7060948) q[1];
sx q[1];
rz(-1.2789627) q[1];
sx q[1];
rz(0.36437558) q[1];
rz(1.4451531) q[2];
sx q[2];
rz(-2.4375765) q[2];
sx q[2];
rz(0.30839105) q[2];
rz(1.2632542) q[3];
sx q[3];
rz(-1.5123541) q[3];
sx q[3];
rz(-2.1337951) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
