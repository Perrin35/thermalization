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
rz(0.24231237) q[0];
sx q[0];
rz(4.0153541) q[0];
sx q[0];
rz(10.913475) q[0];
rz(1.0701264) q[1];
sx q[1];
rz(3.6988255) q[1];
sx q[1];
rz(8.9555376) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37791489) q[0];
sx q[0];
rz(-0.87935424) q[0];
sx q[0];
rz(1.4251644) q[0];
x q[1];
rz(0.69122423) q[2];
sx q[2];
rz(-2.7847053) q[2];
sx q[2];
rz(-1.7920902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8404313) q[1];
sx q[1];
rz(-2.5366548) q[1];
sx q[1];
rz(2.147232) q[1];
x q[2];
rz(1.038299) q[3];
sx q[3];
rz(-2.2275452) q[3];
sx q[3];
rz(-0.19124243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1088408) q[2];
sx q[2];
rz(-2.8665172) q[2];
sx q[2];
rz(-1.8185505) q[2];
rz(0.9465341) q[3];
sx q[3];
rz(-2.5338379) q[3];
sx q[3];
rz(-0.058145903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7846947) q[0];
sx q[0];
rz(-2.8234443) q[0];
sx q[0];
rz(-1.1784877) q[0];
rz(-0.54788852) q[1];
sx q[1];
rz(-1.9375487) q[1];
sx q[1];
rz(1.2335802) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6797076) q[0];
sx q[0];
rz(-1.8543482) q[0];
sx q[0];
rz(0.028282982) q[0];
rz(-pi) q[1];
rz(-1.1217246) q[2];
sx q[2];
rz(-2.9574759) q[2];
sx q[2];
rz(-1.0348233) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4405304) q[1];
sx q[1];
rz(-0.077210285) q[1];
sx q[1];
rz(2.3838897) q[1];
rz(-pi) q[2];
rz(-1.8180394) q[3];
sx q[3];
rz(-1.9056869) q[3];
sx q[3];
rz(3.0680198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2445688) q[2];
sx q[2];
rz(-2.343101) q[2];
sx q[2];
rz(1.7986253) q[2];
rz(3.1161984) q[3];
sx q[3];
rz(-1.7829021) q[3];
sx q[3];
rz(0.079518147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9459166) q[0];
sx q[0];
rz(-1.8657277) q[0];
sx q[0];
rz(-2.6329686) q[0];
rz(1.4397844) q[1];
sx q[1];
rz(-1.6248117) q[1];
sx q[1];
rz(-1.6280828) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17643988) q[0];
sx q[0];
rz(-2.5178268) q[0];
sx q[0];
rz(2.8706615) q[0];
rz(-pi) q[1];
x q[1];
rz(0.150545) q[2];
sx q[2];
rz(-1.7429797) q[2];
sx q[2];
rz(-2.5594222) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.73838961) q[1];
sx q[1];
rz(-1.5753257) q[1];
sx q[1];
rz(0.036935115) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1636038) q[3];
sx q[3];
rz(-1.9137205) q[3];
sx q[3];
rz(-0.69780998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.14579183) q[2];
sx q[2];
rz(-2.3172816) q[2];
sx q[2];
rz(-2.9902747) q[2];
rz(0.96210903) q[3];
sx q[3];
rz(-1.6302707) q[3];
sx q[3];
rz(2.7675653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91604084) q[0];
sx q[0];
rz(-2.2489838) q[0];
sx q[0];
rz(-1.691406) q[0];
rz(-2.080503) q[1];
sx q[1];
rz(-2.7174945) q[1];
sx q[1];
rz(2.012097) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8842402) q[0];
sx q[0];
rz(-1.3761308) q[0];
sx q[0];
rz(0.76567044) q[0];
x q[1];
rz(-1.9607294) q[2];
sx q[2];
rz(-2.0695114) q[2];
sx q[2];
rz(2.5759144) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.40679125) q[1];
sx q[1];
rz(-2.2815721) q[1];
sx q[1];
rz(0.18376499) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2673504) q[3];
sx q[3];
rz(-0.11622322) q[3];
sx q[3];
rz(0.35467142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9811161) q[2];
sx q[2];
rz(-2.1500197) q[2];
sx q[2];
rz(1.1264832) q[2];
rz(-0.22731656) q[3];
sx q[3];
rz(-1.428363) q[3];
sx q[3];
rz(3.0912002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.324976) q[0];
sx q[0];
rz(-2.3351674) q[0];
sx q[0];
rz(-0.6977914) q[0];
rz(-1.6119488) q[1];
sx q[1];
rz(-2.7841214) q[1];
sx q[1];
rz(-2.7059817) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9337845) q[0];
sx q[0];
rz(-1.446604) q[0];
sx q[0];
rz(-1.8343614) q[0];
x q[1];
rz(2.0514652) q[2];
sx q[2];
rz(-0.70106912) q[2];
sx q[2];
rz(-2.4034479) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0406094) q[1];
sx q[1];
rz(-0.48982778) q[1];
sx q[1];
rz(0.71160729) q[1];
rz(-pi) q[2];
rz(-2.8831611) q[3];
sx q[3];
rz(-0.28439097) q[3];
sx q[3];
rz(0.97997682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7464253) q[2];
sx q[2];
rz(-0.94234157) q[2];
sx q[2];
rz(-2.49474) q[2];
rz(-2.4344889) q[3];
sx q[3];
rz(-0.10660684) q[3];
sx q[3];
rz(-2.148237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3498722) q[0];
sx q[0];
rz(-0.40769044) q[0];
sx q[0];
rz(-0.83734018) q[0];
rz(2.1912241) q[1];
sx q[1];
rz(-1.8888387) q[1];
sx q[1];
rz(0.18384917) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35826916) q[0];
sx q[0];
rz(-1.8649432) q[0];
sx q[0];
rz(1.8479895) q[0];
rz(-1.9005132) q[2];
sx q[2];
rz(-1.7779254) q[2];
sx q[2];
rz(-1.1348534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5317418) q[1];
sx q[1];
rz(-2.27503) q[1];
sx q[1];
rz(2.0805712) q[1];
rz(-pi) q[2];
rz(1.4687187) q[3];
sx q[3];
rz(-1.3390171) q[3];
sx q[3];
rz(-1.8382324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6763372) q[2];
sx q[2];
rz(-2.4691171) q[2];
sx q[2];
rz(0.23784168) q[2];
rz(3.0658718) q[3];
sx q[3];
rz(-3.0140311) q[3];
sx q[3];
rz(0.41847509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7521562) q[0];
sx q[0];
rz(-2.6223923) q[0];
sx q[0];
rz(-0.12553781) q[0];
rz(-2.708066) q[1];
sx q[1];
rz(-2.5018689) q[1];
sx q[1];
rz(-1.4406406) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020156064) q[0];
sx q[0];
rz(-1.5227063) q[0];
sx q[0];
rz(1.6318956) q[0];
x q[1];
rz(-2.8464497) q[2];
sx q[2];
rz(-0.50685234) q[2];
sx q[2];
rz(1.7433804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.444601) q[1];
sx q[1];
rz(-1.9755413) q[1];
sx q[1];
rz(2.5974817) q[1];
x q[2];
rz(-2.3976753) q[3];
sx q[3];
rz(-2.1424681) q[3];
sx q[3];
rz(-2.6619892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10616779) q[2];
sx q[2];
rz(-1.1664392) q[2];
sx q[2];
rz(0.20797813) q[2];
rz(2.4584127) q[3];
sx q[3];
rz(-2.4222789) q[3];
sx q[3];
rz(-1.4815559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55027562) q[0];
sx q[0];
rz(-2.2324201) q[0];
sx q[0];
rz(0.43347484) q[0];
rz(-0.47099653) q[1];
sx q[1];
rz(-2.8484671) q[1];
sx q[1];
rz(-0.34117821) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6745673) q[0];
sx q[0];
rz(-1.4816947) q[0];
sx q[0];
rz(-3.0568987) q[0];
rz(-pi) q[1];
rz(2.3287235) q[2];
sx q[2];
rz(-1.4238266) q[2];
sx q[2];
rz(2.9789871) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9097482) q[1];
sx q[1];
rz(-2.1655786) q[1];
sx q[1];
rz(2.9308795) q[1];
x q[2];
rz(-0.59340684) q[3];
sx q[3];
rz(-1.6643502) q[3];
sx q[3];
rz(2.7289651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.53182536) q[2];
sx q[2];
rz(-2.7957081) q[2];
sx q[2];
rz(2.7609265) q[2];
rz(-0.62764132) q[3];
sx q[3];
rz(-2.6142879) q[3];
sx q[3];
rz(-2.2132773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31955433) q[0];
sx q[0];
rz(-1.1975937) q[0];
sx q[0];
rz(2.8354216) q[0];
rz(0.51433688) q[1];
sx q[1];
rz(-0.75175935) q[1];
sx q[1];
rz(0.18246442) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5222591) q[0];
sx q[0];
rz(-0.81237185) q[0];
sx q[0];
rz(-1.4072529) q[0];
x q[1];
rz(2.9391871) q[2];
sx q[2];
rz(-0.53779624) q[2];
sx q[2];
rz(1.1059731) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3300407) q[1];
sx q[1];
rz(-0.43918681) q[1];
sx q[1];
rz(-1.2028834) q[1];
rz(-pi) q[2];
rz(2.6851607) q[3];
sx q[3];
rz(-1.8805537) q[3];
sx q[3];
rz(2.8995958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3719486) q[2];
sx q[2];
rz(-2.9463648) q[2];
sx q[2];
rz(-2.7455184) q[2];
rz(1.9855481) q[3];
sx q[3];
rz(-2.5542185) q[3];
sx q[3];
rz(2.7189232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.7523772) q[0];
sx q[0];
rz(-0.14425819) q[0];
sx q[0];
rz(-0.17917646) q[0];
rz(1.3009118) q[1];
sx q[1];
rz(-2.6866899) q[1];
sx q[1];
rz(0.5295583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1269187) q[0];
sx q[0];
rz(-1.8243941) q[0];
sx q[0];
rz(2.7932667) q[0];
rz(-1.851469) q[2];
sx q[2];
rz(-1.5641128) q[2];
sx q[2];
rz(1.3261507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9147165) q[1];
sx q[1];
rz(-0.90639948) q[1];
sx q[1];
rz(-1.2032582) q[1];
rz(2.8488068) q[3];
sx q[3];
rz(-2.5857627) q[3];
sx q[3];
rz(-1.0852927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35310465) q[2];
sx q[2];
rz(-0.64211923) q[2];
sx q[2];
rz(-0.32969627) q[2];
rz(1.5289791) q[3];
sx q[3];
rz(-1.8776882) q[3];
sx q[3];
rz(0.32040709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0223087) q[0];
sx q[0];
rz(-1.6884463) q[0];
sx q[0];
rz(1.6519188) q[0];
rz(-0.48619167) q[1];
sx q[1];
rz(-1.14569) q[1];
sx q[1];
rz(1.0400269) q[1];
rz(-1.4522105) q[2];
sx q[2];
rz(-0.89761843) q[2];
sx q[2];
rz(0.050326025) q[2];
rz(2.8066951) q[3];
sx q[3];
rz(-0.45968243) q[3];
sx q[3];
rz(-2.2062306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
