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
rz(-2.2678312) q[0];
sx q[0];
rz(1.4886966) q[0];
rz(1.0701264) q[1];
sx q[1];
rz(-2.5843599) q[1];
sx q[1];
rz(-0.4692404) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5375929) q[0];
sx q[0];
rz(-2.4374593) q[0];
sx q[0];
rz(2.9680433) q[0];
rz(-pi) q[1];
rz(-0.27973819) q[2];
sx q[2];
rz(-1.7953897) q[2];
sx q[2];
rz(-0.43817929) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8404313) q[1];
sx q[1];
rz(-0.60493785) q[1];
sx q[1];
rz(0.99436064) q[1];
x q[2];
rz(0.72993536) q[3];
sx q[3];
rz(-1.1570108) q[3];
sx q[3];
rz(-1.034193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0327518) q[2];
sx q[2];
rz(-0.27507541) q[2];
sx q[2];
rz(-1.3230422) q[2];
rz(2.1950586) q[3];
sx q[3];
rz(-0.60775477) q[3];
sx q[3];
rz(3.0834468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35689795) q[0];
sx q[0];
rz(-2.8234443) q[0];
sx q[0];
rz(1.963105) q[0];
rz(0.54788852) q[1];
sx q[1];
rz(-1.204044) q[1];
sx q[1];
rz(1.2335802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57892859) q[0];
sx q[0];
rz(-2.8566716) q[0];
sx q[0];
rz(1.4740545) q[0];
rz(-0.080670653) q[2];
sx q[2];
rz(-1.7364795) q[2];
sx q[2];
rz(-0.57904977) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4602549) q[1];
sx q[1];
rz(-1.6268568) q[1];
sx q[1];
rz(-1.6239151) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3235533) q[3];
sx q[3];
rz(-1.9056869) q[3];
sx q[3];
rz(3.0680198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8970238) q[2];
sx q[2];
rz(-2.343101) q[2];
sx q[2];
rz(-1.7986253) q[2];
rz(3.1161984) q[3];
sx q[3];
rz(-1.7829021) q[3];
sx q[3];
rz(0.079518147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19567604) q[0];
sx q[0];
rz(-1.8657277) q[0];
sx q[0];
rz(-2.6329686) q[0];
rz(-1.4397844) q[1];
sx q[1];
rz(-1.516781) q[1];
sx q[1];
rz(1.5135099) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17643988) q[0];
sx q[0];
rz(-2.5178268) q[0];
sx q[0];
rz(2.8706615) q[0];
rz(-2.9910477) q[2];
sx q[2];
rz(-1.398613) q[2];
sx q[2];
rz(-0.5821705) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.403203) q[1];
sx q[1];
rz(-1.5753257) q[1];
sx q[1];
rz(-3.1046575) q[1];
rz(0.37083309) q[3];
sx q[3];
rz(-1.1885839) q[3];
sx q[3];
rz(-2.1245898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9958008) q[2];
sx q[2];
rz(-2.3172816) q[2];
sx q[2];
rz(0.15131797) q[2];
rz(-2.1794836) q[3];
sx q[3];
rz(-1.6302707) q[3];
sx q[3];
rz(2.7675653) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91604084) q[0];
sx q[0];
rz(-0.89260888) q[0];
sx q[0];
rz(1.4501866) q[0];
rz(-1.0610896) q[1];
sx q[1];
rz(-2.7174945) q[1];
sx q[1];
rz(1.1294956) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8842402) q[0];
sx q[0];
rz(-1.7654618) q[0];
sx q[0];
rz(-2.3759222) q[0];
rz(1.9607294) q[2];
sx q[2];
rz(-1.0720813) q[2];
sx q[2];
rz(2.5759144) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7348014) q[1];
sx q[1];
rz(-0.86002058) q[1];
sx q[1];
rz(2.9578277) q[1];
rz(1.4598386) q[3];
sx q[3];
rz(-1.5361388) q[3];
sx q[3];
rz(2.2269888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.16047655) q[2];
sx q[2];
rz(-2.1500197) q[2];
sx q[2];
rz(2.0151095) q[2];
rz(2.9142761) q[3];
sx q[3];
rz(-1.7132297) q[3];
sx q[3];
rz(-3.0912002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.324976) q[0];
sx q[0];
rz(-0.80642527) q[0];
sx q[0];
rz(-0.6977914) q[0];
rz(1.5296439) q[1];
sx q[1];
rz(-2.7841214) q[1];
sx q[1];
rz(-2.7059817) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067372539) q[0];
sx q[0];
rz(-2.8508503) q[0];
sx q[0];
rz(2.0176397) q[0];
rz(2.0514652) q[2];
sx q[2];
rz(-0.70106912) q[2];
sx q[2];
rz(0.73814476) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.322239) q[1];
sx q[1];
rz(-1.8830918) q[1];
sx q[1];
rz(0.38375591) q[1];
x q[2];
rz(-2.8831611) q[3];
sx q[3];
rz(-0.28439097) q[3];
sx q[3];
rz(-2.1616158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7464253) q[2];
sx q[2];
rz(-0.94234157) q[2];
sx q[2];
rz(-0.64685267) q[2];
rz(-2.4344889) q[3];
sx q[3];
rz(-0.10660684) q[3];
sx q[3];
rz(0.99335563) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7917204) q[0];
sx q[0];
rz(-0.40769044) q[0];
sx q[0];
rz(-0.83734018) q[0];
rz(0.95036858) q[1];
sx q[1];
rz(-1.252754) q[1];
sx q[1];
rz(-2.9577435) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1302251) q[0];
sx q[0];
rz(-1.3058024) q[0];
sx q[0];
rz(-2.8364707) q[0];
rz(-0.21855905) q[2];
sx q[2];
rz(-1.248385) q[2];
sx q[2];
rz(0.36568322) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5317418) q[1];
sx q[1];
rz(-2.27503) q[1];
sx q[1];
rz(-1.0610214) q[1];
rz(-pi) q[2];
rz(-0.40756837) q[3];
sx q[3];
rz(-2.8887013) q[3];
sx q[3];
rz(1.7228218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6763372) q[2];
sx q[2];
rz(-2.4691171) q[2];
sx q[2];
rz(0.23784168) q[2];
rz(3.0658718) q[3];
sx q[3];
rz(-0.12756158) q[3];
sx q[3];
rz(-0.41847509) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38943648) q[0];
sx q[0];
rz(-2.6223923) q[0];
sx q[0];
rz(3.0160548) q[0];
rz(-0.43352661) q[1];
sx q[1];
rz(-2.5018689) q[1];
sx q[1];
rz(1.4406406) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.256973) q[0];
sx q[0];
rz(-0.077736028) q[0];
sx q[0];
rz(-2.2382868) q[0];
x q[1];
rz(-0.29514298) q[2];
sx q[2];
rz(-2.6347403) q[2];
sx q[2];
rz(1.7433804) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.444601) q[1];
sx q[1];
rz(-1.1660514) q[1];
sx q[1];
rz(2.5974817) q[1];
rz(-pi) q[2];
rz(2.3818199) q[3];
sx q[3];
rz(-0.90355146) q[3];
sx q[3];
rz(1.5189288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10616779) q[2];
sx q[2];
rz(-1.9751534) q[2];
sx q[2];
rz(0.20797813) q[2];
rz(2.4584127) q[3];
sx q[3];
rz(-0.71931374) q[3];
sx q[3];
rz(1.4815559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.591317) q[0];
sx q[0];
rz(-2.2324201) q[0];
sx q[0];
rz(-2.7081178) q[0];
rz(0.47099653) q[1];
sx q[1];
rz(-2.8484671) q[1];
sx q[1];
rz(-2.8004144) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0302673) q[0];
sx q[0];
rz(-1.4864392) q[0];
sx q[0];
rz(-1.4813759) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3287235) q[2];
sx q[2];
rz(-1.717766) q[2];
sx q[2];
rz(2.9789871) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5964929) q[1];
sx q[1];
rz(-0.6267281) q[1];
sx q[1];
rz(1.2709446) q[1];
rz(-pi) q[2];
rz(0.16624404) q[3];
sx q[3];
rz(-0.59985929) q[3];
sx q[3];
rz(-1.8458008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6097673) q[2];
sx q[2];
rz(-0.34588459) q[2];
sx q[2];
rz(2.7609265) q[2];
rz(2.5139513) q[3];
sx q[3];
rz(-0.52730477) q[3];
sx q[3];
rz(2.2132773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8220383) q[0];
sx q[0];
rz(-1.1975937) q[0];
sx q[0];
rz(0.30617103) q[0];
rz(0.51433688) q[1];
sx q[1];
rz(-2.3898333) q[1];
sx q[1];
rz(2.9591282) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3838568) q[0];
sx q[0];
rz(-0.77240521) q[0];
sx q[0];
rz(-0.17018233) q[0];
rz(-pi) q[1];
rz(1.4514662) q[2];
sx q[2];
rz(-2.0964604) q[2];
sx q[2];
rz(-1.801072) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.565158) q[1];
sx q[1];
rz(-1.7243313) q[1];
sx q[1];
rz(-1.9839194) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62807758) q[3];
sx q[3];
rz(-2.5961743) q[3];
sx q[3];
rz(-1.884383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.3719486) q[2];
sx q[2];
rz(-0.19522788) q[2];
sx q[2];
rz(2.7455184) q[2];
rz(-1.9855481) q[3];
sx q[3];
rz(-2.5542185) q[3];
sx q[3];
rz(0.42266947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3892155) q[0];
sx q[0];
rz(-2.9973345) q[0];
sx q[0];
rz(-0.17917646) q[0];
rz(-1.3009118) q[1];
sx q[1];
rz(-0.4549028) q[1];
sx q[1];
rz(0.5295583) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4946154) q[0];
sx q[0];
rz(-1.2340607) q[0];
sx q[0];
rz(-1.8398477) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5466724) q[2];
sx q[2];
rz(-0.2807501) q[2];
sx q[2];
rz(-2.8737673) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5770108) q[1];
sx q[1];
rz(-1.2839926) q[1];
sx q[1];
rz(2.4433892) q[1];
x q[2];
rz(2.6050636) q[3];
sx q[3];
rz(-1.4179112) q[3];
sx q[3];
rz(2.9067723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.788488) q[2];
sx q[2];
rz(-2.4994734) q[2];
sx q[2];
rz(2.8118964) q[2];
rz(-1.5289791) q[3];
sx q[3];
rz(-1.8776882) q[3];
sx q[3];
rz(2.8211856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11928398) q[0];
sx q[0];
rz(-1.4531463) q[0];
sx q[0];
rz(-1.4896738) q[0];
rz(0.48619167) q[1];
sx q[1];
rz(-1.9959027) q[1];
sx q[1];
rz(-2.1015658) q[1];
rz(-0.6766161) q[2];
sx q[2];
rz(-1.4781654) q[2];
sx q[2];
rz(-1.5946178) q[2];
rz(2.7042403) q[3];
sx q[3];
rz(-1.4244546) q[3];
sx q[3];
rz(-0.33311346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
