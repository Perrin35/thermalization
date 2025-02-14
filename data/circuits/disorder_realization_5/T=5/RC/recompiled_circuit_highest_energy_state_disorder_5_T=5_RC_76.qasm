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
rz(-2.0714662) q[1];
sx q[1];
rz(-0.5572328) q[1];
sx q[1];
rz(-2.6723523) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0996286) q[0];
sx q[0];
rz(-1.4587741) q[0];
sx q[0];
rz(-0.69667453) q[0];
rz(-0.27973819) q[2];
sx q[2];
rz(-1.7953897) q[2];
sx q[2];
rz(-0.43817929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3011613) q[1];
sx q[1];
rz(-2.5366548) q[1];
sx q[1];
rz(0.99436064) q[1];
x q[2];
rz(-0.58236645) q[3];
sx q[3];
rz(-0.81985788) q[3];
sx q[3];
rz(-0.95891592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0327518) q[2];
sx q[2];
rz(-0.27507541) q[2];
sx q[2];
rz(-1.8185505) q[2];
rz(-2.1950586) q[3];
sx q[3];
rz(-2.5338379) q[3];
sx q[3];
rz(-0.058145903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35689795) q[0];
sx q[0];
rz(-0.31814831) q[0];
sx q[0];
rz(1.1784877) q[0];
rz(-0.54788852) q[1];
sx q[1];
rz(-1.204044) q[1];
sx q[1];
rz(1.9080124) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89900333) q[0];
sx q[0];
rz(-1.5979496) q[0];
sx q[0];
rz(1.287137) q[0];
rz(1.7370102) q[2];
sx q[2];
rz(-1.6503599) q[2];
sx q[2];
rz(-2.1631789) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.11352053) q[1];
sx q[1];
rz(-1.6238316) q[1];
sx q[1];
rz(-0.056139498) q[1];
rz(-pi) q[2];
rz(-0.34459262) q[3];
sx q[3];
rz(-1.8040415) q[3];
sx q[3];
rz(-1.4144554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2445688) q[2];
sx q[2];
rz(-2.343101) q[2];
sx q[2];
rz(1.3429674) q[2];
rz(-3.1161984) q[3];
sx q[3];
rz(-1.3586905) q[3];
sx q[3];
rz(-3.0620745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19567604) q[0];
sx q[0];
rz(-1.275865) q[0];
sx q[0];
rz(0.50862408) q[0];
rz(-1.4397844) q[1];
sx q[1];
rz(-1.516781) q[1];
sx q[1];
rz(-1.6280828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1726097) q[0];
sx q[0];
rz(-1.4138317) q[0];
sx q[0];
rz(0.60625546) q[0];
rz(-pi) q[1];
x q[1];
rz(0.150545) q[2];
sx q[2];
rz(-1.7429797) q[2];
sx q[2];
rz(0.5821705) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83223935) q[1];
sx q[1];
rz(-1.6077311) q[1];
sx q[1];
rz(-1.5753288) q[1];
rz(-pi) q[2];
rz(2.3044488) q[3];
sx q[3];
rz(-2.6154933) q[3];
sx q[3];
rz(0.21080454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14579183) q[2];
sx q[2];
rz(-0.82431102) q[2];
sx q[2];
rz(0.15131797) q[2];
rz(0.96210903) q[3];
sx q[3];
rz(-1.5113219) q[3];
sx q[3];
rz(0.37402737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91604084) q[0];
sx q[0];
rz(-0.89260888) q[0];
sx q[0];
rz(-1.4501866) q[0];
rz(-2.080503) q[1];
sx q[1];
rz(-2.7174945) q[1];
sx q[1];
rz(2.012097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2573525) q[0];
sx q[0];
rz(-1.7654618) q[0];
sx q[0];
rz(-2.3759222) q[0];
x q[1];
rz(-0.60934679) q[2];
sx q[2];
rz(-0.6227535) q[2];
sx q[2];
rz(-1.2755175) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0982588) q[1];
sx q[1];
rz(-1.4318649) q[1];
sx q[1];
rz(2.2899782) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6817541) q[3];
sx q[3];
rz(-1.5361388) q[3];
sx q[3];
rz(-2.2269888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.16047655) q[2];
sx q[2];
rz(-2.1500197) q[2];
sx q[2];
rz(2.0151095) q[2];
rz(-2.9142761) q[3];
sx q[3];
rz(-1.7132297) q[3];
sx q[3];
rz(-0.050392438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8166167) q[0];
sx q[0];
rz(-2.3351674) q[0];
sx q[0];
rz(2.4438013) q[0];
rz(1.5296439) q[1];
sx q[1];
rz(-2.7841214) q[1];
sx q[1];
rz(0.43561092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7451906) q[0];
sx q[0];
rz(-1.8322837) q[0];
sx q[0];
rz(3.0130062) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0514652) q[2];
sx q[2];
rz(-0.70106912) q[2];
sx q[2];
rz(-0.73814476) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81935362) q[1];
sx q[1];
rz(-1.2585009) q[1];
sx q[1];
rz(-0.38375591) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8661679) q[3];
sx q[3];
rz(-1.4990303) q[3];
sx q[3];
rz(-0.83930086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39516732) q[2];
sx q[2];
rz(-0.94234157) q[2];
sx q[2];
rz(0.64685267) q[2];
rz(0.70710373) q[3];
sx q[3];
rz(-0.10660684) q[3];
sx q[3];
rz(0.99335563) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3498722) q[0];
sx q[0];
rz(-0.40769044) q[0];
sx q[0];
rz(0.83734018) q[0];
rz(0.95036858) q[1];
sx q[1];
rz(-1.252754) q[1];
sx q[1];
rz(0.18384917) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0113676) q[0];
sx q[0];
rz(-1.8357903) q[0];
sx q[0];
rz(0.30512198) q[0];
rz(2.1464961) q[2];
sx q[2];
rz(-2.754236) q[2];
sx q[2];
rz(-0.97706276) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.819471) q[1];
sx q[1];
rz(-2.2985704) q[1];
sx q[1];
rz(-2.6202109) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4687187) q[3];
sx q[3];
rz(-1.3390171) q[3];
sx q[3];
rz(1.8382324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4652555) q[2];
sx q[2];
rz(-0.67247552) q[2];
sx q[2];
rz(0.23784168) q[2];
rz(0.075720876) q[3];
sx q[3];
rz(-3.0140311) q[3];
sx q[3];
rz(2.7231176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38943648) q[0];
sx q[0];
rz(-0.51920033) q[0];
sx q[0];
rz(-3.0160548) q[0];
rz(2.708066) q[1];
sx q[1];
rz(-2.5018689) q[1];
sx q[1];
rz(-1.7009521) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.256973) q[0];
sx q[0];
rz(-3.0638566) q[0];
sx q[0];
rz(0.90330584) q[0];
rz(-1.7309181) q[2];
sx q[2];
rz(-2.0537964) q[2];
sx q[2];
rz(-1.4087189) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4384253) q[1];
sx q[1];
rz(-0.66574186) q[1];
sx q[1];
rz(2.4502657) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3976753) q[3];
sx q[3];
rz(-2.1424681) q[3];
sx q[3];
rz(-0.47960348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10616779) q[2];
sx q[2];
rz(-1.9751534) q[2];
sx q[2];
rz(2.9336145) q[2];
rz(0.68317991) q[3];
sx q[3];
rz(-0.71931374) q[3];
sx q[3];
rz(1.6600367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.591317) q[0];
sx q[0];
rz(-0.90917259) q[0];
sx q[0];
rz(-2.7081178) q[0];
rz(2.6705961) q[1];
sx q[1];
rz(-0.29312557) q[1];
sx q[1];
rz(0.34117821) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46702532) q[0];
sx q[0];
rz(-1.4816947) q[0];
sx q[0];
rz(-0.084693935) q[0];
x q[1];
rz(0.20107953) q[2];
sx q[2];
rz(-0.82302457) q[2];
sx q[2];
rz(-1.2704598) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5964929) q[1];
sx q[1];
rz(-2.5148646) q[1];
sx q[1];
rz(1.2709446) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16624404) q[3];
sx q[3];
rz(-2.5417334) q[3];
sx q[3];
rz(1.2957919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.53182536) q[2];
sx q[2];
rz(-2.7957081) q[2];
sx q[2];
rz(0.38066614) q[2];
rz(-0.62764132) q[3];
sx q[3];
rz(-2.6142879) q[3];
sx q[3];
rz(0.92831534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8220383) q[0];
sx q[0];
rz(-1.943999) q[0];
sx q[0];
rz(-0.30617103) q[0];
rz(0.51433688) q[1];
sx q[1];
rz(-0.75175935) q[1];
sx q[1];
rz(0.18246442) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064474551) q[0];
sx q[0];
rz(-1.6892642) q[0];
sx q[0];
rz(-2.3764591) q[0];
rz(2.6128255) q[2];
sx q[2];
rz(-1.4676388) q[2];
sx q[2];
rz(2.8512252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.811552) q[1];
sx q[1];
rz(-0.43918681) q[1];
sx q[1];
rz(-1.9387092) q[1];
x q[2];
rz(2.6851607) q[3];
sx q[3];
rz(-1.2610389) q[3];
sx q[3];
rz(-2.8995958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.3719486) q[2];
sx q[2];
rz(-0.19522788) q[2];
sx q[2];
rz(2.7455184) q[2];
rz(1.1560446) q[3];
sx q[3];
rz(-0.58737415) q[3];
sx q[3];
rz(-0.42266947) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3892155) q[0];
sx q[0];
rz(-0.14425819) q[0];
sx q[0];
rz(-2.9624162) q[0];
rz(-1.3009118) q[1];
sx q[1];
rz(-2.6866899) q[1];
sx q[1];
rz(-0.5295583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64697726) q[0];
sx q[0];
rz(-1.2340607) q[0];
sx q[0];
rz(1.3017449) q[0];
rz(-pi) q[1];
rz(1.5949202) q[2];
sx q[2];
rz(-0.2807501) q[2];
sx q[2];
rz(0.26782537) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5770108) q[1];
sx q[1];
rz(-1.8576) q[1];
sx q[1];
rz(0.69820349) q[1];
rz(-pi) q[2];
rz(-0.53652904) q[3];
sx q[3];
rz(-1.4179112) q[3];
sx q[3];
rz(-0.23482032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.788488) q[2];
sx q[2];
rz(-2.4994734) q[2];
sx q[2];
rz(0.32969627) q[2];
rz(1.6126136) q[3];
sx q[3];
rz(-1.2639045) q[3];
sx q[3];
rz(-2.8211856) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0223087) q[0];
sx q[0];
rz(-1.4531463) q[0];
sx q[0];
rz(-1.4896738) q[0];
rz(-0.48619167) q[1];
sx q[1];
rz(-1.14569) q[1];
sx q[1];
rz(1.0400269) q[1];
rz(-2.4649766) q[2];
sx q[2];
rz(-1.6634273) q[2];
sx q[2];
rz(1.5469748) q[2];
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
