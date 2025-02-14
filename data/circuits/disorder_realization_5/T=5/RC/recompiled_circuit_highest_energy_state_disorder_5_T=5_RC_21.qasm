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
rz(-2.8992803) q[0];
sx q[0];
rz(-0.87376142) q[0];
sx q[0];
rz(-1.4886966) q[0];
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
rz(-1.0996286) q[0];
sx q[0];
rz(-1.4587741) q[0];
sx q[0];
rz(0.69667453) q[0];
rz(-pi) q[1];
rz(-0.27973819) q[2];
sx q[2];
rz(-1.7953897) q[2];
sx q[2];
rz(2.7034134) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1716071) q[1];
sx q[1];
rz(-2.0678221) q[1];
sx q[1];
rz(-2.7812048) q[1];
rz(-pi) q[2];
rz(2.4116573) q[3];
sx q[3];
rz(-1.1570108) q[3];
sx q[3];
rz(-2.1073996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0327518) q[2];
sx q[2];
rz(-2.8665172) q[2];
sx q[2];
rz(1.3230422) q[2];
rz(0.9465341) q[3];
sx q[3];
rz(-2.5338379) q[3];
sx q[3];
rz(3.0834468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35689795) q[0];
sx q[0];
rz(-0.31814831) q[0];
sx q[0];
rz(1.963105) q[0];
rz(2.5937041) q[1];
sx q[1];
rz(-1.9375487) q[1];
sx q[1];
rz(-1.9080124) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2425893) q[0];
sx q[0];
rz(-1.543643) q[0];
sx q[0];
rz(1.8544556) q[0];
rz(-pi) q[1];
rz(3.060922) q[2];
sx q[2];
rz(-1.4051132) q[2];
sx q[2];
rz(-2.5625429) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4405304) q[1];
sx q[1];
rz(-0.077210285) q[1];
sx q[1];
rz(-2.3838897) q[1];
rz(-pi) q[2];
rz(-1.3235533) q[3];
sx q[3];
rz(-1.2359058) q[3];
sx q[3];
rz(3.0680198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2445688) q[2];
sx q[2];
rz(-2.343101) q[2];
sx q[2];
rz(1.3429674) q[2];
rz(3.1161984) q[3];
sx q[3];
rz(-1.7829021) q[3];
sx q[3];
rz(-3.0620745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.9459166) q[0];
sx q[0];
rz(-1.8657277) q[0];
sx q[0];
rz(-0.50862408) q[0];
rz(1.7018082) q[1];
sx q[1];
rz(-1.6248117) q[1];
sx q[1];
rz(1.6280828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17643988) q[0];
sx q[0];
rz(-2.5178268) q[0];
sx q[0];
rz(-0.2709312) q[0];
rz(-pi) q[1];
rz(2.2824487) q[2];
sx q[2];
rz(-0.22822389) q[2];
sx q[2];
rz(1.3068401) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.95437223) q[1];
sx q[1];
rz(-0.037211671) q[1];
sx q[1];
rz(3.0195435) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9779889) q[3];
sx q[3];
rz(-1.9137205) q[3];
sx q[3];
rz(2.4437827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14579183) q[2];
sx q[2];
rz(-2.3172816) q[2];
sx q[2];
rz(-2.9902747) q[2];
rz(0.96210903) q[3];
sx q[3];
rz(-1.5113219) q[3];
sx q[3];
rz(0.37402737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91604084) q[0];
sx q[0];
rz(-2.2489838) q[0];
sx q[0];
rz(1.4501866) q[0];
rz(2.080503) q[1];
sx q[1];
rz(-2.7174945) q[1];
sx q[1];
rz(1.1294956) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6442959) q[0];
sx q[0];
rz(-0.82312778) q[0];
sx q[0];
rz(-1.3038365) q[0];
rz(-pi) q[1];
rz(2.5322459) q[2];
sx q[2];
rz(-2.5188392) q[2];
sx q[2];
rz(-1.8660752) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7348014) q[1];
sx q[1];
rz(-2.2815721) q[1];
sx q[1];
rz(-0.18376499) q[1];
rz(-1.4598386) q[3];
sx q[3];
rz(-1.5361388) q[3];
sx q[3];
rz(-2.2269888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9811161) q[2];
sx q[2];
rz(-2.1500197) q[2];
sx q[2];
rz(-1.1264832) q[2];
rz(0.22731656) q[3];
sx q[3];
rz(-1.7132297) q[3];
sx q[3];
rz(-0.050392438) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.324976) q[0];
sx q[0];
rz(-2.3351674) q[0];
sx q[0];
rz(0.6977914) q[0];
rz(-1.5296439) q[1];
sx q[1];
rz(-2.7841214) q[1];
sx q[1];
rz(-0.43561092) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742201) q[0];
sx q[0];
rz(-0.29074235) q[0];
sx q[0];
rz(-2.0176397) q[0];
rz(-0.37211352) q[2];
sx q[2];
rz(-0.96192322) q[2];
sx q[2];
rz(0.13932589) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2667234) q[1];
sx q[1];
rz(-1.2064955) q[1];
sx q[1];
rz(1.2357388) q[1];
rz(-pi) q[2];
rz(2.8661679) q[3];
sx q[3];
rz(-1.4990303) q[3];
sx q[3];
rz(-2.3022918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39516732) q[2];
sx q[2];
rz(-2.1992511) q[2];
sx q[2];
rz(0.64685267) q[2];
rz(-0.70710373) q[3];
sx q[3];
rz(-3.0349858) q[3];
sx q[3];
rz(-2.148237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7917204) q[0];
sx q[0];
rz(-2.7339022) q[0];
sx q[0];
rz(0.83734018) q[0];
rz(-2.1912241) q[1];
sx q[1];
rz(-1.8888387) q[1];
sx q[1];
rz(2.9577435) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0113676) q[0];
sx q[0];
rz(-1.8357903) q[0];
sx q[0];
rz(-0.30512198) q[0];
rz(1.2410795) q[2];
sx q[2];
rz(-1.7779254) q[2];
sx q[2];
rz(-1.1348534) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5279416) q[1];
sx q[1];
rz(-1.9518449) q[1];
sx q[1];
rz(-2.3697229) q[1];
x q[2];
rz(-0.23294874) q[3];
sx q[3];
rz(-1.4714575) q[3];
sx q[3];
rz(0.24391045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4652555) q[2];
sx q[2];
rz(-2.4691171) q[2];
sx q[2];
rz(-2.903751) q[2];
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
rz(-pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38943648) q[0];
sx q[0];
rz(-0.51920033) q[0];
sx q[0];
rz(3.0160548) q[0];
rz(-0.43352661) q[1];
sx q[1];
rz(-2.5018689) q[1];
sx q[1];
rz(-1.7009521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88461965) q[0];
sx q[0];
rz(-0.077736028) q[0];
sx q[0];
rz(0.90330584) q[0];
rz(-0.48831482) q[2];
sx q[2];
rz(-1.7124694) q[2];
sx q[2];
rz(-0.087208793) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0338965) q[1];
sx q[1];
rz(-2.0667065) q[1];
sx q[1];
rz(2.0350083) q[1];
rz(-pi) q[2];
rz(-0.85236211) q[3];
sx q[3];
rz(-0.96489513) q[3];
sx q[3];
rz(2.5123799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10616779) q[2];
sx q[2];
rz(-1.9751534) q[2];
sx q[2];
rz(-2.9336145) q[2];
rz(-0.68317991) q[3];
sx q[3];
rz(-0.71931374) q[3];
sx q[3];
rz(-1.6600367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55027562) q[0];
sx q[0];
rz(-2.2324201) q[0];
sx q[0];
rz(-2.7081178) q[0];
rz(-2.6705961) q[1];
sx q[1];
rz(-2.8484671) q[1];
sx q[1];
rz(-2.8004144) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0302673) q[0];
sx q[0];
rz(-1.4864392) q[0];
sx q[0];
rz(-1.4813759) q[0];
rz(2.3287235) q[2];
sx q[2];
rz(-1.4238266) q[2];
sx q[2];
rz(-0.16260553) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.78032683) q[1];
sx q[1];
rz(-1.396679) q[1];
sx q[1];
rz(2.1760015) q[1];
rz(-pi) q[2];
rz(-1.6834931) q[3];
sx q[3];
rz(-0.98033614) q[3];
sx q[3];
rz(1.0952418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6097673) q[2];
sx q[2];
rz(-2.7957081) q[2];
sx q[2];
rz(-2.7609265) q[2];
rz(2.5139513) q[3];
sx q[3];
rz(-2.6142879) q[3];
sx q[3];
rz(0.92831534) q[3];
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
sx q[2];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7577359) q[0];
sx q[0];
rz(-2.3691874) q[0];
sx q[0];
rz(0.17018233) q[0];
x q[1];
rz(2.6128255) q[2];
sx q[2];
rz(-1.6739539) q[2];
sx q[2];
rz(0.29036748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.57643465) q[1];
sx q[1];
rz(-1.7243313) q[1];
sx q[1];
rz(1.9839194) q[1];
rz(1.2282886) q[3];
sx q[3];
rz(-2.0039903) q[3];
sx q[3];
rz(1.9613645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7696441) q[2];
sx q[2];
rz(-0.19522788) q[2];
sx q[2];
rz(-0.39607421) q[2];
rz(1.1560446) q[3];
sx q[3];
rz(-2.5542185) q[3];
sx q[3];
rz(0.42266947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7523772) q[0];
sx q[0];
rz(-0.14425819) q[0];
sx q[0];
rz(0.17917646) q[0];
rz(-1.3009118) q[1];
sx q[1];
rz(-0.4549028) q[1];
sx q[1];
rz(0.5295583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1269187) q[0];
sx q[0];
rz(-1.8243941) q[0];
sx q[0];
rz(0.34832592) q[0];
rz(1.5949202) q[2];
sx q[2];
rz(-2.8608426) q[2];
sx q[2];
rz(2.8737673) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22687616) q[1];
sx q[1];
rz(-2.2351932) q[1];
sx q[1];
rz(1.2032582) q[1];
rz(2.6050636) q[3];
sx q[3];
rz(-1.7236815) q[3];
sx q[3];
rz(-2.9067723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.788488) q[2];
sx q[2];
rz(-2.4994734) q[2];
sx q[2];
rz(-2.8118964) q[2];
rz(-1.6126136) q[3];
sx q[3];
rz(-1.2639045) q[3];
sx q[3];
rz(-0.32040709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0223087) q[0];
sx q[0];
rz(-1.4531463) q[0];
sx q[0];
rz(-1.4896738) q[0];
rz(0.48619167) q[1];
sx q[1];
rz(-1.9959027) q[1];
sx q[1];
rz(-2.1015658) q[1];
rz(0.6766161) q[2];
sx q[2];
rz(-1.6634273) q[2];
sx q[2];
rz(1.5469748) q[2];
rz(-0.33489758) q[3];
sx q[3];
rz(-0.45968243) q[3];
sx q[3];
rz(-2.2062306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
