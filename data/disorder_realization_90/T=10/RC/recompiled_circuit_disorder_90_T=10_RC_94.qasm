OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(3.0545711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10843033) q[0];
sx q[0];
rz(-0.96956367) q[0];
sx q[0];
rz(0.45494672) q[0];
rz(-3.0084228) q[2];
sx q[2];
rz(-0.77169092) q[2];
sx q[2];
rz(2.0193677) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4359856) q[1];
sx q[1];
rz(-1.4325303) q[1];
sx q[1];
rz(-2.0558753) q[1];
x q[2];
rz(-0.41214715) q[3];
sx q[3];
rz(-1.6133568) q[3];
sx q[3];
rz(-1.3679078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0027851) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(0.37386093) q[2];
rz(0.3368245) q[3];
sx q[3];
rz(-1.5461494) q[3];
sx q[3];
rz(-0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8841298) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.9467547) q[0];
rz(3.0589814) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(-3.1412178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1966239) q[0];
sx q[0];
rz(-1.9163016) q[0];
sx q[0];
rz(-1.1254805) q[0];
rz(-2.2194901) q[2];
sx q[2];
rz(-0.32425913) q[2];
sx q[2];
rz(2.1622554) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.61987585) q[1];
sx q[1];
rz(-0.65345018) q[1];
sx q[1];
rz(2.7481649) q[1];
x q[2];
rz(1.5156636) q[3];
sx q[3];
rz(-2.030636) q[3];
sx q[3];
rz(-2.9210747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.80883819) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(0.17671281) q[2];
rz(-0.79408944) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(-0.55364048) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90213838) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(2.7368271) q[0];
rz(1.2813214) q[1];
sx q[1];
rz(-2.5679913) q[1];
sx q[1];
rz(-1.3084897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1256516) q[0];
sx q[0];
rz(-2.0309272) q[0];
sx q[0];
rz(-2.5815651) q[0];
rz(-2.6042125) q[2];
sx q[2];
rz(-1.3668622) q[2];
sx q[2];
rz(2.1539719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3209826) q[1];
sx q[1];
rz(-2.2690986) q[1];
sx q[1];
rz(0.33336158) q[1];
rz(-2.7029196) q[3];
sx q[3];
rz(-1.6524501) q[3];
sx q[3];
rz(-2.4077533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9222766) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(-0.02040872) q[2];
rz(-2.0698047) q[3];
sx q[3];
rz(-1.1276378) q[3];
sx q[3];
rz(2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(-2.2256057) q[0];
rz(-0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(1.0571009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0650478) q[0];
sx q[0];
rz(-1.1335982) q[0];
sx q[0];
rz(-1.5989499) q[0];
rz(2.1787203) q[2];
sx q[2];
rz(-2.0553556) q[2];
sx q[2];
rz(-1.0819266) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.033356655) q[1];
sx q[1];
rz(-1.5844966) q[1];
sx q[1];
rz(-1.800888) q[1];
rz(-1.7980868) q[3];
sx q[3];
rz(-2.3152581) q[3];
sx q[3];
rz(-2.8518761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2477734) q[2];
sx q[2];
rz(-0.63988581) q[2];
sx q[2];
rz(-2.7988953) q[2];
rz(1.6977067) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(2.4263583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7856359) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(-3.0551531) q[0];
rz(-1.3899639) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(-0.46868971) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2738004) q[0];
sx q[0];
rz(-0.81875728) q[0];
sx q[0];
rz(-1.0206945) q[0];
rz(-pi) q[1];
rz(-2.6195171) q[2];
sx q[2];
rz(-0.35208382) q[2];
sx q[2];
rz(-0.075369518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9710755) q[1];
sx q[1];
rz(-1.7726388) q[1];
sx q[1];
rz(1.438373) q[1];
rz(0.080805578) q[3];
sx q[3];
rz(-2.3256133) q[3];
sx q[3];
rz(-0.24085837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.1420574) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(0.86432499) q[2];
rz(2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(-0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7140759) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(-0.19700225) q[0];
rz(-1.3621832) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(2.8053455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4649481) q[0];
sx q[0];
rz(-1.1332382) q[0];
sx q[0];
rz(-0.69533555) q[0];
rz(-0.7746081) q[2];
sx q[2];
rz(-2.3850072) q[2];
sx q[2];
rz(-2.2677383) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0974478) q[1];
sx q[1];
rz(-0.98173414) q[1];
sx q[1];
rz(-1.6277905) q[1];
rz(-pi) q[2];
rz(-1.398596) q[3];
sx q[3];
rz(-2.2599054) q[3];
sx q[3];
rz(0.25758753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(-2.2952648) q[2];
rz(-1.9019295) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7476615) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(2.8826662) q[0];
rz(-1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-2.0475725) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31539311) q[0];
sx q[0];
rz(-2.0059735) q[0];
sx q[0];
rz(0.3776334) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3482434) q[2];
sx q[2];
rz(-1.5828653) q[2];
sx q[2];
rz(-2.0633069) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3257521) q[1];
sx q[1];
rz(-3.0430803) q[1];
sx q[1];
rz(-1.2835531) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.098071531) q[3];
sx q[3];
rz(-1.5489998) q[3];
sx q[3];
rz(-1.9779713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9672433) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(-0.32067498) q[2];
rz(0.43618068) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(-0.44803739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86673474) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(1.0048237) q[0];
rz(-2.5514305) q[1];
sx q[1];
rz(-0.90548038) q[1];
sx q[1];
rz(-0.80387962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76669508) q[0];
sx q[0];
rz(-2.156213) q[0];
sx q[0];
rz(-1.5234408) q[0];
rz(-pi) q[1];
rz(0.29381807) q[2];
sx q[2];
rz(-0.57379469) q[2];
sx q[2];
rz(2.3125355) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4623973) q[1];
sx q[1];
rz(-3.052366) q[1];
sx q[1];
rz(-1.6855082) q[1];
rz(-2.478066) q[3];
sx q[3];
rz(-1.6420206) q[3];
sx q[3];
rz(-2.0957259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8089495) q[2];
sx q[2];
rz(-2.2968764) q[2];
sx q[2];
rz(2.1311029) q[2];
rz(2.5381952) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5378961) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(0.4075152) q[0];
rz(-2.852476) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(0.75072748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.663561) q[0];
sx q[0];
rz(-0.41304092) q[0];
sx q[0];
rz(-0.77126276) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6523916) q[2];
sx q[2];
rz(-1.2921234) q[2];
sx q[2];
rz(2.1747053) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0093065) q[1];
sx q[1];
rz(-0.3711776) q[1];
sx q[1];
rz(2.1585141) q[1];
x q[2];
rz(1.5969855) q[3];
sx q[3];
rz(-1.8720172) q[3];
sx q[3];
rz(-1.477004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0728545) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(1.9753974) q[2];
rz(1.7769622) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(-0.021818074) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5831379) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(1.7512084) q[0];
rz(-0.33069262) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.4601382) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6351785) q[0];
sx q[0];
rz(-1.1513452) q[0];
sx q[0];
rz(2.123453) q[0];
x q[1];
rz(-0.023256217) q[2];
sx q[2];
rz(-1.0668945) q[2];
sx q[2];
rz(-2.0397253) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3223443) q[1];
sx q[1];
rz(-1.3867612) q[1];
sx q[1];
rz(-0.40477246) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7695997) q[3];
sx q[3];
rz(-1.0255314) q[3];
sx q[3];
rz(1.256497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7918487) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(1.6798518) q[2];
rz(1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(0.45564836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6256975) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(-1.3810146) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(-2.7881651) q[2];
sx q[2];
rz(-2.3218487) q[2];
sx q[2];
rz(-0.29816366) q[2];
rz(-1.9588884) q[3];
sx q[3];
rz(-2.2065065) q[3];
sx q[3];
rz(1.828215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];