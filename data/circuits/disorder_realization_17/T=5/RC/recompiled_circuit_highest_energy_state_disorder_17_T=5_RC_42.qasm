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
rz(2.0789335) q[0];
sx q[0];
rz(3.9570424) q[0];
sx q[0];
rz(11.856356) q[0];
rz(2.8275936) q[1];
sx q[1];
rz(-2.2065838) q[1];
sx q[1];
rz(1.3318055) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0667257) q[0];
sx q[0];
rz(-0.35860379) q[0];
sx q[0];
rz(-0.68645607) q[0];
x q[1];
rz(1.8237782) q[2];
sx q[2];
rz(-1.2377146) q[2];
sx q[2];
rz(-2.2673504) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.097447473) q[1];
sx q[1];
rz(-2.4821979) q[1];
sx q[1];
rz(1.161518) q[1];
x q[2];
rz(-1.0053357) q[3];
sx q[3];
rz(-2.4255803) q[3];
sx q[3];
rz(2.211192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4189202) q[2];
sx q[2];
rz(-1.2166497) q[2];
sx q[2];
rz(0.81532064) q[2];
rz(-2.3319862) q[3];
sx q[3];
rz(-1.4987192) q[3];
sx q[3];
rz(-1.0838375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0631436) q[0];
sx q[0];
rz(-2.0922631) q[0];
sx q[0];
rz(-0.11058841) q[0];
rz(-0.94509205) q[1];
sx q[1];
rz(-0.4393622) q[1];
sx q[1];
rz(1.5792712) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5131683) q[0];
sx q[0];
rz(-1.1187176) q[0];
sx q[0];
rz(-2.3356209) q[0];
x q[1];
rz(1.1223824) q[2];
sx q[2];
rz(-1.9539991) q[2];
sx q[2];
rz(2.2926083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89797276) q[1];
sx q[1];
rz(-2.0531516) q[1];
sx q[1];
rz(-1.4528758) q[1];
rz(-2.3195295) q[3];
sx q[3];
rz(-1.2044319) q[3];
sx q[3];
rz(-1.0613393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9940146) q[2];
sx q[2];
rz(-2.9958604) q[2];
sx q[2];
rz(-0.4501403) q[2];
rz(2.0077997) q[3];
sx q[3];
rz(-1.4530051) q[3];
sx q[3];
rz(-2.1639737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56871539) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(0.28133389) q[0];
rz(1.7549134) q[1];
sx q[1];
rz(-1.4298871) q[1];
sx q[1];
rz(-0.6001572) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.068014) q[0];
sx q[0];
rz(-1.2805443) q[0];
sx q[0];
rz(1.6192379) q[0];
rz(1.6515031) q[2];
sx q[2];
rz(-2.5332402) q[2];
sx q[2];
rz(2.0474912) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0827209) q[1];
sx q[1];
rz(-1.8959672) q[1];
sx q[1];
rz(2.5356958) q[1];
x q[2];
rz(2.8801877) q[3];
sx q[3];
rz(-2.5525744) q[3];
sx q[3];
rz(-0.95688577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0448138) q[2];
sx q[2];
rz(-1.1324977) q[2];
sx q[2];
rz(1.0463932) q[2];
rz(0.3977631) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(3.0090295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1716877) q[0];
sx q[0];
rz(-0.02709087) q[0];
sx q[0];
rz(2.0890253) q[0];
rz(1.196208) q[1];
sx q[1];
rz(-1.9320678) q[1];
sx q[1];
rz(-0.41950163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013166817) q[0];
sx q[0];
rz(-1.3072398) q[0];
sx q[0];
rz(2.5530784) q[0];
x q[1];
rz(0.93650903) q[2];
sx q[2];
rz(-1.7542931) q[2];
sx q[2];
rz(1.2681792) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4070523) q[1];
sx q[1];
rz(-2.190935) q[1];
sx q[1];
rz(-2.6728515) q[1];
x q[2];
rz(2.6330248) q[3];
sx q[3];
rz(-1.3150104) q[3];
sx q[3];
rz(0.05318197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9516051) q[2];
sx q[2];
rz(-1.1185458) q[2];
sx q[2];
rz(-1.0901964) q[2];
rz(0.53127855) q[3];
sx q[3];
rz(-2.4254906) q[3];
sx q[3];
rz(1.5268911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4215609) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(1.3979727) q[0];
rz(0.13294237) q[1];
sx q[1];
rz(-1.839919) q[1];
sx q[1];
rz(-0.001210777) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5842297) q[0];
sx q[0];
rz(-1.5546636) q[0];
sx q[0];
rz(0.38059142) q[0];
x q[1];
rz(2.9820061) q[2];
sx q[2];
rz(-1.8529466) q[2];
sx q[2];
rz(-0.42195502) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9494755) q[1];
sx q[1];
rz(-1.2671952) q[1];
sx q[1];
rz(-0.66302256) q[1];
x q[2];
rz(-0.026592606) q[3];
sx q[3];
rz(-1.7204434) q[3];
sx q[3];
rz(-0.95434556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2289537) q[2];
sx q[2];
rz(-1.8043844) q[2];
sx q[2];
rz(-0.21052989) q[2];
rz(1.0362961) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(1.2150631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1933111) q[0];
sx q[0];
rz(-1.9183777) q[0];
sx q[0];
rz(-0.40570983) q[0];
rz(1.7871208) q[1];
sx q[1];
rz(-2.3190119) q[1];
sx q[1];
rz(-0.92232651) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8278481) q[0];
sx q[0];
rz(-1.2822312) q[0];
sx q[0];
rz(-0.11103481) q[0];
rz(0.90409235) q[2];
sx q[2];
rz(-1.4897926) q[2];
sx q[2];
rz(-1.9503649) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74444425) q[1];
sx q[1];
rz(-1.6849298) q[1];
sx q[1];
rz(0.95790095) q[1];
x q[2];
rz(-2.5821463) q[3];
sx q[3];
rz(-1.2370951) q[3];
sx q[3];
rz(0.73742407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7102082) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(-0.40194884) q[2];
rz(-1.8534144) q[3];
sx q[3];
rz(-2.2737019) q[3];
sx q[3];
rz(-1.2791546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5227018) q[0];
sx q[0];
rz(-1.1033449) q[0];
sx q[0];
rz(-0.064362137) q[0];
rz(2.3035658) q[1];
sx q[1];
rz(-0.82632724) q[1];
sx q[1];
rz(1.8109969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143599) q[0];
sx q[0];
rz(-1.6961251) q[0];
sx q[0];
rz(-2.5663239) q[0];
rz(-pi) q[1];
rz(-1.072791) q[2];
sx q[2];
rz(-2.118022) q[2];
sx q[2];
rz(0.71820532) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76112642) q[1];
sx q[1];
rz(-2.5674324) q[1];
sx q[1];
rz(0.4131743) q[1];
rz(-pi) q[2];
rz(-2.9693611) q[3];
sx q[3];
rz(-2.1759998) q[3];
sx q[3];
rz(-1.3394522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4055206) q[2];
sx q[2];
rz(-1.6330999) q[2];
sx q[2];
rz(0.68354496) q[2];
rz(-2.4077967) q[3];
sx q[3];
rz(-1.7283231) q[3];
sx q[3];
rz(2.2939513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7965294) q[0];
sx q[0];
rz(-2.0579484) q[0];
sx q[0];
rz(-1.9684568) q[0];
rz(2.7660811) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(-1.0367905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79085892) q[0];
sx q[0];
rz(-3.1111801) q[0];
sx q[0];
rz(2.4615672) q[0];
rz(-pi) q[1];
rz(-2.1027742) q[2];
sx q[2];
rz(-2.8237282) q[2];
sx q[2];
rz(1.3587855) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1634961) q[1];
sx q[1];
rz(-1.4304203) q[1];
sx q[1];
rz(2.3771493) q[1];
rz(-pi) q[2];
rz(2.310318) q[3];
sx q[3];
rz(-0.6704125) q[3];
sx q[3];
rz(1.1408653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57572395) q[2];
sx q[2];
rz(-2.2278991) q[2];
sx q[2];
rz(2.4110528) q[2];
rz(-2.9324487) q[3];
sx q[3];
rz(-0.52339619) q[3];
sx q[3];
rz(1.2701579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4850979) q[0];
sx q[0];
rz(-0.42335835) q[0];
sx q[0];
rz(-3.094161) q[0];
rz(0.221953) q[1];
sx q[1];
rz(-1.0739948) q[1];
sx q[1];
rz(-0.91112959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61924261) q[0];
sx q[0];
rz(-1.6043092) q[0];
sx q[0];
rz(-1.5438118) q[0];
x q[1];
rz(-0.014737562) q[2];
sx q[2];
rz(-0.92916766) q[2];
sx q[2];
rz(2.3126938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2173913) q[1];
sx q[1];
rz(-1.5581521) q[1];
sx q[1];
rz(2.7740757) q[1];
rz(-pi) q[2];
rz(-2.9373475) q[3];
sx q[3];
rz(-1.0139272) q[3];
sx q[3];
rz(-2.000252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.86965108) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(2.9730049) q[2];
rz(-0.21909675) q[3];
sx q[3];
rz(-2.2897661) q[3];
sx q[3];
rz(-1.3271837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.99437) q[0];
sx q[0];
rz(-0.34039012) q[0];
sx q[0];
rz(0.45743531) q[0];
rz(-1.9153197) q[1];
sx q[1];
rz(-0.33718449) q[1];
sx q[1];
rz(-1.7785243) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4922614) q[0];
sx q[0];
rz(-0.96423414) q[0];
sx q[0];
rz(-1.4973634) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8149557) q[2];
sx q[2];
rz(-0.95521046) q[2];
sx q[2];
rz(-0.89060874) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.05039617) q[1];
sx q[1];
rz(-1.1938057) q[1];
sx q[1];
rz(1.3891549) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0567597) q[3];
sx q[3];
rz(-2.2515577) q[3];
sx q[3];
rz(0.57455243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3739796) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(-2.6895831) q[2];
rz(2.7376145) q[3];
sx q[3];
rz(-1.890506) q[3];
sx q[3];
rz(2.0154791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44611888) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(2.6192464) q[1];
sx q[1];
rz(-0.53032395) q[1];
sx q[1];
rz(1.2293336) q[1];
rz(-2.4056388) q[2];
sx q[2];
rz(-0.37307586) q[2];
sx q[2];
rz(1.697344) q[2];
rz(1.7200851) q[3];
sx q[3];
rz(-1.8315157) q[3];
sx q[3];
rz(-0.45890148) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
