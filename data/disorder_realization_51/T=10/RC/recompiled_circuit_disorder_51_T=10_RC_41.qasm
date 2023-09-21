OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5470619) q[0];
sx q[0];
rz(4.2630258) q[0];
sx q[0];
rz(6.1017258) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(1.0531309) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7979413) q[0];
sx q[0];
rz(-2.6290253) q[0];
sx q[0];
rz(-1.1287862) q[0];
rz(-pi) q[1];
rz(2.5187831) q[2];
sx q[2];
rz(-1.8138759) q[2];
sx q[2];
rz(-2.958975) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92703544) q[1];
sx q[1];
rz(-1.7662449) q[1];
sx q[1];
rz(-2.3656225) q[1];
x q[2];
rz(-2.9795322) q[3];
sx q[3];
rz(-1.9063623) q[3];
sx q[3];
rz(0.1048564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9782605) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(2.1286428) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(-0.071454123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.098009) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(-2.5464771) q[0];
rz(1.0455421) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(1.5140623) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058115) q[0];
sx q[0];
rz(-2.0683056) q[0];
sx q[0];
rz(2.6742629) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5397649) q[2];
sx q[2];
rz(-2.4040262) q[2];
sx q[2];
rz(2.8013128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.086386911) q[1];
sx q[1];
rz(-1.7654164) q[1];
sx q[1];
rz(-0.17533949) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4684832) q[3];
sx q[3];
rz(-1.521109) q[3];
sx q[3];
rz(1.697584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65511584) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(-2.3454323) q[2];
rz(-0.97186175) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(0.17175737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765091) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(0.31072101) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(-1.4583189) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38760936) q[0];
sx q[0];
rz(-1.4235272) q[0];
sx q[0];
rz(3.0940042) q[0];
x q[1];
rz(-2.4086191) q[2];
sx q[2];
rz(-0.79967116) q[2];
sx q[2];
rz(3.1250931) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.872936) q[1];
sx q[1];
rz(-1.5169414) q[1];
sx q[1];
rz(2.7320646) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5536669) q[3];
sx q[3];
rz(-1.2560085) q[3];
sx q[3];
rz(-2.0103612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98214275) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(1.3726161) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(-2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963592) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(2.0171719) q[0];
rz(1.2202948) q[1];
sx q[1];
rz(-1.9492457) q[1];
sx q[1];
rz(1.6569998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9254018) q[0];
sx q[0];
rz(-1.6214217) q[0];
sx q[0];
rz(-3.1111858) q[0];
rz(-pi) q[1];
rz(-0.19214432) q[2];
sx q[2];
rz(-1.0612744) q[2];
sx q[2];
rz(2.0083049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93814072) q[1];
sx q[1];
rz(-1.005299) q[1];
sx q[1];
rz(2.6184665) q[1];
rz(-pi) q[2];
rz(-2.7709511) q[3];
sx q[3];
rz(-2.3834043) q[3];
sx q[3];
rz(1.7565808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1241887) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(2.1253288) q[2];
rz(-1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50773412) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(-2.741709) q[0];
rz(-1.9790861) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-2.9679325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0608873) q[0];
sx q[0];
rz(-1.4982491) q[0];
sx q[0];
rz(1.5425496) q[0];
rz(-pi) q[1];
rz(-0.47541754) q[2];
sx q[2];
rz(-0.43736514) q[2];
sx q[2];
rz(1.6944483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9160737) q[1];
sx q[1];
rz(-1.542463) q[1];
sx q[1];
rz(-1.6760992) q[1];
rz(-pi) q[2];
rz(-1.3624304) q[3];
sx q[3];
rz(-1.2456129) q[3];
sx q[3];
rz(1.9735826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48012039) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(2.373467) q[2];
rz(0.85401946) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(2.1693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.1059234) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.4703898) q[0];
rz(0.51180965) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(-1.9981729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4758159) q[0];
sx q[0];
rz(-1.4205298) q[0];
sx q[0];
rz(1.4861602) q[0];
x q[1];
rz(-2.5383699) q[2];
sx q[2];
rz(-1.8458741) q[2];
sx q[2];
rz(3.0657257) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.032635078) q[1];
sx q[1];
rz(-2.7320478) q[1];
sx q[1];
rz(1.1213379) q[1];
rz(0.30541909) q[3];
sx q[3];
rz(-1.050204) q[3];
sx q[3];
rz(0.49314317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4986971) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(1.6112304) q[2];
rz(-1.6879843) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11944184) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(1.7013593) q[0];
rz(-2.4123736) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(1.1332606) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92111174) q[0];
sx q[0];
rz(-2.8163914) q[0];
sx q[0];
rz(0.13334206) q[0];
rz(-pi) q[1];
rz(-0.015958162) q[2];
sx q[2];
rz(-2.4292813) q[2];
sx q[2];
rz(2.6238837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3164191) q[1];
sx q[1];
rz(-2.2146041) q[1];
sx q[1];
rz(-0.34367798) q[1];
rz(1.728119) q[3];
sx q[3];
rz(-1.8295349) q[3];
sx q[3];
rz(1.3820005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7363654) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(2.6531632) q[2];
rz(1.3119665) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(0.64129889) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(3.0704165) q[0];
rz(0.03216234) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(1.2088998) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29584822) q[0];
sx q[0];
rz(-0.83202067) q[0];
sx q[0];
rz(-0.31006281) q[0];
x q[1];
rz(0.9135984) q[2];
sx q[2];
rz(-0.29619869) q[2];
sx q[2];
rz(1.2072472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3370918) q[1];
sx q[1];
rz(-0.13131222) q[1];
sx q[1];
rz(1.4485703) q[1];
rz(-pi) q[2];
rz(2.5197221) q[3];
sx q[3];
rz(-1.8727881) q[3];
sx q[3];
rz(2.5627476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20748392) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(0.87336826) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(-0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97380012) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(2.8299676) q[0];
rz(2.3198126) q[1];
sx q[1];
rz(-2.5506134) q[1];
sx q[1];
rz(1.6315546) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25699297) q[0];
sx q[0];
rz(-1.25367) q[0];
sx q[0];
rz(1.7811799) q[0];
rz(1.0636343) q[2];
sx q[2];
rz(-0.38804752) q[2];
sx q[2];
rz(-1.0840814) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.65731243) q[1];
sx q[1];
rz(-1.3687951) q[1];
sx q[1];
rz(-2.5217418) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4771677) q[3];
sx q[3];
rz(-0.29237745) q[3];
sx q[3];
rz(-1.2053306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2733549) q[2];
sx q[2];
rz(-2.6028825) q[2];
sx q[2];
rz(-2.3256425) q[2];
rz(2.6319035) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(-1.9036487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2508535) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(0.62581217) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(0.65840107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9068245) q[0];
sx q[0];
rz(-1.0848197) q[0];
sx q[0];
rz(-0.66396873) q[0];
x q[1];
rz(-2.7294331) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(1.621643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60724466) q[1];
sx q[1];
rz(-1.736728) q[1];
sx q[1];
rz(-2.0545309) q[1];
rz(-1.9202616) q[3];
sx q[3];
rz(-1.06171) q[3];
sx q[3];
rz(1.6125319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4663503) q[2];
sx q[2];
rz(-0.8224951) q[2];
sx q[2];
rz(2.8038483) q[2];
rz(-1.0234458) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(-1.9256928) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52453775) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(1.9032003) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(1.3255618) q[2];
sx q[2];
rz(-0.64074466) q[2];
sx q[2];
rz(-1.7111039) q[2];
rz(0.57237207) q[3];
sx q[3];
rz(-1.5170245) q[3];
sx q[3];
rz(-1.7046884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];