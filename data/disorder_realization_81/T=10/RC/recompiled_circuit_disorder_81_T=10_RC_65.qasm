OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(-1.3614549) q[0];
sx q[0];
rz(1.7629495) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(-2.690697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3987797) q[0];
sx q[0];
rz(-3.0516041) q[0];
sx q[0];
rz(2.9773657) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0499174) q[2];
sx q[2];
rz(-1.6699104) q[2];
sx q[2];
rz(1.5168158) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7887468) q[1];
sx q[1];
rz(-2.5874918) q[1];
sx q[1];
rz(0.43177859) q[1];
rz(-pi) q[2];
rz(-2.2545635) q[3];
sx q[3];
rz(-1.4966045) q[3];
sx q[3];
rz(2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(-2.297304) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(0.26309183) q[0];
rz(0.94353765) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(1.1862322) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0378368) q[0];
sx q[0];
rz(-3.0825966) q[0];
sx q[0];
rz(0.32365139) q[0];
rz(-1.6329174) q[2];
sx q[2];
rz(-1.3717692) q[2];
sx q[2];
rz(1.3028499) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7344208) q[1];
sx q[1];
rz(-0.91033519) q[1];
sx q[1];
rz(2.9829626) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9119592) q[3];
sx q[3];
rz(-2.4335055) q[3];
sx q[3];
rz(-1.1585483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(1.1594695) q[2];
rz(-0.37108478) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(-0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.7611258) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(2.3348715) q[0];
rz(0.21356788) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8504282) q[0];
sx q[0];
rz(-0.69520742) q[0];
sx q[0];
rz(-1.3948963) q[0];
rz(-pi) q[1];
rz(-2.9224612) q[2];
sx q[2];
rz(-2.0543155) q[2];
sx q[2];
rz(0.52106524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.324675) q[1];
sx q[1];
rz(-1.1500689) q[1];
sx q[1];
rz(-2.7540728) q[1];
x q[2];
rz(1.1060171) q[3];
sx q[3];
rz(-1.6764063) q[3];
sx q[3];
rz(1.6184023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.31072581) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(2.2107928) q[2];
rz(0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61313066) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(2.7032734) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(-1.8211676) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5116918) q[0];
sx q[0];
rz(-1.9803847) q[0];
sx q[0];
rz(0.011261777) q[0];
rz(-pi) q[1];
rz(0.40789149) q[2];
sx q[2];
rz(-0.48569277) q[2];
sx q[2];
rz(-0.25024589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8561613) q[1];
sx q[1];
rz(-2.397575) q[1];
sx q[1];
rz(-0.39399778) q[1];
rz(-pi) q[2];
rz(-2.1774125) q[3];
sx q[3];
rz(-2.2607431) q[3];
sx q[3];
rz(2.6026158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1057672) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(0.34238112) q[2];
rz(-2.9648182) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(-1.140973) q[3];
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
rz(2.8300366) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(-2.2763021) q[0];
rz(1.9150437) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(1.8409761) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39869719) q[0];
sx q[0];
rz(-1.5960176) q[0];
sx q[0];
rz(-2.7874649) q[0];
rz(2.956203) q[2];
sx q[2];
rz(-2.7692147) q[2];
sx q[2];
rz(0.34639726) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.091150065) q[1];
sx q[1];
rz(-1.0346518) q[1];
sx q[1];
rz(-1.1863143) q[1];
x q[2];
rz(-2.8793094) q[3];
sx q[3];
rz(-1.4296921) q[3];
sx q[3];
rz(1.9067681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0098003) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(-0.47719964) q[2];
rz(-2.9495083) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.3451097) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(3.1298424) q[0];
rz(-2.5911962) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(1.5884429) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9827305) q[0];
sx q[0];
rz(-1.2356865) q[0];
sx q[0];
rz(1.1984675) q[0];
x q[1];
rz(3.1197238) q[2];
sx q[2];
rz(-1.2049897) q[2];
sx q[2];
rz(2.0621698) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3192056) q[1];
sx q[1];
rz(-2.2777595) q[1];
sx q[1];
rz(1.4591401) q[1];
rz(1.6586967) q[3];
sx q[3];
rz(-2.4368736) q[3];
sx q[3];
rz(0.3375012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6340296) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(1.2825512) q[2];
rz(1.7717308) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(1.1184568) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(2.4643331) q[0];
rz(-2.989785) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(2.1645434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2718186) q[0];
sx q[0];
rz(-0.5973814) q[0];
sx q[0];
rz(-3.0074044) q[0];
x q[1];
rz(3.0095519) q[2];
sx q[2];
rz(-1.5706976) q[2];
sx q[2];
rz(1.4591109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8057101) q[1];
sx q[1];
rz(-1.4937595) q[1];
sx q[1];
rz(-1.9105934) q[1];
rz(-2.7988775) q[3];
sx q[3];
rz(-1.2667155) q[3];
sx q[3];
rz(-2.8976687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7523505) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(2.1288669) q[2];
rz(1.9536473) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0174719) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.2493856) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65803618) q[0];
sx q[0];
rz(-1.5244966) q[0];
sx q[0];
rz(2.9306843) q[0];
x q[1];
rz(-2.0935358) q[2];
sx q[2];
rz(-1.9881696) q[2];
sx q[2];
rz(1.0440895) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82259761) q[1];
sx q[1];
rz(-1.5929475) q[1];
sx q[1];
rz(-1.5652565) q[1];
x q[2];
rz(0.93805712) q[3];
sx q[3];
rz(-1.8055827) q[3];
sx q[3];
rz(-2.7126922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.32968783) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(1.1784941) q[2];
rz(1.4568436) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4998528) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(-1.8485803) q[0];
rz(1.4216084) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(2.5440149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9874728) q[0];
sx q[0];
rz(-3.0010536) q[0];
sx q[0];
rz(1.2921635) q[0];
rz(-pi) q[1];
rz(2.5848128) q[2];
sx q[2];
rz(-2.0282929) q[2];
sx q[2];
rz(0.58820398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24039195) q[1];
sx q[1];
rz(-1.8933834) q[1];
sx q[1];
rz(2.8278973) q[1];
x q[2];
rz(1.4937917) q[3];
sx q[3];
rz(-2.1553851) q[3];
sx q[3];
rz(-2.0921752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(1.9082327) q[2];
rz(-2.2402066) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(-1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937061) q[0];
sx q[0];
rz(-0.77195764) q[0];
sx q[0];
rz(3.1179324) q[0];
rz(-0.95611447) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-2.4694209) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4155054) q[0];
sx q[0];
rz(-0.011517631) q[0];
sx q[0];
rz(-2.0477717) q[0];
rz(-pi) q[1];
rz(2.6300738) q[2];
sx q[2];
rz(-1.5788955) q[2];
sx q[2];
rz(-2.9715003) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2939261) q[1];
sx q[1];
rz(-2.3617509) q[1];
sx q[1];
rz(-2.6508209) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7087595) q[3];
sx q[3];
rz(-1.7107309) q[3];
sx q[3];
rz(-2.7452552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6293634) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(-2.771634) q[2];
rz(1.6379179) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(-1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(-2.4178986) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(2.9818515) q[2];
sx q[2];
rz(-1.9651056) q[2];
sx q[2];
rz(-0.38412487) q[2];
rz(0.95009347) q[3];
sx q[3];
rz(-0.48662574) q[3];
sx q[3];
rz(1.1403198) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];