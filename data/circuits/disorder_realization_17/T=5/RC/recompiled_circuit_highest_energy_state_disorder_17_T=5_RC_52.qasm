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
rz(-1.0626592) q[0];
sx q[0];
rz(-0.81544977) q[0];
sx q[0];
rz(-2.4315779) q[0];
rz(2.8275936) q[1];
sx q[1];
rz(-2.2065838) q[1];
sx q[1];
rz(-1.8097872) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074866991) q[0];
sx q[0];
rz(-2.7829889) q[0];
sx q[0];
rz(2.4551366) q[0];
x q[1];
rz(-2.5153036) q[2];
sx q[2];
rz(-2.7262027) q[2];
sx q[2];
rz(-2.9363869) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.097447473) q[1];
sx q[1];
rz(-2.4821979) q[1];
sx q[1];
rz(-1.9800746) q[1];
rz(0.43621938) q[3];
sx q[3];
rz(-2.1582104) q[3];
sx q[3];
rz(0.23107108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4189202) q[2];
sx q[2];
rz(-1.2166497) q[2];
sx q[2];
rz(-2.326272) q[2];
rz(-0.80960649) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(-1.0838375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.078449) q[0];
sx q[0];
rz(-1.0493295) q[0];
sx q[0];
rz(-0.11058841) q[0];
rz(-0.94509205) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(1.5623215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5131683) q[0];
sx q[0];
rz(-1.1187176) q[0];
sx q[0];
rz(2.3356209) q[0];
rz(-2.3198691) q[2];
sx q[2];
rz(-0.58124776) q[2];
sx q[2];
rz(1.7591214) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.64792127) q[1];
sx q[1];
rz(-2.646138) q[1];
sx q[1];
rz(-2.9205771) q[1];
rz(-pi) q[2];
rz(-0.82206313) q[3];
sx q[3];
rz(-1.9371607) q[3];
sx q[3];
rz(2.0802534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.14757806) q[2];
sx q[2];
rz(-0.14573228) q[2];
sx q[2];
rz(-0.4501403) q[2];
rz(-1.133793) q[3];
sx q[3];
rz(-1.4530051) q[3];
sx q[3];
rz(-2.1639737) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56871539) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(-0.28133389) q[0];
rz(-1.7549134) q[1];
sx q[1];
rz(-1.7117056) q[1];
sx q[1];
rz(2.5414355) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2413797) q[0];
sx q[0];
rz(-0.29415392) q[0];
sx q[0];
rz(-2.9808729) q[0];
rz(-2.1776206) q[2];
sx q[2];
rz(-1.5247048) q[2];
sx q[2];
rz(0.54296903) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8474947) q[1];
sx q[1];
rz(-2.1408242) q[1];
sx q[1];
rz(-1.1815726) q[1];
rz(-pi) q[2];
rz(0.5731606) q[3];
sx q[3];
rz(-1.7148682) q[3];
sx q[3];
rz(-0.39502963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0448138) q[2];
sx q[2];
rz(-1.1324977) q[2];
sx q[2];
rz(1.0463932) q[2];
rz(-0.3977631) q[3];
sx q[3];
rz(-2.4989276) q[3];
sx q[3];
rz(-0.1325632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9699049) q[0];
sx q[0];
rz(-3.1145018) q[0];
sx q[0];
rz(2.0890253) q[0];
rz(-1.9453847) q[1];
sx q[1];
rz(-1.9320678) q[1];
sx q[1];
rz(2.722091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2118111) q[0];
sx q[0];
rz(-0.63840862) q[0];
sx q[0];
rz(-0.45244658) q[0];
rz(0.2264465) q[2];
sx q[2];
rz(-0.94880494) q[2];
sx q[2];
rz(0.43606191) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1238112) q[1];
sx q[1];
rz(-0.7582802) q[1];
sx q[1];
rz(-1.0067382) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49286263) q[3];
sx q[3];
rz(-0.56418428) q[3];
sx q[3];
rz(1.9436556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9516051) q[2];
sx q[2];
rz(-2.0230468) q[2];
sx q[2];
rz(-2.0513963) q[2];
rz(-0.53127855) q[3];
sx q[3];
rz(-0.71610206) q[3];
sx q[3];
rz(-1.6147015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7200318) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(-1.3979727) q[0];
rz(0.13294237) q[1];
sx q[1];
rz(-1.3016737) q[1];
sx q[1];
rz(-3.1403819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5842297) q[0];
sx q[0];
rz(-1.5869291) q[0];
sx q[0];
rz(2.7610012) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2852077) q[2];
sx q[2];
rz(-1.7240217) q[2];
sx q[2];
rz(1.1040579) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.397314) q[1];
sx q[1];
rz(-0.71956735) q[1];
sx q[1];
rz(0.47081635) q[1];
rz(-3.115) q[3];
sx q[3];
rz(-1.7204434) q[3];
sx q[3];
rz(-2.1872471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2289537) q[2];
sx q[2];
rz(-1.8043844) q[2];
sx q[2];
rz(-0.21052989) q[2];
rz(-1.0362961) q[3];
sx q[3];
rz(-0.784289) q[3];
sx q[3];
rz(-1.9265296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1933111) q[0];
sx q[0];
rz(-1.9183777) q[0];
sx q[0];
rz(-0.40570983) q[0];
rz(-1.7871208) q[1];
sx q[1];
rz(-0.82258075) q[1];
sx q[1];
rz(-0.92232651) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3137445) q[0];
sx q[0];
rz(-1.8593615) q[0];
sx q[0];
rz(0.11103481) q[0];
rz(3.0386557) q[2];
sx q[2];
rz(-2.2349226) q[2];
sx q[2];
rz(2.8256106) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74444425) q[1];
sx q[1];
rz(-1.4566629) q[1];
sx q[1];
rz(2.1836917) q[1];
rz(0.55944632) q[3];
sx q[3];
rz(-1.9044975) q[3];
sx q[3];
rz(2.4041686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4313844) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(2.7396438) q[2];
rz(-1.2881783) q[3];
sx q[3];
rz(-2.2737019) q[3];
sx q[3];
rz(1.2791546) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5227018) q[0];
sx q[0];
rz(-1.1033449) q[0];
sx q[0];
rz(-3.0772305) q[0];
rz(-0.83802682) q[1];
sx q[1];
rz(-0.82632724) q[1];
sx q[1];
rz(-1.3305957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2272328) q[0];
sx q[0];
rz(-1.6961251) q[0];
sx q[0];
rz(2.5663239) q[0];
rz(-pi) q[1];
rz(2.4766972) q[2];
sx q[2];
rz(-0.722363) q[2];
sx q[2];
rz(-1.5255873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3804662) q[1];
sx q[1];
rz(-0.57416026) q[1];
sx q[1];
rz(0.4131743) q[1];
x q[2];
rz(-1.8136386) q[3];
sx q[3];
rz(-2.5153219) q[3];
sx q[3];
rz(2.0988718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4055206) q[2];
sx q[2];
rz(-1.5084927) q[2];
sx q[2];
rz(-0.68354496) q[2];
rz(-2.4077967) q[3];
sx q[3];
rz(-1.4132696) q[3];
sx q[3];
rz(0.84764135) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34506327) q[0];
sx q[0];
rz(-1.0836443) q[0];
sx q[0];
rz(1.1731359) q[0];
rz(-2.7660811) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(1.0367905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0309854) q[0];
sx q[0];
rz(-1.5471503) q[0];
sx q[0];
rz(-1.551669) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0388184) q[2];
sx q[2];
rz(-2.8237282) q[2];
sx q[2];
rz(1.7828072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4040096) q[1];
sx q[1];
rz(-0.77464691) q[1];
sx q[1];
rz(-2.9402016) q[1];
x q[2];
rz(0.83127465) q[3];
sx q[3];
rz(-2.4711802) q[3];
sx q[3];
rz(1.1408653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.57572395) q[2];
sx q[2];
rz(-2.2278991) q[2];
sx q[2];
rz(2.4110528) q[2];
rz(0.20914397) q[3];
sx q[3];
rz(-0.52339619) q[3];
sx q[3];
rz(1.2701579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65649477) q[0];
sx q[0];
rz(-0.42335835) q[0];
sx q[0];
rz(-0.04743162) q[0];
rz(0.221953) q[1];
sx q[1];
rz(-1.0739948) q[1];
sx q[1];
rz(2.2304631) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95064934) q[0];
sx q[0];
rz(-1.5977657) q[0];
sx q[0];
rz(0.033525056) q[0];
x q[1];
rz(-1.59052) q[2];
sx q[2];
rz(-2.4998186) q[2];
sx q[2];
rz(2.3373147) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.641727) q[1];
sx q[1];
rz(-1.9382825) q[1];
sx q[1];
rz(-1.5843452) q[1];
rz(-pi) q[2];
rz(-1.0044596) q[3];
sx q[3];
rz(-1.7438423) q[3];
sx q[3];
rz(-0.53849788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86965108) q[2];
sx q[2];
rz(-2.7615774) q[2];
sx q[2];
rz(-2.9730049) q[2];
rz(-0.21909675) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(-1.8144089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.99437) q[0];
sx q[0];
rz(-0.34039012) q[0];
sx q[0];
rz(2.6841573) q[0];
rz(-1.9153197) q[1];
sx q[1];
rz(-0.33718449) q[1];
sx q[1];
rz(1.3630684) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4922614) q[0];
sx q[0];
rz(-2.1773585) q[0];
sx q[0];
rz(-1.4973634) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62984939) q[2];
sx q[2];
rz(-1.3721264) q[2];
sx q[2];
rz(0.82306403) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41238989) q[1];
sx q[1];
rz(-0.41657428) q[1];
sx q[1];
rz(-2.7135486) q[1];
x q[2];
rz(0.084832974) q[3];
sx q[3];
rz(-0.89003497) q[3];
sx q[3];
rz(2.5670402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7676131) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(2.6895831) q[2];
rz(0.40397817) q[3];
sx q[3];
rz(-1.2510866) q[3];
sx q[3];
rz(-1.1261136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44611888) q[0];
sx q[0];
rz(-1.6071381) q[0];
sx q[0];
rz(0.18679609) q[0];
rz(-2.6192464) q[1];
sx q[1];
rz(-2.6112687) q[1];
sx q[1];
rz(-1.912259) q[1];
rz(-0.28235565) q[2];
sx q[2];
rz(-1.8179802) q[2];
sx q[2];
rz(2.5674934) q[2];
rz(1.4215076) q[3];
sx q[3];
rz(-1.310077) q[3];
sx q[3];
rz(2.6826912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
