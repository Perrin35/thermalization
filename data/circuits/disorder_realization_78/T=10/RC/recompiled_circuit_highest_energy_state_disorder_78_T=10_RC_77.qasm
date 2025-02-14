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
rz(-1.7595093) q[0];
sx q[0];
rz(-2.3926662) q[0];
sx q[0];
rz(1.393526) q[0];
rz(-2.3387609) q[1];
sx q[1];
rz(-2.3454911) q[1];
sx q[1];
rz(2.610745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323468) q[0];
sx q[0];
rz(-1.0826001) q[0];
sx q[0];
rz(-2.7335579) q[0];
rz(-pi) q[1];
rz(-0.16678236) q[2];
sx q[2];
rz(-1.5878526) q[2];
sx q[2];
rz(-0.38558233) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4208274) q[1];
sx q[1];
rz(-1.6457575) q[1];
sx q[1];
rz(1.0383181) q[1];
rz(-3.0396456) q[3];
sx q[3];
rz(-1.0539712) q[3];
sx q[3];
rz(-2.068813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.406245) q[2];
sx q[2];
rz(-2.2082059) q[2];
sx q[2];
rz(0.99785844) q[2];
rz(0.88456279) q[3];
sx q[3];
rz(-1.9622784) q[3];
sx q[3];
rz(2.4295889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3957921) q[0];
sx q[0];
rz(-0.53676787) q[0];
sx q[0];
rz(2.0641548) q[0];
rz(-0.61915818) q[1];
sx q[1];
rz(-1.2666603) q[1];
sx q[1];
rz(-1.223863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.773945) q[0];
sx q[0];
rz(-1.2804693) q[0];
sx q[0];
rz(2.9591857) q[0];
rz(-pi) q[1];
rz(1.5406088) q[2];
sx q[2];
rz(-1.5689625) q[2];
sx q[2];
rz(0.48170127) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.465816) q[1];
sx q[1];
rz(-0.73938939) q[1];
sx q[1];
rz(0.63111102) q[1];
rz(-pi) q[2];
rz(2.2961462) q[3];
sx q[3];
rz(-0.92638141) q[3];
sx q[3];
rz(3.0385142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18664843) q[2];
sx q[2];
rz(-0.84543219) q[2];
sx q[2];
rz(2.1050982) q[2];
rz(1.9519818) q[3];
sx q[3];
rz(-1.2396953) q[3];
sx q[3];
rz(0.069124393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67967296) q[0];
sx q[0];
rz(-2.9000977) q[0];
sx q[0];
rz(-1.6756206) q[0];
rz(-1.8849323) q[1];
sx q[1];
rz(-2.4938221) q[1];
sx q[1];
rz(0.57674903) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.758848) q[0];
sx q[0];
rz(-0.85475105) q[0];
sx q[0];
rz(-2.3780253) q[0];
rz(-pi) q[1];
rz(-1.6079712) q[2];
sx q[2];
rz(-2.0293183) q[2];
sx q[2];
rz(-2.9165845) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4119775) q[1];
sx q[1];
rz(-1.5190647) q[1];
sx q[1];
rz(1.1572786) q[1];
x q[2];
rz(1.7347832) q[3];
sx q[3];
rz(-1.4832368) q[3];
sx q[3];
rz(1.5045741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3999148) q[2];
sx q[2];
rz(-2.8027813) q[2];
sx q[2];
rz(2.3728288) q[2];
rz(2.9703043) q[3];
sx q[3];
rz(-1.4408709) q[3];
sx q[3];
rz(0.35458529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.6866368) q[0];
sx q[0];
rz(-2.2577715) q[0];
sx q[0];
rz(-0.53632847) q[0];
rz(2.3388011) q[1];
sx q[1];
rz(-1.8575467) q[1];
sx q[1];
rz(3.0435496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8756384) q[0];
sx q[0];
rz(-1.2438626) q[0];
sx q[0];
rz(0.083766706) q[0];
x q[1];
rz(-2.5478159) q[2];
sx q[2];
rz(-2.2662518) q[2];
sx q[2];
rz(-1.8107294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72934276) q[1];
sx q[1];
rz(-2.5855484) q[1];
sx q[1];
rz(-1.9623318) q[1];
x q[2];
rz(-0.81252247) q[3];
sx q[3];
rz(-0.56216633) q[3];
sx q[3];
rz(2.6443329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98196491) q[2];
sx q[2];
rz(-1.030913) q[2];
sx q[2];
rz(2.8221596) q[2];
rz(0.56644136) q[3];
sx q[3];
rz(-2.3947377) q[3];
sx q[3];
rz(-1.9802861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40120688) q[0];
sx q[0];
rz(-1.3196608) q[0];
sx q[0];
rz(-0.55289406) q[0];
rz(-2.3623908) q[1];
sx q[1];
rz(-2.3190277) q[1];
sx q[1];
rz(-0.78831569) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32892431) q[0];
sx q[0];
rz(-1.5018437) q[0];
sx q[0];
rz(-0.2608932) q[0];
rz(-pi) q[1];
rz(-3.0296828) q[2];
sx q[2];
rz(-2.9315278) q[2];
sx q[2];
rz(3.1386167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3723046) q[1];
sx q[1];
rz(-1.8471698) q[1];
sx q[1];
rz(0.45120542) q[1];
rz(-1.3230927) q[3];
sx q[3];
rz(-1.2799529) q[3];
sx q[3];
rz(-1.820562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.39764443) q[2];
sx q[2];
rz(-1.1279736) q[2];
sx q[2];
rz(-2.9359342) q[2];
rz(-1.3501984) q[3];
sx q[3];
rz(-0.74068991) q[3];
sx q[3];
rz(-0.010644309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.3288997) q[0];
sx q[0];
rz(-2.6943272) q[0];
sx q[0];
rz(1.6269667) q[0];
rz(2.336592) q[1];
sx q[1];
rz(-1.4470419) q[1];
sx q[1];
rz(0.31594333) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8522911) q[0];
sx q[0];
rz(-1.9123494) q[0];
sx q[0];
rz(0.41595398) q[0];
rz(-pi) q[1];
rz(-2.6348389) q[2];
sx q[2];
rz(-2.4656418) q[2];
sx q[2];
rz(1.6981704) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8711612) q[1];
sx q[1];
rz(-1.3140895) q[1];
sx q[1];
rz(1.8553599) q[1];
rz(-2.6942433) q[3];
sx q[3];
rz(-2.7998689) q[3];
sx q[3];
rz(1.1319834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67779764) q[2];
sx q[2];
rz(-2.4424876) q[2];
sx q[2];
rz(-0.15764906) q[2];
rz(2.6080103) q[3];
sx q[3];
rz(-1.6505417) q[3];
sx q[3];
rz(-1.5177479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59876281) q[0];
sx q[0];
rz(-2.6451126) q[0];
sx q[0];
rz(-2.1606523) q[0];
rz(2.6988103) q[1];
sx q[1];
rz(-1.1863703) q[1];
sx q[1];
rz(-2.669899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169051) q[0];
sx q[0];
rz(-1.1393095) q[0];
sx q[0];
rz(-0.3226852) q[0];
rz(-pi) q[1];
rz(-1.4002789) q[2];
sx q[2];
rz(-1.2954752) q[2];
sx q[2];
rz(-1.6570651) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.737094) q[1];
sx q[1];
rz(-1.5570159) q[1];
sx q[1];
rz(3.0062129) q[1];
x q[2];
rz(2.8011462) q[3];
sx q[3];
rz(-1.3967112) q[3];
sx q[3];
rz(2.749884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20785759) q[2];
sx q[2];
rz(-2.2722878) q[2];
sx q[2];
rz(-3.0233439) q[2];
rz(0.56784981) q[3];
sx q[3];
rz(-1.3563145) q[3];
sx q[3];
rz(-0.80287272) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0826913) q[0];
sx q[0];
rz(-1.7401798) q[0];
sx q[0];
rz(-2.4543104) q[0];
rz(-1.2673238) q[1];
sx q[1];
rz(-1.2143538) q[1];
sx q[1];
rz(0.6257239) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1835577) q[0];
sx q[0];
rz(-1.5661217) q[0];
sx q[0];
rz(1.5497394) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56628801) q[2];
sx q[2];
rz(-2.1293961) q[2];
sx q[2];
rz(-1.2142912) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3*pi/4) q[1];
sx q[1];
rz(-0.57481261) q[1];
sx q[1];
rz(-1.5936127) q[1];
rz(-pi) q[2];
rz(0.3097624) q[3];
sx q[3];
rz(-2.1283669) q[3];
sx q[3];
rz(2.6738809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9132793) q[2];
sx q[2];
rz(-1.2484756) q[2];
sx q[2];
rz(2.0466364) q[2];
rz(2.6878808) q[3];
sx q[3];
rz(-0.79115051) q[3];
sx q[3];
rz(-2.9554844) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5472645) q[0];
sx q[0];
rz(-0.493395) q[0];
sx q[0];
rz(-3.0739947) q[0];
rz(-1.0293845) q[1];
sx q[1];
rz(-1.8600347) q[1];
sx q[1];
rz(0.13551113) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9477959) q[0];
sx q[0];
rz(-2.3289177) q[0];
sx q[0];
rz(-0.77085797) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4082528) q[2];
sx q[2];
rz(-0.83627273) q[2];
sx q[2];
rz(-0.61406174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1780336) q[1];
sx q[1];
rz(-1.8521223) q[1];
sx q[1];
rz(0.60124361) q[1];
rz(-pi) q[2];
rz(1.816808) q[3];
sx q[3];
rz(-1.8455077) q[3];
sx q[3];
rz(1.6149855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4542666) q[2];
sx q[2];
rz(-0.25298515) q[2];
sx q[2];
rz(2.6542286) q[2];
rz(-0.38960114) q[3];
sx q[3];
rz(-0.95249683) q[3];
sx q[3];
rz(3.0756557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49366632) q[0];
sx q[0];
rz(-1.0856029) q[0];
sx q[0];
rz(1.7568463) q[0];
rz(-1.0549649) q[1];
sx q[1];
rz(-2.2672548) q[1];
sx q[1];
rz(-0.25996444) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1459634) q[0];
sx q[0];
rz(-0.71486799) q[0];
sx q[0];
rz(1.174182) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2222338) q[2];
sx q[2];
rz(-0.16460379) q[2];
sx q[2];
rz(2.2155264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1321226) q[1];
sx q[1];
rz(-0.94317052) q[1];
sx q[1];
rz(-0.043379003) q[1];
rz(-pi) q[2];
rz(-1.7694951) q[3];
sx q[3];
rz(-1.7883781) q[3];
sx q[3];
rz(-2.8189903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42183033) q[2];
sx q[2];
rz(-0.43785849) q[2];
sx q[2];
rz(1.3502632) q[2];
rz(2.0566025) q[3];
sx q[3];
rz(-1.9271873) q[3];
sx q[3];
rz(-2.0124281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0709025) q[0];
sx q[0];
rz(-0.68690837) q[0];
sx q[0];
rz(2.6223781) q[0];
rz(-1.6593973) q[1];
sx q[1];
rz(-1.2980325) q[1];
sx q[1];
rz(1.4549805) q[1];
rz(0.71375511) q[2];
sx q[2];
rz(-2.426385) q[2];
sx q[2];
rz(-0.89520988) q[2];
rz(-2.0706035) q[3];
sx q[3];
rz(-1.8355814) q[3];
sx q[3];
rz(1.7818835) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
