OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.06935057) q[0];
sx q[0];
rz(-2.0877593) q[0];
sx q[0];
rz(-1.2487489) q[0];
rz(1.9034003) q[1];
sx q[1];
rz(-0.44013953) q[1];
sx q[1];
rz(1.2425437) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1093501) q[0];
sx q[0];
rz(-2.1324131) q[0];
sx q[0];
rz(-2.0294782) q[0];
x q[1];
rz(-0.094775188) q[2];
sx q[2];
rz(-2.3412626) q[2];
sx q[2];
rz(0.86343414) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93852167) q[1];
sx q[1];
rz(-1.2483828) q[1];
sx q[1];
rz(-1.5702269) q[1];
rz(2.5322024) q[3];
sx q[3];
rz(-1.8577788) q[3];
sx q[3];
rz(3.0521986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.093546346) q[2];
sx q[2];
rz(-0.28154937) q[2];
sx q[2];
rz(1.9264889) q[2];
rz(-1.9563227) q[3];
sx q[3];
rz(-0.90648854) q[3];
sx q[3];
rz(-1.5426481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95618653) q[0];
sx q[0];
rz(-0.39919272) q[0];
sx q[0];
rz(1.4584374) q[0];
rz(-0.1419119) q[1];
sx q[1];
rz(-1.2071573) q[1];
sx q[1];
rz(-1.9292319) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4282517) q[0];
sx q[0];
rz(-1.9104092) q[0];
sx q[0];
rz(2.2560512) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59213068) q[2];
sx q[2];
rz(-1.5685602) q[2];
sx q[2];
rz(2.1606902) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72313655) q[1];
sx q[1];
rz(-2.2424556) q[1];
sx q[1];
rz(-3.0211012) q[1];
x q[2];
rz(-1.2207915) q[3];
sx q[3];
rz(-0.98549609) q[3];
sx q[3];
rz(0.096297527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.045018) q[2];
sx q[2];
rz(-2.9587032) q[2];
sx q[2];
rz(0.96192399) q[2];
rz(2.9436881) q[3];
sx q[3];
rz(-1.3062545) q[3];
sx q[3];
rz(1.4977247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69979954) q[0];
sx q[0];
rz(-0.72023359) q[0];
sx q[0];
rz(1.066347) q[0];
rz(2.9329246) q[1];
sx q[1];
rz(-0.51869789) q[1];
sx q[1];
rz(2.4934703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84635669) q[0];
sx q[0];
rz(-1.6235678) q[0];
sx q[0];
rz(-2.1480302) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0096613) q[2];
sx q[2];
rz(-2.6474617) q[2];
sx q[2];
rz(0.021683387) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3274182) q[1];
sx q[1];
rz(-1.2035554) q[1];
sx q[1];
rz(0.7805853) q[1];
rz(0.1005248) q[3];
sx q[3];
rz(-0.70954743) q[3];
sx q[3];
rz(-1.2840934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2842747) q[2];
sx q[2];
rz(-1.6782574) q[2];
sx q[2];
rz(2.5939482) q[2];
rz(1.0037496) q[3];
sx q[3];
rz(-2.818483) q[3];
sx q[3];
rz(-2.9799262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169287) q[0];
sx q[0];
rz(-1.1268317) q[0];
sx q[0];
rz(2.774985) q[0];
rz(1.8560575) q[1];
sx q[1];
rz(-1.5204241) q[1];
sx q[1];
rz(-2.3191648) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7912589) q[0];
sx q[0];
rz(-2.9295577) q[0];
sx q[0];
rz(-0.37104443) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9738178) q[2];
sx q[2];
rz(-1.3047403) q[2];
sx q[2];
rz(-0.38420907) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.3318204) q[1];
sx q[1];
rz(-1.0882411) q[1];
sx q[1];
rz(0.81565522) q[1];
rz(-pi) q[2];
rz(0.27536686) q[3];
sx q[3];
rz(-1.5726701) q[3];
sx q[3];
rz(3.0089859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6377247) q[2];
sx q[2];
rz(-0.34054264) q[2];
sx q[2];
rz(1.5857504) q[2];
rz(1.2268892) q[3];
sx q[3];
rz(-2.5033958) q[3];
sx q[3];
rz(-1.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.9521088) q[0];
sx q[0];
rz(-2.0229078) q[0];
sx q[0];
rz(2.191191) q[0];
rz(0.19418007) q[1];
sx q[1];
rz(-1.2545398) q[1];
sx q[1];
rz(0.28087428) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9683198) q[0];
sx q[0];
rz(-1.313198) q[0];
sx q[0];
rz(2.2031242) q[0];
rz(-1.5581649) q[2];
sx q[2];
rz(-1.4213741) q[2];
sx q[2];
rz(1.0557501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5944526) q[1];
sx q[1];
rz(-3.0277589) q[1];
sx q[1];
rz(-0.74233858) q[1];
x q[2];
rz(-2.0858795) q[3];
sx q[3];
rz(-0.56929811) q[3];
sx q[3];
rz(1.4610987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6516271) q[2];
sx q[2];
rz(-1.6941035) q[2];
sx q[2];
rz(0.70072407) q[2];
rz(2.0476332) q[3];
sx q[3];
rz(-2.14812) q[3];
sx q[3];
rz(2.3195364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.1489498) q[0];
sx q[0];
rz(-1.0473017) q[0];
sx q[0];
rz(2.6348422) q[0];
rz(-2.0172334) q[1];
sx q[1];
rz(-1.9729112) q[1];
sx q[1];
rz(-0.36875025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4294037) q[0];
sx q[0];
rz(-2.2656239) q[0];
sx q[0];
rz(0.24526986) q[0];
rz(-pi) q[1];
rz(2.314744) q[2];
sx q[2];
rz(-1.294181) q[2];
sx q[2];
rz(-2.6630304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1146093) q[1];
sx q[1];
rz(-0.31216808) q[1];
sx q[1];
rz(1.8725558) q[1];
rz(-pi) q[2];
rz(-2.3193377) q[3];
sx q[3];
rz(-2.5516627) q[3];
sx q[3];
rz(2.024789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2345978) q[2];
sx q[2];
rz(-2.5556421) q[2];
sx q[2];
rz(0.21729812) q[2];
rz(1.4761188) q[3];
sx q[3];
rz(-2.3264591) q[3];
sx q[3];
rz(0.34172094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9076964) q[0];
sx q[0];
rz(-2.8192769) q[0];
sx q[0];
rz(0.79793683) q[0];
rz(-1.7012874) q[1];
sx q[1];
rz(-2.1013575) q[1];
sx q[1];
rz(0.61680102) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26825702) q[0];
sx q[0];
rz(-2.5028879) q[0];
sx q[0];
rz(3.0144342) q[0];
rz(0.5011214) q[2];
sx q[2];
rz(-0.62868147) q[2];
sx q[2];
rz(-1.4528265) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4176158) q[1];
sx q[1];
rz(-0.25088746) q[1];
sx q[1];
rz(-2.1881359) q[1];
rz(-2.2374898) q[3];
sx q[3];
rz(-0.98424235) q[3];
sx q[3];
rz(-1.3400638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0133609) q[2];
sx q[2];
rz(-1.1527454) q[2];
sx q[2];
rz(-2.4398003) q[2];
rz(-2.2026786) q[3];
sx q[3];
rz(-2.2434668) q[3];
sx q[3];
rz(-1.2452589) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.668648) q[0];
sx q[0];
rz(-2.9845147) q[0];
sx q[0];
rz(0.12338403) q[0];
rz(-0.75048796) q[1];
sx q[1];
rz(-1.4742943) q[1];
sx q[1];
rz(-2.1231245) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.930542) q[0];
sx q[0];
rz(-1.3500542) q[0];
sx q[0];
rz(3.0342191) q[0];
rz(-pi) q[1];
rz(-0.46393259) q[2];
sx q[2];
rz(-2.0407157) q[2];
sx q[2];
rz(1.7352833) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97216304) q[1];
sx q[1];
rz(-1.0854682) q[1];
sx q[1];
rz(1.5281488) q[1];
rz(0.090224548) q[3];
sx q[3];
rz(-2.1975747) q[3];
sx q[3];
rz(2.9240328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5156775) q[2];
sx q[2];
rz(-1.4695784) q[2];
sx q[2];
rz(0.53543004) q[2];
rz(1.2299906) q[3];
sx q[3];
rz(-2.1401236) q[3];
sx q[3];
rz(3.1297704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54315058) q[0];
sx q[0];
rz(-2.4813528) q[0];
sx q[0];
rz(-1.862233) q[0];
rz(1.2641501) q[1];
sx q[1];
rz(-1.8708355) q[1];
sx q[1];
rz(0.15170161) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21148602) q[0];
sx q[0];
rz(-1.9686342) q[0];
sx q[0];
rz(-1.5087288) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13387605) q[2];
sx q[2];
rz(-1.5508442) q[2];
sx q[2];
rz(-2.6559115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7599185) q[1];
sx q[1];
rz(-1.2592788) q[1];
sx q[1];
rz(2.1550234) q[1];
rz(-pi) q[2];
rz(-2.5429498) q[3];
sx q[3];
rz(-0.81889048) q[3];
sx q[3];
rz(0.18695565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5558527) q[2];
sx q[2];
rz(-0.93583217) q[2];
sx q[2];
rz(1.9110511) q[2];
rz(-1.3452283) q[3];
sx q[3];
rz(-1.3627005) q[3];
sx q[3];
rz(1.8432553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1718488) q[0];
sx q[0];
rz(-1.3405223) q[0];
sx q[0];
rz(0.44542435) q[0];
rz(-2.716966) q[1];
sx q[1];
rz(-1.5716962) q[1];
sx q[1];
rz(0.62579036) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9741083) q[0];
sx q[0];
rz(-1.0667598) q[0];
sx q[0];
rz(0.91668769) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0883415) q[2];
sx q[2];
rz(-1.4212593) q[2];
sx q[2];
rz(0.47433269) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6546331) q[1];
sx q[1];
rz(-0.59096293) q[1];
sx q[1];
rz(-1.1019568) q[1];
x q[2];
rz(1.5324235) q[3];
sx q[3];
rz(-1.4202446) q[3];
sx q[3];
rz(2.2398219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17973913) q[2];
sx q[2];
rz(-2.1135795) q[2];
sx q[2];
rz(-1.6801768) q[2];
rz(3.0840868) q[3];
sx q[3];
rz(-1.4534566) q[3];
sx q[3];
rz(-2.1632975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4089324) q[0];
sx q[0];
rz(-1.6608149) q[0];
sx q[0];
rz(0.59857359) q[0];
rz(0.21845017) q[1];
sx q[1];
rz(-1.598806) q[1];
sx q[1];
rz(1.7576408) q[1];
rz(1.4870395) q[2];
sx q[2];
rz(-1.0261921) q[2];
sx q[2];
rz(1.5941317) q[2];
rz(0.045513734) q[3];
sx q[3];
rz(-0.78452605) q[3];
sx q[3];
rz(1.8922643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
