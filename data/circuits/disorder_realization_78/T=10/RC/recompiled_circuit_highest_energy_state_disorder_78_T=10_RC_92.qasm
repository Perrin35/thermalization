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
rz(3.8905191) q[0];
sx q[0];
rz(10.818304) q[0];
rz(-2.3387609) q[1];
sx q[1];
rz(-2.3454911) q[1];
sx q[1];
rz(2.610745) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323468) q[0];
sx q[0];
rz(-2.0589925) q[0];
sx q[0];
rz(-0.40803473) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5535001) q[2];
sx q[2];
rz(-1.7375542) q[2];
sx q[2];
rz(1.1823428) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2474607) q[1];
sx q[1];
rz(-2.1016205) q[1];
sx q[1];
rz(-0.08695072) q[1];
x q[2];
rz(-2.0898607) q[3];
sx q[3];
rz(-1.4822019) q[3];
sx q[3];
rz(-0.54852329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.406245) q[2];
sx q[2];
rz(-0.9333868) q[2];
sx q[2];
rz(-2.1437342) q[2];
rz(-2.2570299) q[3];
sx q[3];
rz(-1.9622784) q[3];
sx q[3];
rz(2.4295889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74580055) q[0];
sx q[0];
rz(-0.53676787) q[0];
sx q[0];
rz(1.0774379) q[0];
rz(-2.5224345) q[1];
sx q[1];
rz(-1.2666603) q[1];
sx q[1];
rz(-1.9177297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9400459) q[0];
sx q[0];
rz(-0.34149656) q[0];
sx q[0];
rz(-1.0251371) q[0];
rz(1.6009839) q[2];
sx q[2];
rz(-1.5689625) q[2];
sx q[2];
rz(-0.48170127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.465816) q[1];
sx q[1];
rz(-2.4022033) q[1];
sx q[1];
rz(2.5104816) q[1];
rz(-pi) q[2];
rz(2.3540865) q[3];
sx q[3];
rz(-1.0117784) q[3];
sx q[3];
rz(-2.1632568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.18664843) q[2];
sx q[2];
rz(-2.2961605) q[2];
sx q[2];
rz(2.1050982) q[2];
rz(-1.1896108) q[3];
sx q[3];
rz(-1.2396953) q[3];
sx q[3];
rz(-3.0724683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67967296) q[0];
sx q[0];
rz(-0.24149495) q[0];
sx q[0];
rz(-1.4659721) q[0];
rz(-1.8849323) q[1];
sx q[1];
rz(-2.4938221) q[1];
sx q[1];
rz(-2.5648436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58698049) q[0];
sx q[0];
rz(-0.99444992) q[0];
sx q[0];
rz(-0.89926855) q[0];
rz(-pi) q[1];
rz(2.6827963) q[2];
sx q[2];
rz(-1.5374628) q[2];
sx q[2];
rz(-1.3622487) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4119775) q[1];
sx q[1];
rz(-1.5190647) q[1];
sx q[1];
rz(-1.9843141) q[1];
rz(-3.0528487) q[3];
sx q[3];
rz(-1.4074433) q[3];
sx q[3];
rz(0.08069144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74167788) q[2];
sx q[2];
rz(-0.33881131) q[2];
sx q[2];
rz(-2.3728288) q[2];
rz(0.1712884) q[3];
sx q[3];
rz(-1.7007217) q[3];
sx q[3];
rz(-2.7870074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6866368) q[0];
sx q[0];
rz(-0.88382116) q[0];
sx q[0];
rz(0.53632847) q[0];
rz(-2.3388011) q[1];
sx q[1];
rz(-1.8575467) q[1];
sx q[1];
rz(0.098043052) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2778846) q[0];
sx q[0];
rz(-1.4914762) q[0];
sx q[0];
rz(1.2427928) q[0];
x q[1];
rz(-2.1613902) q[2];
sx q[2];
rz(-0.88104311) q[2];
sx q[2];
rz(0.99926567) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6373972) q[1];
sx q[1];
rz(-1.7736083) q[1];
sx q[1];
rz(-2.0921973) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9997659) q[3];
sx q[3];
rz(-1.9460856) q[3];
sx q[3];
rz(0.39794014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1596277) q[2];
sx q[2];
rz(-1.030913) q[2];
sx q[2];
rz(0.31943303) q[2];
rz(-2.5751513) q[3];
sx q[3];
rz(-2.3947377) q[3];
sx q[3];
rz(1.1613065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(2.7403858) q[0];
sx q[0];
rz(-1.3196608) q[0];
sx q[0];
rz(0.55289406) q[0];
rz(-0.7792019) q[1];
sx q[1];
rz(-0.82256493) q[1];
sx q[1];
rz(2.353277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9181128) q[0];
sx q[0];
rz(-1.3105375) q[0];
sx q[0];
rz(1.4994367) q[0];
rz(-pi) q[1];
rz(1.5946024) q[2];
sx q[2];
rz(-1.362065) q[2];
sx q[2];
rz(-3.0301651) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32994575) q[1];
sx q[1];
rz(-1.1378985) q[1];
sx q[1];
rz(-1.2654773) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8184999) q[3];
sx q[3];
rz(-1.8616397) q[3];
sx q[3];
rz(1.3210306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.39764443) q[2];
sx q[2];
rz(-2.0136191) q[2];
sx q[2];
rz(0.20565847) q[2];
rz(1.7913943) q[3];
sx q[3];
rz(-2.4009027) q[3];
sx q[3];
rz(0.010644309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3288997) q[0];
sx q[0];
rz(-0.44726547) q[0];
sx q[0];
rz(1.6269667) q[0];
rz(2.336592) q[1];
sx q[1];
rz(-1.4470419) q[1];
sx q[1];
rz(0.31594333) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.508873) q[0];
sx q[0];
rz(-2.609786) q[0];
sx q[0];
rz(0.72152941) q[0];
rz(-pi) q[1];
rz(-2.6348389) q[2];
sx q[2];
rz(-0.67595081) q[2];
sx q[2];
rz(1.4434222) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.915357) q[1];
sx q[1];
rz(-1.8457883) q[1];
sx q[1];
rz(-0.26696856) q[1];
rz(-0.44734939) q[3];
sx q[3];
rz(-2.7998689) q[3];
sx q[3];
rz(2.0096092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.463795) q[2];
sx q[2];
rz(-0.69910502) q[2];
sx q[2];
rz(0.15764906) q[2];
rz(0.53358233) q[3];
sx q[3];
rz(-1.6505417) q[3];
sx q[3];
rz(-1.6238448) q[3];
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
rz(-2.5428298) q[0];
sx q[0];
rz(-0.49648008) q[0];
sx q[0];
rz(-0.98094034) q[0];
rz(0.44278231) q[1];
sx q[1];
rz(-1.1863703) q[1];
sx q[1];
rz(2.669899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9424635) q[0];
sx q[0];
rz(-0.53265306) q[0];
sx q[0];
rz(-0.96766242) q[0];
rz(-pi) q[1];
rz(0.54092714) q[2];
sx q[2];
rz(-2.8188883) q[2];
sx q[2];
rz(0.91992119) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.973418) q[1];
sx q[1];
rz(-1.7061632) q[1];
sx q[1];
rz(-1.5568887) q[1];
x q[2];
rz(2.8011462) q[3];
sx q[3];
rz(-1.3967112) q[3];
sx q[3];
rz(-0.3917087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9337351) q[2];
sx q[2];
rz(-2.2722878) q[2];
sx q[2];
rz(0.11824879) q[2];
rz(-2.5737428) q[3];
sx q[3];
rz(-1.7852781) q[3];
sx q[3];
rz(0.80287272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0826913) q[0];
sx q[0];
rz(-1.7401798) q[0];
sx q[0];
rz(0.68728224) q[0];
rz(1.8742689) q[1];
sx q[1];
rz(-1.9272389) q[1];
sx q[1];
rz(2.5158688) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61266292) q[0];
sx q[0];
rz(-1.5497396) q[0];
sx q[0];
rz(3.1369169) q[0];
rz(0.93333107) q[2];
sx q[2];
rz(-1.0984761) q[2];
sx q[2];
rz(2.460091) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81258147) q[1];
sx q[1];
rz(-0.99615232) q[1];
sx q[1];
rz(-0.01477764) q[1];
x q[2];
rz(-0.3097624) q[3];
sx q[3];
rz(-1.0132257) q[3];
sx q[3];
rz(-0.46771177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9132793) q[2];
sx q[2];
rz(-1.8931171) q[2];
sx q[2];
rz(1.0949562) q[2];
rz(2.6878808) q[3];
sx q[3];
rz(-2.3504421) q[3];
sx q[3];
rz(2.9554844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.5472645) q[0];
sx q[0];
rz(-0.493395) q[0];
sx q[0];
rz(-0.067597978) q[0];
rz(-1.0293845) q[1];
sx q[1];
rz(-1.2815579) q[1];
sx q[1];
rz(3.0060815) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96582039) q[0];
sx q[0];
rz(-1.0403353) q[0];
sx q[0];
rz(0.6483174) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9642815) q[2];
sx q[2];
rz(-2.3925892) q[2];
sx q[2];
rz(0.85402358) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.20457224) q[1];
sx q[1];
rz(-2.1453152) q[1];
sx q[1];
rz(-1.9078698) q[1];
rz(-pi) q[2];
rz(-2.8587946) q[3];
sx q[3];
rz(-1.334189) q[3];
sx q[3];
rz(0.11219926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4542666) q[2];
sx q[2];
rz(-2.8886075) q[2];
sx q[2];
rz(-0.48736408) q[2];
rz(2.7519915) q[3];
sx q[3];
rz(-2.1890958) q[3];
sx q[3];
rz(0.065936955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(2.6479263) q[0];
sx q[0];
rz(-1.0856029) q[0];
sx q[0];
rz(-1.3847463) q[0];
rz(2.0866277) q[1];
sx q[1];
rz(-0.87433785) q[1];
sx q[1];
rz(0.25996444) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11851234) q[0];
sx q[0];
rz(-1.8268088) q[0];
sx q[0];
rz(0.89566083) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4394664) q[2];
sx q[2];
rz(-1.6703147) q[2];
sx q[2];
rz(-0.00015774568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.205906) q[1];
sx q[1];
rz(-2.5126713) q[1];
sx q[1];
rz(-1.511093) q[1];
x q[2];
rz(1.7694951) q[3];
sx q[3];
rz(-1.3532146) q[3];
sx q[3];
rz(-2.8189903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42183033) q[2];
sx q[2];
rz(-0.43785849) q[2];
sx q[2];
rz(-1.3502632) q[2];
rz(-1.0849902) q[3];
sx q[3];
rz(-1.2144054) q[3];
sx q[3];
rz(2.0124281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0706901) q[0];
sx q[0];
rz(-2.4546843) q[0];
sx q[0];
rz(-0.51921459) q[0];
rz(1.4821953) q[1];
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
