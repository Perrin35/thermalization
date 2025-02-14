OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.6264853) q[0];
sx q[0];
rz(6.1942979) q[0];
sx q[0];
rz(11.789378) q[0];
rz(0.69500336) q[1];
sx q[1];
rz(-3.0106795) q[1];
sx q[1];
rz(0.34832365) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94960864) q[0];
sx q[0];
rz(-0.79553793) q[0];
sx q[0];
rz(2.863671) q[0];
rz(-2.8856002) q[2];
sx q[2];
rz(-2.8446994) q[2];
sx q[2];
rz(0.96742899) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6435223) q[1];
sx q[1];
rz(-0.76979945) q[1];
sx q[1];
rz(0.78242015) q[1];
x q[2];
rz(0.11669604) q[3];
sx q[3];
rz(-2.7435997) q[3];
sx q[3];
rz(-1.7577049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4942887) q[2];
sx q[2];
rz(-0.56607294) q[2];
sx q[2];
rz(1.0757793) q[2];
rz(-3.0268269) q[3];
sx q[3];
rz(-1.4817295) q[3];
sx q[3];
rz(-2.4858294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7455416) q[0];
sx q[0];
rz(-0.93463722) q[0];
sx q[0];
rz(2.4387687) q[0];
rz(-1.5952236) q[1];
sx q[1];
rz(-2.399235) q[1];
sx q[1];
rz(0.57483086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84687656) q[0];
sx q[0];
rz(-2.623577) q[0];
sx q[0];
rz(2.4344492) q[0];
rz(-pi) q[1];
rz(-1.0444591) q[2];
sx q[2];
rz(-1.0246201) q[2];
sx q[2];
rz(2.1072497) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3177111) q[1];
sx q[1];
rz(-2.1924824) q[1];
sx q[1];
rz(-1.409338) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63268707) q[3];
sx q[3];
rz(-1.1345769) q[3];
sx q[3];
rz(-1.3496646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90291643) q[2];
sx q[2];
rz(-1.6409589) q[2];
sx q[2];
rz(2.404876) q[2];
rz(-2.369407) q[3];
sx q[3];
rz(-0.18852791) q[3];
sx q[3];
rz(0.34672117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18948983) q[0];
sx q[0];
rz(-1.075407) q[0];
sx q[0];
rz(1.1601675) q[0];
rz(-0.53933764) q[1];
sx q[1];
rz(-2.5181006) q[1];
sx q[1];
rz(0.070405237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6269913) q[0];
sx q[0];
rz(-0.87341269) q[0];
sx q[0];
rz(0.91430362) q[0];
rz(-pi) q[1];
rz(-0.83924509) q[2];
sx q[2];
rz(-0.23942023) q[2];
sx q[2];
rz(-3.1175304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16500948) q[1];
sx q[1];
rz(-1.1900468) q[1];
sx q[1];
rz(-1.6409207) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0892404) q[3];
sx q[3];
rz(-0.97739554) q[3];
sx q[3];
rz(-0.39077047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.77014852) q[2];
sx q[2];
rz(-1.4263209) q[2];
sx q[2];
rz(-2.9798129) q[2];
rz(0.77156228) q[3];
sx q[3];
rz(-3.0050889) q[3];
sx q[3];
rz(-2.0989044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55018392) q[0];
sx q[0];
rz(-0.78224459) q[0];
sx q[0];
rz(2.8915306) q[0];
rz(1.0391883) q[1];
sx q[1];
rz(-0.74165529) q[1];
sx q[1];
rz(2.340462) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88177365) q[0];
sx q[0];
rz(-1.7823813) q[0];
sx q[0];
rz(-1.5863281) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.807937) q[2];
sx q[2];
rz(-1.8225428) q[2];
sx q[2];
rz(2.8636914) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5591079) q[1];
sx q[1];
rz(-1.6402157) q[1];
sx q[1];
rz(0.51377929) q[1];
x q[2];
rz(1.9362373) q[3];
sx q[3];
rz(-1.1574928) q[3];
sx q[3];
rz(-0.45729056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4074126) q[2];
sx q[2];
rz(-1.4449747) q[2];
sx q[2];
rz(0.59047353) q[2];
rz(-2.7967795) q[3];
sx q[3];
rz(-0.36546388) q[3];
sx q[3];
rz(1.4709681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.19584) q[0];
sx q[0];
rz(-1.7029637) q[0];
sx q[0];
rz(-1.6666743) q[0];
rz(1.6189812) q[1];
sx q[1];
rz(-1.5983862) q[1];
sx q[1];
rz(-2.7326857) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58137888) q[0];
sx q[0];
rz(-1.5239851) q[0];
sx q[0];
rz(-1.5792094) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7402503) q[2];
sx q[2];
rz(-2.231488) q[2];
sx q[2];
rz(-1.5398538) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.489456) q[1];
sx q[1];
rz(-2.5083275) q[1];
sx q[1];
rz(-0.79852028) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6650531) q[3];
sx q[3];
rz(-2.2426212) q[3];
sx q[3];
rz(-0.49956027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54395479) q[2];
sx q[2];
rz(-0.13263098) q[2];
sx q[2];
rz(-2.4908716) q[2];
rz(-1.4607653) q[3];
sx q[3];
rz(-1.7284349) q[3];
sx q[3];
rz(-0.4167324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8553218) q[0];
sx q[0];
rz(-1.0690419) q[0];
sx q[0];
rz(-1.8977813) q[0];
rz(2.6142201) q[1];
sx q[1];
rz(-1.6970044) q[1];
sx q[1];
rz(-1.3810371) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6731411) q[0];
sx q[0];
rz(-1.8126789) q[0];
sx q[0];
rz(-1.5833447) q[0];
rz(-pi) q[1];
rz(2.5041298) q[2];
sx q[2];
rz(-0.12345498) q[2];
sx q[2];
rz(0.055094624) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3917249) q[1];
sx q[1];
rz(-2.0473745) q[1];
sx q[1];
rz(0.99443545) q[1];
rz(-2.2562669) q[3];
sx q[3];
rz(-2.4866085) q[3];
sx q[3];
rz(1.8670435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1018299) q[2];
sx q[2];
rz(-1.7888174) q[2];
sx q[2];
rz(1.9977894) q[2];
rz(2.3666798) q[3];
sx q[3];
rz(-1.1760005) q[3];
sx q[3];
rz(-2.9162143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5932805) q[0];
sx q[0];
rz(-2.6052124) q[0];
sx q[0];
rz(-0.29658741) q[0];
rz(-2.2071154) q[1];
sx q[1];
rz(-2.2521033) q[1];
sx q[1];
rz(2.8737822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2598617) q[0];
sx q[0];
rz(-1.7825883) q[0];
sx q[0];
rz(0.21881108) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61106793) q[2];
sx q[2];
rz(-1.7253118) q[2];
sx q[2];
rz(0.25777486) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0632656) q[1];
sx q[1];
rz(-0.73042471) q[1];
sx q[1];
rz(2.6759381) q[1];
rz(-pi) q[2];
rz(0.15800627) q[3];
sx q[3];
rz(-1.4337162) q[3];
sx q[3];
rz(-2.7003845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3725793) q[2];
sx q[2];
rz(-2.09435) q[2];
sx q[2];
rz(0.11865842) q[2];
rz(0.81698051) q[3];
sx q[3];
rz(-1.1840362) q[3];
sx q[3];
rz(-0.55678862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7939222) q[0];
sx q[0];
rz(-0.477328) q[0];
sx q[0];
rz(-2.0678066) q[0];
rz(0.99717623) q[1];
sx q[1];
rz(-1.8938096) q[1];
sx q[1];
rz(1.006385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1365161) q[0];
sx q[0];
rz(-3.0126326) q[0];
sx q[0];
rz(1.6411247) q[0];
x q[1];
rz(-1.2362376) q[2];
sx q[2];
rz(-1.2090956) q[2];
sx q[2];
rz(1.6702056) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.83676618) q[1];
sx q[1];
rz(-1.2427034) q[1];
sx q[1];
rz(-1.7734709) q[1];
rz(-pi) q[2];
rz(-2.6615922) q[3];
sx q[3];
rz(-2.1566803) q[3];
sx q[3];
rz(0.65174499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5239079) q[2];
sx q[2];
rz(-2.0670117) q[2];
sx q[2];
rz(0.96134821) q[2];
rz(-0.62054595) q[3];
sx q[3];
rz(-0.70182645) q[3];
sx q[3];
rz(-2.3294241) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0470225) q[0];
sx q[0];
rz(-0.71806878) q[0];
sx q[0];
rz(-0.72599232) q[0];
rz(-2.831931) q[1];
sx q[1];
rz(-2.1084712) q[1];
sx q[1];
rz(2.9594701) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5512729) q[0];
sx q[0];
rz(-2.3819551) q[0];
sx q[0];
rz(-2.9653427) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51515963) q[2];
sx q[2];
rz(-2.023989) q[2];
sx q[2];
rz(0.31441385) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88090522) q[1];
sx q[1];
rz(-1.9753078) q[1];
sx q[1];
rz(-1.717156) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5043132) q[3];
sx q[3];
rz(-1.2361119) q[3];
sx q[3];
rz(-0.94663564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2187659) q[2];
sx q[2];
rz(-1.4309692) q[2];
sx q[2];
rz(1.1440841) q[2];
rz(-2.4036582) q[3];
sx q[3];
rz(-0.72665557) q[3];
sx q[3];
rz(1.1206333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69796193) q[0];
sx q[0];
rz(-1.8101036) q[0];
sx q[0];
rz(2.4586082) q[0];
rz(-1.5525275) q[1];
sx q[1];
rz(-2.2582159) q[1];
sx q[1];
rz(-0.89673269) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67499298) q[0];
sx q[0];
rz(-1.3764125) q[0];
sx q[0];
rz(0.43793508) q[0];
rz(-pi) q[1];
rz(2.9994869) q[2];
sx q[2];
rz(-2.1143205) q[2];
sx q[2];
rz(2.6212111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7234133) q[1];
sx q[1];
rz(-1.4825943) q[1];
sx q[1];
rz(-3.0107016) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0192782) q[3];
sx q[3];
rz(-2.0297673) q[3];
sx q[3];
rz(0.77536264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5985742) q[2];
sx q[2];
rz(-1.2335346) q[2];
sx q[2];
rz(0.42178556) q[2];
rz(-0.56454033) q[3];
sx q[3];
rz(-0.83344236) q[3];
sx q[3];
rz(-0.88258266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9660139) q[0];
sx q[0];
rz(-1.5400374) q[0];
sx q[0];
rz(-1.2651428) q[0];
rz(1.9059146) q[1];
sx q[1];
rz(-1.0503678) q[1];
sx q[1];
rz(-2.8138524) q[1];
rz(-1.3597455) q[2];
sx q[2];
rz(-1.2086972) q[2];
sx q[2];
rz(0.093314485) q[2];
rz(1.0698224) q[3];
sx q[3];
rz(-1.9961431) q[3];
sx q[3];
rz(2.3186984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
