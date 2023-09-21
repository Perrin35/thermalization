OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274347) q[0];
sx q[0];
rz(-0.43570575) q[0];
sx q[0];
rz(-2.2154007) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(-2.7690673) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6719138) q[0];
sx q[0];
rz(-2.1436474) q[0];
sx q[0];
rz(1.8125305) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63106545) q[2];
sx q[2];
rz(-1.6562781) q[2];
sx q[2];
rz(0.98110547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26973154) q[1];
sx q[1];
rz(-1.6850123) q[1];
sx q[1];
rz(1.3807339) q[1];
rz(1.7872693) q[3];
sx q[3];
rz(-1.6811922) q[3];
sx q[3];
rz(-2.20761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42840502) q[2];
sx q[2];
rz(-1.5822072) q[2];
sx q[2];
rz(-2.2170846) q[2];
rz(-1.472578) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(-1.6424461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532042) q[0];
sx q[0];
rz(-3.0792455) q[0];
sx q[0];
rz(-1.5361319) q[0];
rz(2.9470782) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(-3.0867192) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0957735) q[0];
sx q[0];
rz(-1.4574483) q[0];
sx q[0];
rz(3.0485831) q[0];
rz(-2.1129677) q[2];
sx q[2];
rz(-2.1193868) q[2];
sx q[2];
rz(0.86663914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.730719) q[1];
sx q[1];
rz(-1.5667856) q[1];
sx q[1];
rz(1.9133168) q[1];
rz(2.4684286) q[3];
sx q[3];
rz(-1.737397) q[3];
sx q[3];
rz(-1.8002321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.43869552) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(1.4820209) q[2];
rz(-2.1510018) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.3668485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7822587) q[0];
sx q[0];
rz(-2.6222836) q[0];
sx q[0];
rz(-2.7666132) q[0];
rz(-0.24770501) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(-1.8050271) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96391812) q[0];
sx q[0];
rz(-1.5128711) q[0];
sx q[0];
rz(1.792568) q[0];
rz(-pi) q[1];
rz(0.61632421) q[2];
sx q[2];
rz(-0.68683544) q[2];
sx q[2];
rz(-0.64885215) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3414351) q[1];
sx q[1];
rz(-2.2812124) q[1];
sx q[1];
rz(0.40409778) q[1];
rz(-2.7258354) q[3];
sx q[3];
rz(-2.3299179) q[3];
sx q[3];
rz(-0.31879253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1027362) q[2];
sx q[2];
rz(-0.83013022) q[2];
sx q[2];
rz(-1.0245727) q[2];
rz(-1.864795) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17722002) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(-0.033551034) q[0];
rz(1.230348) q[1];
sx q[1];
rz(-2.3760445) q[1];
sx q[1];
rz(-1.7376602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0516292) q[0];
sx q[0];
rz(-1.7647867) q[0];
sx q[0];
rz(2.7882663) q[0];
rz(-pi) q[1];
rz(-0.83547445) q[2];
sx q[2];
rz(-2.2611141) q[2];
sx q[2];
rz(0.26963216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85642203) q[1];
sx q[1];
rz(-0.31177786) q[1];
sx q[1];
rz(-1.1986033) q[1];
x q[2];
rz(-2.2622044) q[3];
sx q[3];
rz(-1.5048358) q[3];
sx q[3];
rz(2.0397759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4924865) q[2];
sx q[2];
rz(-2.2968569) q[2];
sx q[2];
rz(2.9283294) q[2];
rz(3.1212741) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(2.7385353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3605109) q[0];
sx q[0];
rz(-2.0251944) q[0];
sx q[0];
rz(-1.0158585) q[0];
rz(1.6332743) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(3.1075081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40092295) q[0];
sx q[0];
rz(-1.4105083) q[0];
sx q[0];
rz(2.7435061) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.836852) q[2];
sx q[2];
rz(-1.012946) q[2];
sx q[2];
rz(-1.7709874) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7923813) q[1];
sx q[1];
rz(-1.5776909) q[1];
sx q[1];
rz(-1.2864119) q[1];
rz(0.50877737) q[3];
sx q[3];
rz(-2.1366589) q[3];
sx q[3];
rz(-2.1600427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4013227) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(1.1711228) q[2];
rz(1.8088388) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(-0.12935054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4532582) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(-0.43584287) q[0];
rz(2.0360937) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(-1.2058535) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1036557) q[0];
sx q[0];
rz(-0.78563848) q[0];
sx q[0];
rz(1.908179) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4756104) q[2];
sx q[2];
rz(-1.5548445) q[2];
sx q[2];
rz(2.8841281) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6916794) q[1];
sx q[1];
rz(-2.0813585) q[1];
sx q[1];
rz(1.540586) q[1];
x q[2];
rz(1.1711575) q[3];
sx q[3];
rz(-1.540278) q[3];
sx q[3];
rz(2.7680824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2445406) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(0.051606027) q[2];
rz(-2.7339325) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5200941) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(-0.67333418) q[0];
rz(-2.3576221) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(0.53692445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4250701) q[0];
sx q[0];
rz(-1.636583) q[0];
sx q[0];
rz(1.487088) q[0];
rz(-pi) q[1];
rz(1.0579254) q[2];
sx q[2];
rz(-0.29209902) q[2];
sx q[2];
rz(-0.7823173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6785477) q[1];
sx q[1];
rz(-1.4765413) q[1];
sx q[1];
rz(2.3386392) q[1];
rz(-pi) q[2];
rz(-0.22646871) q[3];
sx q[3];
rz(-0.60301757) q[3];
sx q[3];
rz(2.9035567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15988222) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(-2.9505777) q[2];
rz(-0.28132176) q[3];
sx q[3];
rz(-1.9502935) q[3];
sx q[3];
rz(-0.34240001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1087588) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(-2.7539745) q[0];
rz(3.0265813) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(0.33755916) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25305504) q[0];
sx q[0];
rz(-2.6153784) q[0];
sx q[0];
rz(2.7387268) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0701596) q[2];
sx q[2];
rz(-2.2087503) q[2];
sx q[2];
rz(-2.273794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.90606373) q[1];
sx q[1];
rz(-1.4396832) q[1];
sx q[1];
rz(2.1436585) q[1];
x q[2];
rz(-1.6345331) q[3];
sx q[3];
rz(-0.81634854) q[3];
sx q[3];
rz(0.07894978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9541624) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(-1.2672651) q[2];
rz(2.3665442) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(-2.3601941) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18701126) q[0];
sx q[0];
rz(-1.0382074) q[0];
sx q[0];
rz(-0.22098456) q[0];
rz(-0.9221319) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(2.1386713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.636165) q[0];
sx q[0];
rz(-1.822585) q[0];
sx q[0];
rz(-1.7136784) q[0];
x q[1];
rz(0.50311627) q[2];
sx q[2];
rz(-2.0965577) q[2];
sx q[2];
rz(-0.8612649) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3993966) q[1];
sx q[1];
rz(-0.70397111) q[1];
sx q[1];
rz(2.3854371) q[1];
rz(-pi) q[2];
rz(0.88463155) q[3];
sx q[3];
rz(-1.7836708) q[3];
sx q[3];
rz(0.5205982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6691436) q[2];
sx q[2];
rz(-0.74206918) q[2];
sx q[2];
rz(2.9830902) q[2];
rz(0.45378271) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(0.84428549) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52699387) q[0];
sx q[0];
rz(-0.90899962) q[0];
sx q[0];
rz(0.19113834) q[0];
rz(0.29516164) q[1];
sx q[1];
rz(-0.8894397) q[1];
sx q[1];
rz(0.89231649) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44360456) q[0];
sx q[0];
rz(-0.83570489) q[0];
sx q[0];
rz(2.9025335) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92992444) q[2];
sx q[2];
rz(-1.495549) q[2];
sx q[2];
rz(-1.4571232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82980624) q[1];
sx q[1];
rz(-1.808842) q[1];
sx q[1];
rz(-1.5029552) q[1];
rz(-2.7097706) q[3];
sx q[3];
rz(-1.8027657) q[3];
sx q[3];
rz(-2.1524129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3698547) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(2.2820293) q[2];
rz(1.9178948) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(-2.3971476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201465) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(1.6620811) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(-1.4798726) q[2];
sx q[2];
rz(-0.29279136) q[2];
sx q[2];
rz(1.472483) q[2];
rz(-2.0966093) q[3];
sx q[3];
rz(-0.079418728) q[3];
sx q[3];
rz(3.1374745) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
