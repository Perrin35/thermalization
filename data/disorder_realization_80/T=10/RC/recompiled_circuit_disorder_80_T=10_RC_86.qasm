OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39188448) q[0];
sx q[0];
rz(-0.19667974) q[0];
sx q[0];
rz(-1.952202) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(-2.7639311) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1188287) q[0];
sx q[0];
rz(-1.4979935) q[0];
sx q[0];
rz(1.4490436) q[0];
x q[1];
rz(-3.112552) q[2];
sx q[2];
rz(-1.2872868) q[2];
sx q[2];
rz(-3.0344506) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3440173) q[1];
sx q[1];
rz(-2.0098149) q[1];
sx q[1];
rz(-2.3637191) q[1];
rz(-pi) q[2];
rz(-2.7226726) q[3];
sx q[3];
rz(-1.6994611) q[3];
sx q[3];
rz(-0.045623771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.010014023) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(2.4689891) q[2];
rz(0.16942313) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(-1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.910903) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(-0.22856523) q[0];
rz(0.16054343) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(-2.8536318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44293091) q[0];
sx q[0];
rz(-1.7611815) q[0];
sx q[0];
rz(-2.9605244) q[0];
rz(-pi) q[1];
rz(-1.0238043) q[2];
sx q[2];
rz(-0.66811251) q[2];
sx q[2];
rz(-2.748865) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60625633) q[1];
sx q[1];
rz(-1.9793946) q[1];
sx q[1];
rz(-2.4819083) q[1];
x q[2];
rz(3.037022) q[3];
sx q[3];
rz(-0.35211709) q[3];
sx q[3];
rz(0.62192384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59445375) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(1.6681558) q[2];
rz(0.93747059) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32524747) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(0.26741272) q[0];
rz(1.422241) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(-2.1898988) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62943447) q[0];
sx q[0];
rz(-0.21244563) q[0];
sx q[0];
rz(1.1266707) q[0];
x q[1];
rz(2.0864262) q[2];
sx q[2];
rz(-1.7679169) q[2];
sx q[2];
rz(-2.4899763) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0905076) q[1];
sx q[1];
rz(-1.8790434) q[1];
sx q[1];
rz(-2.3329263) q[1];
rz(-0.53593105) q[3];
sx q[3];
rz(-2.1246353) q[3];
sx q[3];
rz(1.2371847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(0.34417957) q[2];
rz(2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13720559) q[0];
sx q[0];
rz(-1.4973649) q[0];
sx q[0];
rz(1.244506) q[0];
rz(0.1291153) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(-0.37277645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9739649) q[0];
sx q[0];
rz(-0.81866696) q[0];
sx q[0];
rz(1.7276093) q[0];
rz(0.77809019) q[2];
sx q[2];
rz(-1.2568682) q[2];
sx q[2];
rz(-2.2997466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.43416926) q[1];
sx q[1];
rz(-1.3278828) q[1];
sx q[1];
rz(-1.0711819) q[1];
rz(-pi) q[2];
rz(2.8126206) q[3];
sx q[3];
rz(-1.1165459) q[3];
sx q[3];
rz(-2.8153552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63561511) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(2.2367031) q[2];
rz(0.44057009) q[3];
sx q[3];
rz(-0.4959271) q[3];
sx q[3];
rz(2.1499965) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47914094) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.2878081) q[0];
rz(1.6632535) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(2.6073661) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9984765) q[0];
sx q[0];
rz(-1.2623708) q[0];
sx q[0];
rz(-2.8240859) q[0];
rz(-2.062837) q[2];
sx q[2];
rz(-1.8740219) q[2];
sx q[2];
rz(0.505503) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5200107) q[1];
sx q[1];
rz(-3.0428805) q[1];
sx q[1];
rz(-2.2732544) q[1];
x q[2];
rz(1.7027431) q[3];
sx q[3];
rz(-2.2852516) q[3];
sx q[3];
rz(0.42657846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6405032) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(0.14349288) q[2];
rz(-1.7701373) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(-2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.7025529) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(2.904536) q[0];
rz(1.9006231) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(-1.3751078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2415337) q[0];
sx q[0];
rz(-0.43474475) q[0];
sx q[0];
rz(-0.66381201) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3926922) q[2];
sx q[2];
rz(-1.6237215) q[2];
sx q[2];
rz(2.5089094) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.314234) q[1];
sx q[1];
rz(-1.5463366) q[1];
sx q[1];
rz(-2.7315061) q[1];
rz(-pi) q[2];
rz(-2.4414805) q[3];
sx q[3];
rz(-0.60986076) q[3];
sx q[3];
rz(-1.8240579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6087626) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(-2.9928845) q[2];
rz(-3.1249629) q[3];
sx q[3];
rz(-2.7745268) q[3];
sx q[3];
rz(-0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6454813) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(-2.2684229) q[0];
rz(-0.9219777) q[1];
sx q[1];
rz(-1.1154122) q[1];
sx q[1];
rz(-3.057664) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5556363) q[0];
sx q[0];
rz(-1.3480098) q[0];
sx q[0];
rz(1.0046093) q[0];
rz(-pi) q[1];
x q[1];
rz(2.601868) q[2];
sx q[2];
rz(-2.735609) q[2];
sx q[2];
rz(-0.27137953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0107207) q[1];
sx q[1];
rz(-0.92492231) q[1];
sx q[1];
rz(-0.801553) q[1];
rz(-pi) q[2];
x q[2];
rz(2.300802) q[3];
sx q[3];
rz(-1.6222686) q[3];
sx q[3];
rz(1.3703025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38368791) q[2];
sx q[2];
rz(-0.44054511) q[2];
sx q[2];
rz(0.54523462) q[2];
rz(-2.7111354) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6298744) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(-1.9301201) q[0];
rz(-2.1684872) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(-2.1957695) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57293939) q[0];
sx q[0];
rz(-1.4927673) q[0];
sx q[0];
rz(0.22139876) q[0];
rz(-pi) q[1];
rz(-2.321645) q[2];
sx q[2];
rz(-2.2689399) q[2];
sx q[2];
rz(-1.1993711) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9981873) q[1];
sx q[1];
rz(-0.37969509) q[1];
sx q[1];
rz(3.0642964) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5787015) q[3];
sx q[3];
rz(-2.2369011) q[3];
sx q[3];
rz(-2.283309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8699845) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(-2.8611709) q[2];
rz(0.60493207) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(0.98208565) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9873001) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(-2.5262685) q[0];
rz(2.212021) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(0.5756793) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79918039) q[0];
sx q[0];
rz(-2.3665782) q[0];
sx q[0];
rz(-0.70223017) q[0];
rz(-1.766878) q[2];
sx q[2];
rz(-1.6819281) q[2];
sx q[2];
rz(-2.9642504) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74734028) q[1];
sx q[1];
rz(-1.5451952) q[1];
sx q[1];
rz(1.4634553) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18746312) q[3];
sx q[3];
rz(-2.1037357) q[3];
sx q[3];
rz(0.77123469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3118887) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(-3.1196307) q[2];
rz(0.17045505) q[3];
sx q[3];
rz(-1.0237834) q[3];
sx q[3];
rz(2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72572529) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(0.6341933) q[0];
rz(-0.028907396) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(2.172487) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4934851) q[0];
sx q[0];
rz(-0.095795184) q[0];
sx q[0];
rz(0.84798725) q[0];
rz(-pi) q[1];
rz(0.74659851) q[2];
sx q[2];
rz(-1.4274297) q[2];
sx q[2];
rz(1.6105086) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2057048) q[1];
sx q[1];
rz(-1.74946) q[1];
sx q[1];
rz(-2.4680069) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5892331) q[3];
sx q[3];
rz(-0.5061572) q[3];
sx q[3];
rz(1.6278704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6897631) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(0.36007145) q[2];
rz(-1.5277956) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(-2.578919) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9344899) q[0];
sx q[0];
rz(-1.5710545) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(3.0974401) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(-1.1883696) q[2];
sx q[2];
rz(-1.2570751) q[2];
sx q[2];
rz(-3.0260069) q[2];
rz(-0.58892693) q[3];
sx q[3];
rz(-2.6116461) q[3];
sx q[3];
rz(-2.6245821) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
