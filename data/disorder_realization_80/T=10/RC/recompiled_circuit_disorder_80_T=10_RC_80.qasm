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
rz(2.9449129) q[0];
sx q[0];
rz(11.37698) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(2.300188) q[1];
sx q[1];
rz(9.0471164) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.057214) q[0];
sx q[0];
rz(-0.14176653) q[0];
sx q[0];
rz(-2.111582) q[0];
rz(-pi) q[1];
x q[1];
rz(1.471465) q[2];
sx q[2];
rz(-0.28495312) q[2];
sx q[2];
rz(-3.1379267) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.764896) q[1];
sx q[1];
rz(-2.2590859) q[1];
sx q[1];
rz(-0.98801686) q[1];
rz(2.8336485) q[3];
sx q[3];
rz(-2.7044798) q[3];
sx q[3];
rz(1.8969769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1315786) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(0.67260355) q[2];
rz(2.9721695) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23068962) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(0.22856523) q[0];
rz(-0.16054343) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(0.28796089) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8156793) q[0];
sx q[0];
rz(-0.26198146) q[0];
sx q[0];
rz(-2.322305) q[0];
rz(-pi) q[1];
rz(0.38950133) q[2];
sx q[2];
rz(-1.0132388) q[2];
sx q[2];
rz(0.26706375) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6504753) q[1];
sx q[1];
rz(-0.75956356) q[1];
sx q[1];
rz(-2.526545) q[1];
rz(3.037022) q[3];
sx q[3];
rz(-2.7894756) q[3];
sx q[3];
rz(2.5196688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(-1.4734369) q[2];
rz(2.2041221) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(0.37500769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32524747) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(2.8741799) q[0];
rz(-1.7193517) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(2.1898988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5121582) q[0];
sx q[0];
rz(-2.929147) q[0];
sx q[0];
rz(-1.1266707) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9159413) q[2];
sx q[2];
rz(-1.0661085) q[2];
sx q[2];
rz(0.80863189) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82745951) q[1];
sx q[1];
rz(-2.3312807) q[1];
sx q[1];
rz(1.1387226) q[1];
rz(-pi) q[2];
rz(-0.94727256) q[3];
sx q[3];
rz(-1.1215278) q[3];
sx q[3];
rz(0.63637966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(-2.7974131) q[2];
rz(0.88642818) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(2.5755431) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13720559) q[0];
sx q[0];
rz(-1.4973649) q[0];
sx q[0];
rz(1.8970867) q[0];
rz(-0.1291153) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(0.37277645) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6308206) q[0];
sx q[0];
rz(-1.456506) q[0];
sx q[0];
rz(-0.75829102) q[0];
rz(-pi) q[1];
rz(1.9984841) q[2];
sx q[2];
rz(-2.3017985) q[2];
sx q[2];
rz(2.1172303) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4200538) q[1];
sx q[1];
rz(-2.5905847) q[1];
sx q[1];
rz(2.0481471) q[1];
x q[2];
rz(0.32897207) q[3];
sx q[3];
rz(-2.0250468) q[3];
sx q[3];
rz(0.3262375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5059775) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(-2.2367031) q[2];
rz(-2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(-2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6624517) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.2878081) q[0];
rz(1.4783391) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(-0.53422654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9984765) q[0];
sx q[0];
rz(-1.2623708) q[0];
sx q[0];
rz(0.31750676) q[0];
x q[1];
rz(2.1557501) q[2];
sx q[2];
rz(-2.5702229) q[2];
sx q[2];
rz(0.55703288) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6492669) q[1];
sx q[1];
rz(-1.5070792) q[1];
sx q[1];
rz(1.6462412) q[1];
rz(-1.7027431) q[3];
sx q[3];
rz(-0.85634106) q[3];
sx q[3];
rz(-2.7150142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(2.9980998) q[2];
rz(1.3714553) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(-0.21970704) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4390398) q[0];
sx q[0];
rz(-2.2655903) q[0];
sx q[0];
rz(-0.23705661) q[0];
rz(1.2409695) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(-1.3751078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8538044) q[0];
sx q[0];
rz(-1.3082936) q[0];
sx q[0];
rz(-0.35065035) q[0];
rz(0.077652046) q[2];
sx q[2];
rz(-2.3911871) q[2];
sx q[2];
rz(2.1466308) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.314234) q[1];
sx q[1];
rz(-1.5463366) q[1];
sx q[1];
rz(0.41008653) q[1];
rz(-pi) q[2];
rz(0.49075134) q[3];
sx q[3];
rz(-1.1928344) q[3];
sx q[3];
rz(-0.85765391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6087626) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(-0.14870816) q[2];
rz(0.016629774) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(-2.9490411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49611133) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(0.87316978) q[0];
rz(0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(0.08392863) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8221995) q[0];
sx q[0];
rz(-2.5376352) q[0];
sx q[0];
rz(-1.1711867) q[0];
rz(-pi) q[1];
x q[1];
rz(2.601868) q[2];
sx q[2];
rz(-2.735609) q[2];
sx q[2];
rz(-0.27137953) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0107207) q[1];
sx q[1];
rz(-2.2166703) q[1];
sx q[1];
rz(0.801553) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6478959) q[3];
sx q[3];
rz(-2.4101083) q[3];
sx q[3];
rz(-0.25792083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38368791) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(-0.54523462) q[2];
rz(0.43045726) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6298744) q[0];
sx q[0];
rz(-3.1261303) q[0];
sx q[0];
rz(1.9301201) q[0];
rz(2.1684872) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(-0.94582311) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5686533) q[0];
sx q[0];
rz(-1.6488254) q[0];
sx q[0];
rz(2.9201939) q[0];
rz(2.2875167) q[2];
sx q[2];
rz(-2.1207361) q[2];
sx q[2];
rz(-2.9727109) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.91499) q[1];
sx q[1];
rz(-1.1922925) q[1];
sx q[1];
rz(-1.539991) q[1];
x q[2];
rz(-2.4754727) q[3];
sx q[3];
rz(-1.564581) q[3];
sx q[3];
rz(2.4339649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2716081) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(0.28042173) q[2];
rz(0.60493207) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(2.159507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1542926) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(0.61532414) q[0];
rz(-2.212021) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(-0.5756793) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3424123) q[0];
sx q[0];
rz(-2.3665782) q[0];
sx q[0];
rz(0.70223017) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0909537) q[2];
sx q[2];
rz(-2.9165604) q[2];
sx q[2];
rz(-0.88423836) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74734028) q[1];
sx q[1];
rz(-1.5451952) q[1];
sx q[1];
rz(1.4634553) q[1];
rz(-pi) q[2];
rz(0.18746312) q[3];
sx q[3];
rz(-1.0378569) q[3];
sx q[3];
rz(-2.370358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3118887) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(0.021961948) q[2];
rz(0.17045505) q[3];
sx q[3];
rz(-1.0237834) q[3];
sx q[3];
rz(2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4158674) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(3.1126853) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(-0.96910563) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983738) q[0];
sx q[0];
rz(-1.507483) q[0];
sx q[0];
rz(-1.4988585) q[0];
rz(-pi) q[1];
rz(-14*pi/15) q[2];
sx q[2];
rz(-2.3839715) q[2];
sx q[2];
rz(2.9486738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9358878) q[1];
sx q[1];
rz(-1.74946) q[1];
sx q[1];
rz(-0.67358576) q[1];
rz(0.44090791) q[3];
sx q[3];
rz(-1.8280067) q[3];
sx q[3];
rz(0.5514901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4518296) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(-2.7815212) q[2];
rz(-1.5277956) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(-0.56267363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2071028) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(0.044152505) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(-2.2864441) q[2];
sx q[2];
rz(-2.6519041) q[2];
sx q[2];
rz(2.3408163) q[2];
rz(0.58892693) q[3];
sx q[3];
rz(-0.52994655) q[3];
sx q[3];
rz(0.51701057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
