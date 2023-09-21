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
rz(1.1893907) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(2.300188) q[1];
sx q[1];
rz(9.0471164) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46086754) q[0];
sx q[0];
rz(-1.4493677) q[0];
sx q[0];
rz(-0.073343883) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8544191) q[2];
sx q[2];
rz(-1.5986773) q[2];
sx q[2];
rz(1.6698128) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1805815) q[1];
sx q[1];
rz(-0.86997021) q[1];
sx q[1];
rz(-2.5518774) q[1];
x q[2];
rz(-0.30794413) q[3];
sx q[3];
rz(-2.7044798) q[3];
sx q[3];
rz(1.8969769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1315786) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(-0.67260355) q[2];
rz(-0.16942313) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(-1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.910903) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(2.9130274) q[0];
rz(0.16054343) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(0.28796089) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8156793) q[0];
sx q[0];
rz(-2.8796112) q[0];
sx q[0];
rz(0.81928763) q[0];
rz(-0.9777074) q[2];
sx q[2];
rz(-1.8988673) q[2];
sx q[2];
rz(-1.5175982) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.60625633) q[1];
sx q[1];
rz(-1.1621981) q[1];
sx q[1];
rz(0.6596843) q[1];
x q[2];
rz(-1.6091299) q[3];
sx q[3];
rz(-1.2206856) q[3];
sx q[3];
rz(2.4083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5471389) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(-1.6681558) q[2];
rz(2.2041221) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8163452) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(0.26741272) q[0];
rz(-1.422241) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(-0.95169383) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6355977) q[0];
sx q[0];
rz(-1.6615168) q[0];
sx q[0];
rz(1.3784301) q[0];
rz(-pi) q[1];
rz(2.9159413) q[2];
sx q[2];
rz(-1.0661085) q[2];
sx q[2];
rz(0.80863189) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.051085) q[1];
sx q[1];
rz(-1.2625492) q[1];
sx q[1];
rz(-2.3329263) q[1];
rz(2.2610407) q[3];
sx q[3];
rz(-2.3909702) q[3];
sx q[3];
rz(-2.7504138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(2.7974131) q[2];
rz(2.2551645) q[3];
sx q[3];
rz(-0.27799806) q[3];
sx q[3];
rz(-0.56604958) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043871) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(1.8970867) q[0];
rz(-0.1291153) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(2.7688162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059797212) q[0];
sx q[0];
rz(-0.76515388) q[0];
sx q[0];
rz(0.16539903) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9984841) q[2];
sx q[2];
rz(-2.3017985) q[2];
sx q[2];
rz(-1.0243624) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4200538) q[1];
sx q[1];
rz(-2.5905847) q[1];
sx q[1];
rz(2.0481471) q[1];
rz(-pi) q[2];
rz(-2.0471441) q[3];
sx q[3];
rz(-1.2762478) q[3];
sx q[3];
rz(1.095872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63561511) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(-0.90488952) q[2];
rz(-0.44057009) q[3];
sx q[3];
rz(-0.4959271) q[3];
sx q[3];
rz(0.99159616) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47914094) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(-1.8537846) q[0];
rz(-1.4783391) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(-2.6073661) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5271082) q[0];
sx q[0];
rz(-1.2687578) q[0];
sx q[0];
rz(-1.8943647) q[0];
x q[1];
rz(-1.0787557) q[2];
sx q[2];
rz(-1.8740219) q[2];
sx q[2];
rz(-0.505503) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2248762) q[1];
sx q[1];
rz(-1.4955048) q[1];
sx q[1];
rz(0.063898357) q[1];
rz(-pi) q[2];
rz(1.7027431) q[3];
sx q[3];
rz(-2.2852516) q[3];
sx q[3];
rz(0.42657846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6405032) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(-0.14349288) q[2];
rz(1.3714553) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(-2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7025529) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-2.904536) q[0];
rz(-1.2409695) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(-1.3751078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2415337) q[0];
sx q[0];
rz(-0.43474475) q[0];
sx q[0];
rz(-2.4777806) q[0];
rz(-pi) q[1];
rz(3.0639406) q[2];
sx q[2];
rz(-2.3911871) q[2];
sx q[2];
rz(0.99496182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.68723893) q[1];
sx q[1];
rz(-0.41077405) q[1];
sx q[1];
rz(-0.06128581) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49075134) q[3];
sx q[3];
rz(-1.1928344) q[3];
sx q[3];
rz(2.2839387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6087626) q[2];
sx q[2];
rz(-0.71981788) q[2];
sx q[2];
rz(0.14870816) q[2];
rz(0.016629774) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49611133) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(2.2684229) q[0];
rz(0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(0.08392863) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31939313) q[0];
sx q[0];
rz(-0.60395741) q[0];
sx q[0];
rz(1.970406) q[0];
rz(1.7882118) q[2];
sx q[2];
rz(-1.916421) q[2];
sx q[2];
rz(2.2923922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9675688) q[1];
sx q[1];
rz(-2.1597383) q[1];
sx q[1];
rz(2.3322361) q[1];
rz(-2.300802) q[3];
sx q[3];
rz(-1.6222686) q[3];
sx q[3];
rz(1.7712902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7579047) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(-0.54523462) q[2];
rz(-0.43045726) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(0.30495131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-3.1261303) q[0];
sx q[0];
rz(-1.9301201) q[0];
rz(0.97310549) q[1];
sx q[1];
rz(-0.5935697) q[1];
sx q[1];
rz(2.1957695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57293939) q[0];
sx q[0];
rz(-1.4927673) q[0];
sx q[0];
rz(2.9201939) q[0];
x q[1];
rz(2.321645) q[2];
sx q[2];
rz(-0.87265271) q[2];
sx q[2];
rz(-1.1993711) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.91499) q[1];
sx q[1];
rz(-1.9493002) q[1];
sx q[1];
rz(1.539991) q[1];
rz(-pi) q[2];
rz(-0.010057851) q[3];
sx q[3];
rz(-0.66614449) q[3];
sx q[3];
rz(-0.87107623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(-2.8611709) q[2];
rz(2.5366606) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9873001) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(-2.5262685) q[0];
rz(0.92957169) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(-2.5659134) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070411365) q[0];
sx q[0];
rz(-1.0072664) q[0];
sx q[0];
rz(1.0067183) q[0];
rz(-1.766878) q[2];
sx q[2];
rz(-1.6819281) q[2];
sx q[2];
rz(-2.9642504) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3153783) q[1];
sx q[1];
rz(-1.678102) q[1];
sx q[1];
rz(3.1158434) q[1];
rz(-pi) q[2];
rz(2.9541295) q[3];
sx q[3];
rz(-1.0378569) q[3];
sx q[3];
rz(-0.77123469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82970396) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(3.1196307) q[2];
rz(2.9711376) q[3];
sx q[3];
rz(-1.0237834) q[3];
sx q[3];
rz(0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4158674) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(2.5073994) q[0];
rz(3.1126853) q[1];
sx q[1];
rz(-2.3529265) q[1];
sx q[1];
rz(0.96910563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6432188) q[0];
sx q[0];
rz(-1.6341097) q[0];
sx q[0];
rz(1.6427342) q[0];
x q[1];
rz(2.3949941) q[2];
sx q[2];
rz(-1.4274297) q[2];
sx q[2];
rz(-1.6105086) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9358878) q[1];
sx q[1];
rz(-1.74946) q[1];
sx q[1];
rz(-0.67358576) q[1];
rz(-pi) q[2];
x q[2];
rz(1.85384) q[3];
sx q[3];
rz(-1.9962365) q[3];
sx q[3];
rz(-2.2417559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6897631) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(2.7815212) q[2];
rz(-1.6137971) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(2.578919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.9344899) q[0];
sx q[0];
rz(-1.5710545) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(-3.0974401) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(-2.2864441) q[2];
sx q[2];
rz(-2.6519041) q[2];
sx q[2];
rz(2.3408163) q[2];
rz(-2.5526657) q[3];
sx q[3];
rz(-0.52994655) q[3];
sx q[3];
rz(0.51701057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
