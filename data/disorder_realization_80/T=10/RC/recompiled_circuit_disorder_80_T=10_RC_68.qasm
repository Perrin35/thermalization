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
rz(0.37766159) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6807251) q[0];
sx q[0];
rz(-1.6922249) q[0];
sx q[0];
rz(-0.073343883) q[0];
x q[1];
rz(1.8544191) q[2];
sx q[2];
rz(-1.5986773) q[2];
sx q[2];
rz(1.6698128) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37669668) q[1];
sx q[1];
rz(-0.88250676) q[1];
sx q[1];
rz(-0.98801686) q[1];
x q[2];
rz(-1.7114867) q[3];
sx q[3];
rz(-1.9860387) q[3];
sx q[3];
rz(1.5822441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.010014023) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(-2.4689891) q[2];
rz(-2.9721695) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23068962) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(-2.9130274) q[0];
rz(-2.9810492) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(-2.8536318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44293091) q[0];
sx q[0];
rz(-1.7611815) q[0];
sx q[0];
rz(2.9605244) q[0];
x q[1];
rz(2.1638853) q[2];
sx q[2];
rz(-1.8988673) q[2];
sx q[2];
rz(1.6239945) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5353363) q[1];
sx q[1];
rz(-1.1621981) q[1];
sx q[1];
rz(2.4819083) q[1];
rz(-pi) q[2];
rz(0.1045707) q[3];
sx q[3];
rz(-2.7894756) q[3];
sx q[3];
rz(0.62192384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59445375) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(-1.4734369) q[2];
rz(0.93747059) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(-2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32524747) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(-1.7193517) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(0.95169383) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0824453) q[0];
sx q[0];
rz(-1.379231) q[0];
sx q[0];
rz(3.0491769) q[0];
rz(1.1859602) q[2];
sx q[2];
rz(-0.54883146) q[2];
sx q[2];
rz(1.8897111) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.051085) q[1];
sx q[1];
rz(-1.2625492) q[1];
sx q[1];
rz(0.80866637) q[1];
x q[2];
rz(2.6056616) q[3];
sx q[3];
rz(-1.0169573) q[3];
sx q[3];
rz(-1.2371847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13720559) q[0];
sx q[0];
rz(-1.4973649) q[0];
sx q[0];
rz(-1.244506) q[0];
rz(3.0124774) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(-2.7688162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059797212) q[0];
sx q[0];
rz(-2.3764388) q[0];
sx q[0];
rz(-2.9761936) q[0];
rz(-pi) q[1];
rz(-0.77809019) q[2];
sx q[2];
rz(-1.8847244) q[2];
sx q[2];
rz(-2.2997466) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8744295) q[1];
sx q[1];
rz(-1.0871372) q[1];
sx q[1];
rz(-2.8664385) q[1];
x q[2];
rz(2.8126206) q[3];
sx q[3];
rz(-2.0250468) q[3];
sx q[3];
rz(2.8153552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63561511) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(-0.90488952) q[2];
rz(-0.44057009) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(-0.99159616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47914094) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(-1.2878081) q[0];
rz(-1.6632535) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(0.53422654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68191093) q[0];
sx q[0];
rz(-0.43897438) q[0];
sx q[0];
rz(-0.79553332) q[0];
rz(2.1557501) q[2];
sx q[2];
rz(-0.57136977) q[2];
sx q[2];
rz(-0.55703288) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4923258) q[1];
sx q[1];
rz(-1.6345134) q[1];
sx q[1];
rz(-1.4953514) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7027431) q[3];
sx q[3];
rz(-2.2852516) q[3];
sx q[3];
rz(0.42657846) q[3];
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
rz(-1.7701373) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(-2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4390398) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-2.904536) q[0];
rz(1.2409695) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(1.3751078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8538044) q[0];
sx q[0];
rz(-1.8332991) q[0];
sx q[0];
rz(-2.7909423) q[0];
rz(-pi) q[1];
rz(-1.4985956) q[2];
sx q[2];
rz(-0.82319665) q[2];
sx q[2];
rz(2.2526134) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.314234) q[1];
sx q[1];
rz(-1.595256) q[1];
sx q[1];
rz(-2.7315061) q[1];
rz(-0.49075134) q[3];
sx q[3];
rz(-1.1928344) q[3];
sx q[3];
rz(0.85765391) q[3];
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
rz(-0.016629774) q[3];
sx q[3];
rz(-2.7745268) q[3];
sx q[3];
rz(0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6454813) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(2.2684229) q[0];
rz(-0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(3.057664) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8221995) q[0];
sx q[0];
rz(-0.60395741) q[0];
sx q[0];
rz(-1.970406) q[0];
rz(-1.3533808) q[2];
sx q[2];
rz(-1.916421) q[2];
sx q[2];
rz(2.2923922) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9675688) q[1];
sx q[1];
rz(-0.98185437) q[1];
sx q[1];
rz(-0.80935652) q[1];
rz(-pi) q[2];
x q[2];
rz(0.069025741) q[3];
sx q[3];
rz(-0.8419753) q[3];
sx q[3];
rz(2.98711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7579047) q[2];
sx q[2];
rz(-0.44054511) q[2];
sx q[2];
rz(0.54523462) q[2];
rz(0.43045726) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(-2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6298744) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(1.2114725) q[0];
rz(0.97310549) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(-2.1957695) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66447542) q[0];
sx q[0];
rz(-0.23453377) q[0];
sx q[0];
rz(0.3420591) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4589355) q[2];
sx q[2];
rz(-2.1652512) q[2];
sx q[2];
rz(-0.97460954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.91499) q[1];
sx q[1];
rz(-1.1922925) q[1];
sx q[1];
rz(-1.539991) q[1];
x q[2];
rz(0.66611992) q[3];
sx q[3];
rz(-1.5770116) q[3];
sx q[3];
rz(0.70762779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(-2.8611709) q[2];
rz(-0.60493207) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1542926) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(-0.61532414) q[0];
rz(2.212021) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(-2.5659134) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0711813) q[0];
sx q[0];
rz(-1.0072664) q[0];
sx q[0];
rz(2.1348743) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11328463) q[2];
sx q[2];
rz(-1.7656529) q[2];
sx q[2];
rz(1.726113) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8262144) q[1];
sx q[1];
rz(-1.678102) q[1];
sx q[1];
rz(-3.1158434) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2647763) q[3];
sx q[3];
rz(-0.56193202) q[3];
sx q[3];
rz(-0.41390536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3118887) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(-3.1196307) q[2];
rz(2.9711376) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72572529) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(-2.5073994) q[0];
rz(-3.1126853) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(-2.172487) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64810753) q[0];
sx q[0];
rz(-3.0457975) q[0];
sx q[0];
rz(0.84798725) q[0];
rz(-pi) q[1];
rz(-2.3949941) q[2];
sx q[2];
rz(-1.714163) q[2];
sx q[2];
rz(1.5310841) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2057048) q[1];
sx q[1];
rz(-1.74946) q[1];
sx q[1];
rz(-2.4680069) q[1];
rz(-1.85384) q[3];
sx q[3];
rz(-1.1453562) q[3];
sx q[3];
rz(-2.2417559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6897631) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(-2.7815212) q[2];
rz(1.5277956) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(-2.578919) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
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
rz(-0.85514851) q[2];
sx q[2];
rz(-0.48968857) q[2];
sx q[2];
rz(-0.80077632) q[2];
rz(-1.885407) q[3];
sx q[3];
rz(-1.1369858) q[3];
sx q[3];
rz(2.9997957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];