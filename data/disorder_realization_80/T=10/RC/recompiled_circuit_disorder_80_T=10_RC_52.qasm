OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7497082) q[0];
sx q[0];
rz(-2.9449129) q[0];
sx q[0];
rz(-1.1893907) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(-2.7639311) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1188287) q[0];
sx q[0];
rz(-1.4979935) q[0];
sx q[0];
rz(1.692549) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.112552) q[2];
sx q[2];
rz(-1.8543058) q[2];
sx q[2];
rz(3.0344506) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9610112) q[1];
sx q[1];
rz(-0.86997021) q[1];
sx q[1];
rz(0.58971528) q[1];
rz(-pi) q[2];
rz(1.430106) q[3];
sx q[3];
rz(-1.9860387) q[3];
sx q[3];
rz(1.5822441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1315786) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(-2.4689891) q[2];
rz(-2.9721695) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.910903) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(0.22856523) q[0];
rz(0.16054343) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(-0.28796089) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1624958) q[0];
sx q[0];
rz(-1.7485577) q[0];
sx q[0];
rz(1.3773247) q[0];
rz(0.38950133) q[2];
sx q[2];
rz(-1.0132388) q[2];
sx q[2];
rz(-2.8745289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60625633) q[1];
sx q[1];
rz(-1.9793946) q[1];
sx q[1];
rz(0.6596843) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.037022) q[3];
sx q[3];
rz(-2.7894756) q[3];
sx q[3];
rz(0.62192384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59445375) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(-1.6681558) q[2];
rz(0.93747059) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8163452) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(1.7193517) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(-2.1898988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62943447) q[0];
sx q[0];
rz(-2.929147) q[0];
sx q[0];
rz(-2.0149219) q[0];
rz(-1.9556324) q[2];
sx q[2];
rz(-0.54883146) q[2];
sx q[2];
rz(1.8897111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9037595) q[1];
sx q[1];
rz(-0.85274285) q[1];
sx q[1];
rz(0.41463931) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94727256) q[3];
sx q[3];
rz(-2.0200649) q[3];
sx q[3];
rz(0.63637966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13016985) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(0.34417957) q[2];
rz(-2.2551645) q[3];
sx q[3];
rz(-0.27799806) q[3];
sx q[3];
rz(-2.5755431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0043871) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(1.244506) q[0];
rz(-3.0124774) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(-2.7688162) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5107721) q[0];
sx q[0];
rz(-1.6850867) q[0];
sx q[0];
rz(-0.75829102) q[0];
rz(1.9984841) q[2];
sx q[2];
rz(-0.83979411) q[2];
sx q[2];
rz(1.0243624) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2671632) q[1];
sx q[1];
rz(-1.0871372) q[1];
sx q[1];
rz(-0.27515414) q[1];
rz(-pi) q[2];
rz(2.0471441) q[3];
sx q[3];
rz(-1.2762478) q[3];
sx q[3];
rz(-1.095872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5059775) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(0.90488952) q[2];
rz(-2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(0.99159616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6624517) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(1.8537846) q[0];
rz(-1.4783391) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(-0.53422654) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5271082) q[0];
sx q[0];
rz(-1.8728349) q[0];
sx q[0];
rz(-1.2472279) q[0];
rz(-0.9858426) q[2];
sx q[2];
rz(-0.57136977) q[2];
sx q[2];
rz(-0.55703288) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91671645) q[1];
sx q[1];
rz(-1.4955048) q[1];
sx q[1];
rz(-0.063898357) q[1];
rz(2.4228135) q[3];
sx q[3];
rz(-1.4712417) q[3];
sx q[3];
rz(1.0574785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6405032) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(2.9980998) q[2];
rz(1.7701373) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(0.21970704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.4390398) q[0];
sx q[0];
rz(-2.2655903) q[0];
sx q[0];
rz(2.904536) q[0];
rz(1.9006231) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(1.7664849) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2415337) q[0];
sx q[0];
rz(-2.7068479) q[0];
sx q[0];
rz(0.66381201) q[0];
rz(-0.077652046) q[2];
sx q[2];
rz(-2.3911871) q[2];
sx q[2];
rz(-2.1466308) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75406972) q[1];
sx q[1];
rz(-1.9807528) q[1];
sx q[1];
rz(-1.5974664) q[1];
rz(-pi) q[2];
rz(0.49075134) q[3];
sx q[3];
rz(-1.1928344) q[3];
sx q[3];
rz(-0.85765391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.53283006) q[2];
sx q[2];
rz(-2.4217748) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49611133) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(-0.87316978) q[0];
rz(0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(-3.057664) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.98691) q[0];
sx q[0];
rz(-1.0202408) q[0];
sx q[0];
rz(-0.26225342) q[0];
rz(-0.53972466) q[2];
sx q[2];
rz(-2.735609) q[2];
sx q[2];
rz(2.8702131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1740239) q[1];
sx q[1];
rz(-2.1597383) q[1];
sx q[1];
rz(2.3322361) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6478959) q[3];
sx q[3];
rz(-0.73148433) q[3];
sx q[3];
rz(-0.25792083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7579047) q[2];
sx q[2];
rz(-0.44054511) q[2];
sx q[2];
rz(0.54523462) q[2];
rz(-2.7111354) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(-0.30495131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6298744) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(1.9301201) q[0];
rz(2.1684872) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(2.1957695) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57293939) q[0];
sx q[0];
rz(-1.6488254) q[0];
sx q[0];
rz(2.9201939) q[0];
rz(-pi) q[1];
rz(0.81994762) q[2];
sx q[2];
rz(-0.87265271) q[2];
sx q[2];
rz(-1.9422216) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7860124) q[1];
sx q[1];
rz(-1.5994206) q[1];
sx q[1];
rz(-2.7629258) q[1];
rz(-pi) q[2];
rz(0.010057851) q[3];
sx q[3];
rz(-0.66614449) q[3];
sx q[3];
rz(-2.2705164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(0.28042173) q[2];
rz(-2.5366606) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9873001) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(-2.5262685) q[0];
rz(-0.92957169) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(2.5659134) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3153148) q[0];
sx q[0];
rz(-1.1018254) q[0];
sx q[0];
rz(2.4995575) q[0];
x q[1];
rz(-3.028308) q[2];
sx q[2];
rz(-1.7656529) q[2];
sx q[2];
rz(-1.4154797) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74734028) q[1];
sx q[1];
rz(-1.5963975) q[1];
sx q[1];
rz(-1.4634553) q[1];
x q[2];
rz(-1.8768164) q[3];
sx q[3];
rz(-0.56193202) q[3];
sx q[3];
rz(-2.7276873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3118887) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(-0.021961948) q[2];
rz(-2.9711376) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(-2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4158674) q[0];
sx q[0];
rz(-0.9265582) q[0];
sx q[0];
rz(-2.5073994) q[0];
rz(0.028907396) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(-2.172487) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6432188) q[0];
sx q[0];
rz(-1.507483) q[0];
sx q[0];
rz(1.6427342) q[0];
rz(0.74659851) q[2];
sx q[2];
rz(-1.4274297) q[2];
sx q[2];
rz(-1.5310841) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22419588) q[1];
sx q[1];
rz(-2.2317413) q[1];
sx q[1];
rz(1.7978653) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44090791) q[3];
sx q[3];
rz(-1.3135859) q[3];
sx q[3];
rz(-0.5514901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4518296) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(0.36007145) q[2];
rz(-1.6137971) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(-2.578919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
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
rz(2.2864441) q[2];
sx q[2];
rz(-0.48968857) q[2];
sx q[2];
rz(-0.80077632) q[2];
rz(-0.45331656) q[3];
sx q[3];
rz(-1.8554056) q[3];
sx q[3];
rz(1.5649395) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
