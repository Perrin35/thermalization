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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6807251) q[0];
sx q[0];
rz(-1.6922249) q[0];
sx q[0];
rz(-3.0682488) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2871735) q[2];
sx q[2];
rz(-1.5986773) q[2];
sx q[2];
rz(1.6698128) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.764896) q[1];
sx q[1];
rz(-0.88250676) q[1];
sx q[1];
rz(2.1535758) q[1];
rz(2.8336485) q[3];
sx q[3];
rz(-0.43711284) q[3];
sx q[3];
rz(-1.8969769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.010014023) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(2.4689891) q[2];
rz(-2.9721695) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(-1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.910903) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(-2.9130274) q[0];
rz(-0.16054343) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(-0.28796089) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1624958) q[0];
sx q[0];
rz(-1.7485577) q[0];
sx q[0];
rz(-1.7642679) q[0];
rz(2.1177883) q[2];
sx q[2];
rz(-2.4734801) q[2];
sx q[2];
rz(2.748865) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6504753) q[1];
sx q[1];
rz(-0.75956356) q[1];
sx q[1];
rz(-0.61504765) q[1];
x q[2];
rz(-1.6091299) q[3];
sx q[3];
rz(-1.2206856) q[3];
sx q[3];
rz(-0.73327524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5471389) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(1.6681558) q[2];
rz(2.2041221) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(-0.37500769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8163452) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(0.26741272) q[0];
rz(-1.422241) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(-0.95169383) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0824453) q[0];
sx q[0];
rz(-1.7623616) q[0];
sx q[0];
rz(0.092415718) q[0];
x q[1];
rz(-2.9159413) q[2];
sx q[2];
rz(-1.0661085) q[2];
sx q[2];
rz(-0.80863189) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23783319) q[1];
sx q[1];
rz(-2.2888498) q[1];
sx q[1];
rz(-0.41463931) q[1];
x q[2];
rz(-2.6056616) q[3];
sx q[3];
rz(-2.1246353) q[3];
sx q[3];
rz(-1.2371847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13016985) q[2];
sx q[2];
rz(-1.6458076) q[2];
sx q[2];
rz(2.7974131) q[2];
rz(-2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043871) q[0];
sx q[0];
rz(-1.4973649) q[0];
sx q[0];
rz(1.244506) q[0];
rz(0.1291153) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(-0.37277645) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0817954) q[0];
sx q[0];
rz(-2.3764388) q[0];
sx q[0];
rz(0.16539903) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7083621) q[2];
sx q[2];
rz(-0.82651143) q[2];
sx q[2];
rz(-2.7162958) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8744295) q[1];
sx q[1];
rz(-1.0871372) q[1];
sx q[1];
rz(-0.27515414) q[1];
rz(-0.9863015) q[3];
sx q[3];
rz(-2.5875475) q[3];
sx q[3];
rz(-0.98741764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5059775) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(0.90488952) q[2];
rz(-0.44057009) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6624517) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.2878081) q[0];
rz(-1.6632535) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(-0.53422654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5271082) q[0];
sx q[0];
rz(-1.8728349) q[0];
sx q[0];
rz(1.8943647) q[0];
rz(-pi) q[1];
x q[1];
rz(2.062837) q[2];
sx q[2];
rz(-1.2675708) q[2];
sx q[2];
rz(-2.6360896) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.91671645) q[1];
sx q[1];
rz(-1.6460878) q[1];
sx q[1];
rz(3.0776943) q[1];
x q[2];
rz(-1.7027431) q[3];
sx q[3];
rz(-0.85634106) q[3];
sx q[3];
rz(-2.7150142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(-2.9980998) q[2];
rz(1.3714553) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(0.21970704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7025529) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-0.23705661) q[0];
rz(1.2409695) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(-1.7664849) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9532167) q[0];
sx q[0];
rz(-1.9089451) q[0];
sx q[0];
rz(1.8494649) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74890045) q[2];
sx q[2];
rz(-1.5178711) q[2];
sx q[2];
rz(2.5089094) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82735862) q[1];
sx q[1];
rz(-1.595256) q[1];
sx q[1];
rz(0.41008653) q[1];
rz(-pi) q[2];
rz(-2.4414805) q[3];
sx q[3];
rz(-2.5317319) q[3];
sx q[3];
rz(1.8240579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6087626) q[2];
sx q[2];
rz(-0.71981788) q[2];
sx q[2];
rz(-2.9928845) q[2];
rz(3.1249629) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(-0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-0.9219777) q[1];
sx q[1];
rz(-1.1154122) q[1];
sx q[1];
rz(-3.057664) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5859563) q[0];
sx q[0];
rz(-1.3480098) q[0];
sx q[0];
rz(-1.0046093) q[0];
rz(1.3533808) q[2];
sx q[2];
rz(-1.916421) q[2];
sx q[2];
rz(0.84920041) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9675688) q[1];
sx q[1];
rz(-0.98185437) q[1];
sx q[1];
rz(-0.80935652) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6478959) q[3];
sx q[3];
rz(-0.73148433) q[3];
sx q[3];
rz(2.8836718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38368791) q[2];
sx q[2];
rz(-0.44054511) q[2];
sx q[2];
rz(-0.54523462) q[2];
rz(-2.7111354) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(0.30495131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.6298744) q[0];
sx q[0];
rz(-3.1261303) q[0];
sx q[0];
rz(1.9301201) q[0];
rz(0.97310549) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(-2.1957695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0154008) q[0];
sx q[0];
rz(-1.7915103) q[0];
sx q[0];
rz(1.4908233) q[0];
rz(0.68265712) q[2];
sx q[2];
rz(-0.97634146) q[2];
sx q[2];
rz(2.1669831) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1434053) q[1];
sx q[1];
rz(-2.7618976) q[1];
sx q[1];
rz(3.0642964) q[1];
rz(3.1315348) q[3];
sx q[3];
rz(-2.4754482) q[3];
sx q[3];
rz(-2.2705164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(0.28042173) q[2];
rz(0.60493207) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9873001) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(-0.61532414) q[0];
rz(0.92957169) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(2.5659134) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8262779) q[0];
sx q[0];
rz(-1.1018254) q[0];
sx q[0];
rz(-2.4995575) q[0];
rz(-3.028308) q[2];
sx q[2];
rz(-1.3759398) q[2];
sx q[2];
rz(1.4154797) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74734028) q[1];
sx q[1];
rz(-1.5451952) q[1];
sx q[1];
rz(1.4634553) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2647763) q[3];
sx q[3];
rz(-2.5796606) q[3];
sx q[3];
rz(0.41390536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82970396) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(-3.1196307) q[2];
rz(-2.9711376) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(0.26486614) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72572529) q[0];
sx q[0];
rz(-0.9265582) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(-3.1126853) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(0.96910563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64810753) q[0];
sx q[0];
rz(-3.0457975) q[0];
sx q[0];
rz(2.2936054) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74659851) q[2];
sx q[2];
rz(-1.714163) q[2];
sx q[2];
rz(1.5310841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2057048) q[1];
sx q[1];
rz(-1.3921326) q[1];
sx q[1];
rz(0.67358576) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7006847) q[3];
sx q[3];
rz(-1.8280067) q[3];
sx q[3];
rz(2.5901026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4518296) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(-2.7815212) q[2];
rz(-1.6137971) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(-2.578919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9344899) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(-0.044152505) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(-1.953223) q[2];
sx q[2];
rz(-1.8845176) q[2];
sx q[2];
rz(0.11558576) q[2];
rz(0.45331656) q[3];
sx q[3];
rz(-1.286187) q[3];
sx q[3];
rz(-1.5766531) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
