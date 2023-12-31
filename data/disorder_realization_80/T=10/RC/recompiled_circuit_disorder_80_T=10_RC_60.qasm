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
rz(-0.84140468) q[1];
sx q[1];
rz(0.37766159) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084378622) q[0];
sx q[0];
rz(-0.14176653) q[0];
sx q[0];
rz(1.0300107) q[0];
rz(-1.471465) q[2];
sx q[2];
rz(-0.28495312) q[2];
sx q[2];
rz(-0.003665912) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.764896) q[1];
sx q[1];
rz(-0.88250676) q[1];
sx q[1];
rz(-0.98801686) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30794413) q[3];
sx q[3];
rz(-0.43711284) q[3];
sx q[3];
rz(1.2446158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1315786) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(-0.67260355) q[2];
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
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23068962) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(-0.22856523) q[0];
rz(0.16054343) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(0.28796089) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8156793) q[0];
sx q[0];
rz(-0.26198146) q[0];
sx q[0];
rz(2.322305) q[0];
x q[1];
rz(-2.1638853) q[2];
sx q[2];
rz(-1.8988673) q[2];
sx q[2];
rz(-1.6239945) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5353363) q[1];
sx q[1];
rz(-1.1621981) q[1];
sx q[1];
rz(0.6596843) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1045707) q[3];
sx q[3];
rz(-0.35211709) q[3];
sx q[3];
rz(-2.5196688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5471389) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(-1.6681558) q[2];
rz(0.93747059) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(0.37500769) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8163452) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(-2.8741799) q[0];
rz(1.7193517) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(0.95169383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0824453) q[0];
sx q[0];
rz(-1.7623616) q[0];
sx q[0];
rz(-0.092415718) q[0];
rz(-1.0551664) q[2];
sx q[2];
rz(-1.3736758) q[2];
sx q[2];
rz(-0.65161639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0905076) q[1];
sx q[1];
rz(-1.8790434) q[1];
sx q[1];
rz(0.80866637) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2610407) q[3];
sx q[3];
rz(-0.75062245) q[3];
sx q[3];
rz(0.39117884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(0.34417957) q[2];
rz(2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-2.5755431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0043871) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(-1.8970867) q[0];
rz(-3.0124774) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(0.37277645) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9739649) q[0];
sx q[0];
rz(-2.3229257) q[0];
sx q[0];
rz(-1.4139834) q[0];
rz(-pi) q[1];
rz(1.9984841) q[2];
sx q[2];
rz(-0.83979411) q[2];
sx q[2];
rz(1.0243624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.72153889) q[1];
sx q[1];
rz(-0.55100793) q[1];
sx q[1];
rz(2.0481471) q[1];
rz(1.0944486) q[3];
sx q[3];
rz(-1.8653449) q[3];
sx q[3];
rz(2.0457207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63561511) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(2.2367031) q[2];
rz(-2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(-2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6624517) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(1.8537846) q[0];
rz(1.6632535) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(2.6073661) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9984765) q[0];
sx q[0];
rz(-1.2623708) q[0];
sx q[0];
rz(-2.8240859) q[0];
rz(-pi) q[1];
x q[1];
rz(2.062837) q[2];
sx q[2];
rz(-1.8740219) q[2];
sx q[2];
rz(-0.505503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6215819) q[1];
sx q[1];
rz(-3.0428805) q[1];
sx q[1];
rz(-2.2732544) q[1];
rz(-1.4388496) q[3];
sx q[3];
rz(-2.2852516) q[3];
sx q[3];
rz(-2.7150142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(-0.14349288) q[2];
rz(1.3714553) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(0.21970704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7025529) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(0.23705661) q[0];
rz(1.2409695) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(1.7664849) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8538044) q[0];
sx q[0];
rz(-1.8332991) q[0];
sx q[0];
rz(-0.35065035) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4985956) q[2];
sx q[2];
rz(-2.318396) q[2];
sx q[2];
rz(-2.2526134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.68723893) q[1];
sx q[1];
rz(-0.41077405) q[1];
sx q[1];
rz(0.06128581) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9938019) q[3];
sx q[3];
rz(-1.1173964) q[3];
sx q[3];
rz(-0.51844937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6087626) q[2];
sx q[2];
rz(-0.71981788) q[2];
sx q[2];
rz(-0.14870816) q[2];
rz(-0.016629774) q[3];
sx q[3];
rz(-2.7745268) q[3];
sx q[3];
rz(-2.9490411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6454813) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(-2.2684229) q[0];
rz(2.219615) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(3.057664) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31939313) q[0];
sx q[0];
rz(-2.5376352) q[0];
sx q[0];
rz(1.1711867) q[0];
rz(-pi) q[1];
rz(-0.35328816) q[2];
sx q[2];
rz(-1.366426) q[2];
sx q[2];
rz(-0.7962966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2578656) q[1];
sx q[1];
rz(-2.1818433) q[1];
sx q[1];
rz(-0.74531598) q[1];
rz(-pi) q[2];
rz(-3.0725669) q[3];
sx q[3];
rz(-0.8419753) q[3];
sx q[3];
rz(-0.1544827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38368791) q[2];
sx q[2];
rz(-0.44054511) q[2];
sx q[2];
rz(-2.596358) q[2];
rz(0.43045726) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(-2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(1.9301201) q[0];
rz(-0.97310549) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(2.1957695) q[1];
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
rz(-0.68265712) q[2];
sx q[2];
rz(-2.1652512) q[2];
sx q[2];
rz(2.1669831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.91499) q[1];
sx q[1];
rz(-1.9493002) q[1];
sx q[1];
rz(-1.6016017) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4754727) q[3];
sx q[3];
rz(-1.5770116) q[3];
sx q[3];
rz(-2.4339649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(2.8611709) q[2];
rz(-0.60493207) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(-2.159507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9873001) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(-0.61532414) q[0];
rz(0.92957169) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(0.5756793) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3424123) q[0];
sx q[0];
rz(-0.77501446) q[0];
sx q[0];
rz(2.4393625) q[0];
rz(-pi) q[1];
rz(2.0909537) q[2];
sx q[2];
rz(-2.9165604) q[2];
sx q[2];
rz(-2.2573543) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5513735) q[1];
sx q[1];
rz(-0.11034036) q[1];
sx q[1];
rz(-1.3361841) q[1];
rz(-pi) q[2];
rz(1.8768164) q[3];
sx q[3];
rz(-2.5796606) q[3];
sx q[3];
rz(-2.7276873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3118887) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(-0.021961948) q[2];
rz(0.17045505) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(-2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72572529) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(2.5073994) q[0];
rz(3.1126853) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(2.172487) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0646107) q[0];
sx q[0];
rz(-1.4990028) q[0];
sx q[0];
rz(3.0781156) q[0];
rz(-2.3949941) q[2];
sx q[2];
rz(-1.4274297) q[2];
sx q[2];
rz(1.6105086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9358878) q[1];
sx q[1];
rz(-1.74946) q[1];
sx q[1];
rz(0.67358576) q[1];
rz(-1.2877527) q[3];
sx q[3];
rz(-1.9962365) q[3];
sx q[3];
rz(0.89983672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4518296) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(-0.36007145) q[2];
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
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2071028) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(3.0974401) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(-1.953223) q[2];
sx q[2];
rz(-1.8845176) q[2];
sx q[2];
rz(0.11558576) q[2];
rz(1.885407) q[3];
sx q[3];
rz(-2.0046069) q[3];
sx q[3];
rz(-0.14179695) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
