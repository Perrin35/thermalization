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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6807251) q[0];
sx q[0];
rz(-1.6922249) q[0];
sx q[0];
rz(3.0682488) q[0];
rz(-pi) q[1];
rz(1.8544191) q[2];
sx q[2];
rz(-1.5429153) q[2];
sx q[2];
rz(1.4717799) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.764896) q[1];
sx q[1];
rz(-2.2590859) q[1];
sx q[1];
rz(-0.98801686) q[1];
x q[2];
rz(1.430106) q[3];
sx q[3];
rz(-1.9860387) q[3];
sx q[3];
rz(1.5822441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.010014023) q[2];
sx q[2];
rz(-2.6476314) q[2];
sx q[2];
rz(-2.4689891) q[2];
rz(0.16942313) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23068962) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(-0.22856523) q[0];
rz(-0.16054343) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(-2.8536318) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32591336) q[0];
sx q[0];
rz(-0.26198146) q[0];
sx q[0];
rz(2.322305) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9777074) q[2];
sx q[2];
rz(-1.8988673) q[2];
sx q[2];
rz(1.5175982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2634695) q[1];
sx q[1];
rz(-2.1681004) q[1];
sx q[1];
rz(-1.0695446) q[1];
x q[2];
rz(-1.6091299) q[3];
sx q[3];
rz(-1.2206856) q[3];
sx q[3];
rz(2.4083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(1.4734369) q[2];
rz(-0.93747059) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(-2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8163452) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(0.26741272) q[0];
rz(1.422241) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(-2.1898988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6355977) q[0];
sx q[0];
rz(-1.6615168) q[0];
sx q[0];
rz(1.7631625) q[0];
rz(-pi) q[1];
rz(-1.0551664) q[2];
sx q[2];
rz(-1.7679169) q[2];
sx q[2];
rz(-2.4899763) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9037595) q[1];
sx q[1];
rz(-2.2888498) q[1];
sx q[1];
rz(-0.41463931) q[1];
rz(-pi) q[2];
rz(0.94727256) q[3];
sx q[3];
rz(-2.0200649) q[3];
sx q[3];
rz(-2.505213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
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
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
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
rz(-1.244506) q[0];
rz(3.0124774) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(-2.7688162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9739649) q[0];
sx q[0];
rz(-0.81866696) q[0];
sx q[0];
rz(1.7276093) q[0];
rz(-pi) q[1];
rz(2.3635025) q[2];
sx q[2];
rz(-1.2568682) q[2];
sx q[2];
rz(-0.84184605) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43416926) q[1];
sx q[1];
rz(-1.8137099) q[1];
sx q[1];
rz(1.0711819) q[1];
rz(-2.0471441) q[3];
sx q[3];
rz(-1.2762478) q[3];
sx q[3];
rz(-2.0457207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5059775) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(0.90488952) q[2];
rz(2.7010226) q[3];
sx q[3];
rz(-0.4959271) q[3];
sx q[3];
rz(0.99159616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6624517) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(1.8537846) q[0];
rz(1.4783391) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(-0.53422654) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4596817) q[0];
sx q[0];
rz(-2.7026183) q[0];
sx q[0];
rz(2.3460593) q[0];
rz(-pi) q[1];
rz(-2.8004831) q[2];
sx q[2];
rz(-1.1030536) q[2];
sx q[2];
rz(1.2240006) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6215819) q[1];
sx q[1];
rz(-0.098712155) q[1];
sx q[1];
rz(-2.2732544) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4388496) q[3];
sx q[3];
rz(-2.2852516) q[3];
sx q[3];
rz(0.42657846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6405032) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(0.14349288) q[2];
rz(-1.7701373) q[3];
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
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4390398) q[0];
sx q[0];
rz(-2.2655903) q[0];
sx q[0];
rz(0.23705661) q[0];
rz(1.2409695) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(1.7664849) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8538044) q[0];
sx q[0];
rz(-1.3082936) q[0];
sx q[0];
rz(0.35065035) q[0];
rz(-pi) q[1];
rz(-1.6429971) q[2];
sx q[2];
rz(-2.318396) q[2];
sx q[2];
rz(-0.88897926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4543537) q[1];
sx q[1];
rz(-0.41077405) q[1];
sx q[1];
rz(-3.0803068) q[1];
rz(2.6508413) q[3];
sx q[3];
rz(-1.9487582) q[3];
sx q[3];
rz(2.2839387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53283006) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6454813) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(-0.87316978) q[0];
rz(2.219615) q[1];
sx q[1];
rz(-1.1154122) q[1];
sx q[1];
rz(-3.057664) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15468266) q[0];
sx q[0];
rz(-2.1213518) q[0];
sx q[0];
rz(0.26225342) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3533808) q[2];
sx q[2];
rz(-1.916421) q[2];
sx q[2];
rz(2.2923922) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2578656) q[1];
sx q[1];
rz(-2.1818433) q[1];
sx q[1];
rz(0.74531598) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4936968) q[3];
sx q[3];
rz(-2.4101083) q[3];
sx q[3];
rz(-2.8836718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38368791) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(-2.596358) q[2];
rz(-0.43045726) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(-2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6298744) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(-1.2114725) q[0];
rz(-2.1684872) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(-2.1957695) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0154008) q[0];
sx q[0];
rz(-1.7915103) q[0];
sx q[0];
rz(-1.6507694) q[0];
x q[1];
rz(-0.81994762) q[2];
sx q[2];
rz(-2.2689399) q[2];
sx q[2];
rz(-1.9422216) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3555803) q[1];
sx q[1];
rz(-1.5994206) q[1];
sx q[1];
rz(-2.7629258) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66611992) q[3];
sx q[3];
rz(-1.5770116) q[3];
sx q[3];
rz(-2.4339649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8699845) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(2.8611709) q[2];
rz(-0.60493207) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(2.159507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1542926) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(-0.61532414) q[0];
rz(-0.92957169) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(0.5756793) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070411365) q[0];
sx q[0];
rz(-1.0072664) q[0];
sx q[0];
rz(2.1348743) q[0];
rz(1.3747146) q[2];
sx q[2];
rz(-1.6819281) q[2];
sx q[2];
rz(-2.9642504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5513735) q[1];
sx q[1];
rz(-3.0312523) q[1];
sx q[1];
rz(1.8054086) q[1];
x q[2];
rz(-2.9541295) q[3];
sx q[3];
rz(-1.0378569) q[3];
sx q[3];
rz(-2.370358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4158674) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(-3.1126853) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(0.96910563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4934851) q[0];
sx q[0];
rz(-3.0457975) q[0];
sx q[0];
rz(-0.84798725) q[0];
rz(-2.3949941) q[2];
sx q[2];
rz(-1.4274297) q[2];
sx q[2];
rz(1.6105086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5574054) q[1];
sx q[1];
rz(-2.4483042) q[1];
sx q[1];
rz(2.8597945) q[1];
x q[2];
rz(-0.5523596) q[3];
sx q[3];
rz(-0.5061572) q[3];
sx q[3];
rz(-1.5137223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6897631) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(2.7815212) q[2];
rz(1.6137971) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(-0.56267363) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2071028) q[0];
sx q[0];
rz(-1.5710545) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(0.044152505) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(-1.953223) q[2];
sx q[2];
rz(-1.8845176) q[2];
sx q[2];
rz(0.11558576) q[2];
rz(-2.6882761) q[3];
sx q[3];
rz(-1.286187) q[3];
sx q[3];
rz(-1.5766531) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];