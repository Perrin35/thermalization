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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46086754) q[0];
sx q[0];
rz(-1.4493677) q[0];
sx q[0];
rz(0.073343883) q[0];
rz(-1.471465) q[2];
sx q[2];
rz(-0.28495312) q[2];
sx q[2];
rz(3.1379267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9610112) q[1];
sx q[1];
rz(-2.2716224) q[1];
sx q[1];
rz(0.58971528) q[1];
rz(-pi) q[2];
rz(0.41892003) q[3];
sx q[3];
rz(-1.6994611) q[3];
sx q[3];
rz(-0.045623771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1315786) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(-2.4689891) q[2];
rz(0.16942313) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-0.23068962) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(-0.22856523) q[0];
rz(-0.16054343) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(2.8536318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1624958) q[0];
sx q[0];
rz(-1.393035) q[0];
sx q[0];
rz(-1.7642679) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1177883) q[2];
sx q[2];
rz(-2.4734801) q[2];
sx q[2];
rz(-0.39272768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5353363) q[1];
sx q[1];
rz(-1.1621981) q[1];
sx q[1];
rz(0.6596843) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1045707) q[3];
sx q[3];
rz(-2.7894756) q[3];
sx q[3];
rz(-2.5196688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59445375) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(-1.6681558) q[2];
rz(-2.2041221) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8163452) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(-1.422241) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(-2.1898988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5121582) q[0];
sx q[0];
rz(-2.929147) q[0];
sx q[0];
rz(1.1266707) q[0];
x q[1];
rz(1.0551664) q[2];
sx q[2];
rz(-1.3736758) q[2];
sx q[2];
rz(0.65161639) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.051085) q[1];
sx q[1];
rz(-1.2625492) q[1];
sx q[1];
rz(-0.80866637) q[1];
x q[2];
rz(0.94727256) q[3];
sx q[3];
rz(-1.1215278) q[3];
sx q[3];
rz(-0.63637966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13016985) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(2.7974131) q[2];
rz(-0.88642818) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-2.5755431) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043871) q[0];
sx q[0];
rz(-1.4973649) q[0];
sx q[0];
rz(1.244506) q[0];
rz(3.0124774) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(-0.37277645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9739649) q[0];
sx q[0];
rz(-0.81866696) q[0];
sx q[0];
rz(-1.4139834) q[0];
x q[1];
rz(0.77809019) q[2];
sx q[2];
rz(-1.8847244) q[2];
sx q[2];
rz(-0.84184605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4200538) q[1];
sx q[1];
rz(-2.5905847) q[1];
sx q[1];
rz(1.0934456) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0471441) q[3];
sx q[3];
rz(-1.2762478) q[3];
sx q[3];
rz(-1.095872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.63561511) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(-0.90488952) q[2];
rz(0.44057009) q[3];
sx q[3];
rz(-0.4959271) q[3];
sx q[3];
rz(2.1499965) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6624517) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(1.2878081) q[0];
rz(-1.6632535) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(-0.53422654) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6144845) q[0];
sx q[0];
rz(-1.8728349) q[0];
sx q[0];
rz(1.2472279) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.062837) q[2];
sx q[2];
rz(-1.2675708) q[2];
sx q[2];
rz(-0.505503) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5200107) q[1];
sx q[1];
rz(-0.098712155) q[1];
sx q[1];
rz(2.2732544) q[1];
x q[2];
rz(0.15054536) q[3];
sx q[3];
rz(-2.4171722) q[3];
sx q[3];
rz(-0.62643334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(-2.9980998) q[2];
rz(1.7701373) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(-0.21970704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.4390398) q[0];
sx q[0];
rz(-2.2655903) q[0];
sx q[0];
rz(0.23705661) q[0];
rz(1.2409695) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(-1.7664849) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.188376) q[0];
sx q[0];
rz(-1.9089451) q[0];
sx q[0];
rz(1.2921278) q[0];
rz(1.6429971) q[2];
sx q[2];
rz(-0.82319665) q[2];
sx q[2];
rz(2.2526134) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82735862) q[1];
sx q[1];
rz(-1.5463366) q[1];
sx q[1];
rz(-0.41008653) q[1];
x q[2];
rz(1.9938019) q[3];
sx q[3];
rz(-1.1173964) q[3];
sx q[3];
rz(2.6231433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.53283006) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(2.9928845) q[2];
rz(3.1249629) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(2.9490411) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49611133) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(-0.87316978) q[0];
rz(0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(-3.057664) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5556363) q[0];
sx q[0];
rz(-1.3480098) q[0];
sx q[0];
rz(2.1369834) q[0];
x q[1];
rz(-1.7882118) q[2];
sx q[2];
rz(-1.2251717) q[2];
sx q[2];
rz(-0.84920041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1740239) q[1];
sx q[1];
rz(-2.1597383) q[1];
sx q[1];
rz(0.80935652) q[1];
rz(-1.4936968) q[3];
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
rz(0.38368791) q[2];
sx q[2];
rz(-2.7010475) q[2];
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
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(-1.2114725) q[0];
rz(0.97310549) q[1];
sx q[1];
rz(-0.5935697) q[1];
sx q[1];
rz(2.1957695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1261919) q[0];
sx q[0];
rz(-1.7915103) q[0];
sx q[0];
rz(1.6507694) q[0];
x q[1];
rz(-2.2875167) q[2];
sx q[2];
rz(-1.0208566) q[2];
sx q[2];
rz(-2.9727109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1434053) q[1];
sx q[1];
rz(-2.7618976) q[1];
sx q[1];
rz(3.0642964) q[1];
rz(0.010057851) q[3];
sx q[3];
rz(-0.66614449) q[3];
sx q[3];
rz(-2.2705164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(2.8611709) q[2];
rz(0.60493207) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(-2.159507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.1542926) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(2.5262685) q[0];
rz(-0.92957169) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(0.5756793) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0711813) q[0];
sx q[0];
rz(-2.1343263) q[0];
sx q[0];
rz(1.0067183) q[0];
rz(-1.766878) q[2];
sx q[2];
rz(-1.4596645) q[2];
sx q[2];
rz(-0.17734222) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5513735) q[1];
sx q[1];
rz(-0.11034036) q[1];
sx q[1];
rz(1.8054086) q[1];
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
rz(-pi) q[1];
rz(0.82970396) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(-0.021961948) q[2];
rz(-2.9711376) q[3];
sx q[3];
rz(-1.0237834) q[3];
sx q[3];
rz(-0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4158674) q[0];
sx q[0];
rz(-0.9265582) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(3.1126853) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(2.172487) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983738) q[0];
sx q[0];
rz(-1.6341097) q[0];
sx q[0];
rz(-1.6427342) q[0];
rz(-pi) q[1];
rz(1.3766039) q[2];
sx q[2];
rz(-0.83364928) q[2];
sx q[2];
rz(-0.091723524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.58418729) q[1];
sx q[1];
rz(-2.4483042) q[1];
sx q[1];
rz(0.28179817) q[1];
rz(-pi) q[2];
rz(0.44090791) q[3];
sx q[3];
rz(-1.8280067) q[3];
sx q[3];
rz(0.5514901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6897631) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(0.36007145) q[2];
rz(-1.5277956) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(0.56267363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9344899) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(-3.0974401) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(1.1883696) q[2];
sx q[2];
rz(-1.8845176) q[2];
sx q[2];
rz(0.11558576) q[2];
rz(-1.2561856) q[3];
sx q[3];
rz(-2.0046069) q[3];
sx q[3];
rz(-0.14179695) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
