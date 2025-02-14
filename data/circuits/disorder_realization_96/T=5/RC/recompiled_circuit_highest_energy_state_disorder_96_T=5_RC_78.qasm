OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.30492914) q[0];
sx q[0];
rz(-0.066216901) q[0];
sx q[0];
rz(-2.9856292) q[0];
rz(1.1671542) q[1];
sx q[1];
rz(-2.4894297) q[1];
sx q[1];
rz(-0.48643938) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4733361) q[0];
sx q[0];
rz(-2.121068) q[0];
sx q[0];
rz(0.86007287) q[0];
rz(2.0912284) q[2];
sx q[2];
rz(-1.2005521) q[2];
sx q[2];
rz(0.56132078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6156494) q[1];
sx q[1];
rz(-1.2541703) q[1];
sx q[1];
rz(1.8773882) q[1];
rz(2.9286195) q[3];
sx q[3];
rz(-0.83847441) q[3];
sx q[3];
rz(-0.71668437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.277694) q[2];
sx q[2];
rz(-0.70964491) q[2];
sx q[2];
rz(-0.42897439) q[2];
rz(3.0947558) q[3];
sx q[3];
rz(-2.7497079) q[3];
sx q[3];
rz(2.528791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.866339) q[0];
sx q[0];
rz(-2.1493122) q[0];
sx q[0];
rz(-2.476165) q[0];
rz(-1.1412507) q[1];
sx q[1];
rz(-1.3803866) q[1];
sx q[1];
rz(0.64249396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35303799) q[0];
sx q[0];
rz(-0.73563535) q[0];
sx q[0];
rz(3.0658196) q[0];
rz(-pi) q[1];
rz(-2.7879509) q[2];
sx q[2];
rz(-2.414915) q[2];
sx q[2];
rz(0.35219372) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.68243945) q[1];
sx q[1];
rz(-2.2687468) q[1];
sx q[1];
rz(-1.5345011) q[1];
x q[2];
rz(1.2824545) q[3];
sx q[3];
rz(-2.5817462) q[3];
sx q[3];
rz(2.0324872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3147754) q[2];
sx q[2];
rz(-2.3544957) q[2];
sx q[2];
rz(-2.5891506) q[2];
rz(-1.6636482) q[3];
sx q[3];
rz(-1.2976357) q[3];
sx q[3];
rz(0.67599952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8058572) q[0];
sx q[0];
rz(-2.4215846) q[0];
sx q[0];
rz(2.3190401) q[0];
rz(-2.3303253) q[1];
sx q[1];
rz(-0.30458105) q[1];
sx q[1];
rz(1.4623581) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050348076) q[0];
sx q[0];
rz(-1.1093265) q[0];
sx q[0];
rz(1.9965003) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41303582) q[2];
sx q[2];
rz(-1.6390642) q[2];
sx q[2];
rz(0.78536805) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84356704) q[1];
sx q[1];
rz(-2.0223534) q[1];
sx q[1];
rz(2.8140912) q[1];
rz(-pi) q[2];
rz(0.30841804) q[3];
sx q[3];
rz(-1.2756516) q[3];
sx q[3];
rz(-2.4663188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2751969) q[2];
sx q[2];
rz(-2.3123645) q[2];
sx q[2];
rz(2.5466476) q[2];
rz(-2.7368937) q[3];
sx q[3];
rz(-1.1070822) q[3];
sx q[3];
rz(1.8428724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054656) q[0];
sx q[0];
rz(-1.8647702) q[0];
sx q[0];
rz(2.8986616) q[0];
rz(0.88492197) q[1];
sx q[1];
rz(-0.72598571) q[1];
sx q[1];
rz(-3.1387709) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3784175) q[0];
sx q[0];
rz(-0.88930128) q[0];
sx q[0];
rz(-0.50294089) q[0];
rz(0.94945939) q[2];
sx q[2];
rz(-1.1295801) q[2];
sx q[2];
rz(0.37829933) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.62880781) q[1];
sx q[1];
rz(-1.50042) q[1];
sx q[1];
rz(2.2712399) q[1];
rz(-pi) q[2];
rz(-2.5266227) q[3];
sx q[3];
rz(-1.0677862) q[3];
sx q[3];
rz(0.78395432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5367077) q[2];
sx q[2];
rz(-1.1515836) q[2];
sx q[2];
rz(-0.73337698) q[2];
rz(0.65507656) q[3];
sx q[3];
rz(-0.30787444) q[3];
sx q[3];
rz(-0.56183279) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68375278) q[0];
sx q[0];
rz(-2.8003052) q[0];
sx q[0];
rz(0.14232464) q[0];
rz(1.3617474) q[1];
sx q[1];
rz(-0.35265499) q[1];
sx q[1];
rz(2.6432162) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55492586) q[0];
sx q[0];
rz(-2.2827671) q[0];
sx q[0];
rz(0.081454309) q[0];
rz(-3.0875838) q[2];
sx q[2];
rz(-2.3819792) q[2];
sx q[2];
rz(-2.1646966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.315334) q[1];
sx q[1];
rz(-0.79763258) q[1];
sx q[1];
rz(-0.41528292) q[1];
rz(-pi) q[2];
rz(-2.3617414) q[3];
sx q[3];
rz(-1.3289467) q[3];
sx q[3];
rz(2.4132501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.132823) q[2];
sx q[2];
rz(-1.2558197) q[2];
sx q[2];
rz(-0.69457561) q[2];
rz(0.83235598) q[3];
sx q[3];
rz(-2.0280301) q[3];
sx q[3];
rz(-2.8913403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4271127) q[0];
sx q[0];
rz(-0.50943333) q[0];
sx q[0];
rz(0.66202778) q[0];
rz(-1.9561249) q[1];
sx q[1];
rz(-2.1123501) q[1];
sx q[1];
rz(-1.9756636) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2625354) q[0];
sx q[0];
rz(-1.3692807) q[0];
sx q[0];
rz(0.67964566) q[0];
rz(-pi) q[1];
rz(-1.7244965) q[2];
sx q[2];
rz(-2.1734093) q[2];
sx q[2];
rz(0.034687925) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5422029) q[1];
sx q[1];
rz(-1.6205317) q[1];
sx q[1];
rz(-1.6433783) q[1];
rz(-pi) q[2];
rz(0.77591578) q[3];
sx q[3];
rz(-2.3563623) q[3];
sx q[3];
rz(1.531158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.379443) q[2];
sx q[2];
rz(-1.9337485) q[2];
sx q[2];
rz(0.70972788) q[2];
rz(-0.48192561) q[3];
sx q[3];
rz(-0.47526264) q[3];
sx q[3];
rz(-3.1221534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703281) q[0];
sx q[0];
rz(-2.1605991) q[0];
sx q[0];
rz(1.574466) q[0];
rz(1.6471242) q[1];
sx q[1];
rz(-0.44691214) q[1];
sx q[1];
rz(2.2692197) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0013638) q[0];
sx q[0];
rz(-1.2622381) q[0];
sx q[0];
rz(2.9635327) q[0];
rz(-0.72231897) q[2];
sx q[2];
rz(-1.4294942) q[2];
sx q[2];
rz(-1.6091122) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1885729) q[1];
sx q[1];
rz(-1.6011878) q[1];
sx q[1];
rz(-1.4646444) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6029583) q[3];
sx q[3];
rz(-2.4866606) q[3];
sx q[3];
rz(0.82603588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3288154) q[2];
sx q[2];
rz(-2.2810563) q[2];
sx q[2];
rz(-1.1278197) q[2];
rz(2.7382216) q[3];
sx q[3];
rz(-1.6962467) q[3];
sx q[3];
rz(-2.4283714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23680747) q[0];
sx q[0];
rz(-1.9984364) q[0];
sx q[0];
rz(1.2971725) q[0];
rz(-2.6203268) q[1];
sx q[1];
rz(-1.2626941) q[1];
sx q[1];
rz(-3.0788132) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32048827) q[0];
sx q[0];
rz(-2.0135499) q[0];
sx q[0];
rz(0.073320079) q[0];
x q[1];
rz(-0.56030484) q[2];
sx q[2];
rz(-2.7231999) q[2];
sx q[2];
rz(0.095730893) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.081689) q[1];
sx q[1];
rz(-1.4270596) q[1];
sx q[1];
rz(1.3215076) q[1];
rz(-pi) q[2];
rz(2.217123) q[3];
sx q[3];
rz(-1.1271584) q[3];
sx q[3];
rz(0.2547338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2433743) q[2];
sx q[2];
rz(-0.90460193) q[2];
sx q[2];
rz(-2.173219) q[2];
rz(2.6280256) q[3];
sx q[3];
rz(-1.3526724) q[3];
sx q[3];
rz(-2.8106522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6043855) q[0];
sx q[0];
rz(-0.66050118) q[0];
sx q[0];
rz(0.78395098) q[0];
rz(1.9369269) q[1];
sx q[1];
rz(-2.7684863) q[1];
sx q[1];
rz(1.7663667) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7989144) q[0];
sx q[0];
rz(-1.3564137) q[0];
sx q[0];
rz(0.27002295) q[0];
rz(-pi) q[1];
rz(1.398574) q[2];
sx q[2];
rz(-2.9841514) q[2];
sx q[2];
rz(-0.96913183) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32574496) q[1];
sx q[1];
rz(-2.1195852) q[1];
sx q[1];
rz(2.4239848) q[1];
rz(-pi) q[2];
rz(1.2581539) q[3];
sx q[3];
rz(-1.4015159) q[3];
sx q[3];
rz(1.1723684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2969926) q[2];
sx q[2];
rz(-2.832909) q[2];
sx q[2];
rz(-2.6958579) q[2];
rz(0.60278696) q[3];
sx q[3];
rz(-1.876839) q[3];
sx q[3];
rz(0.84356892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8730901) q[0];
sx q[0];
rz(-2.6188445) q[0];
sx q[0];
rz(-0.52421808) q[0];
rz(1.2007319) q[1];
sx q[1];
rz(-1.8914765) q[1];
sx q[1];
rz(-1.9573617) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5843274) q[0];
sx q[0];
rz(-1.3818437) q[0];
sx q[0];
rz(2.7897631) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8866059) q[2];
sx q[2];
rz(-1.6355875) q[2];
sx q[2];
rz(-3.0805317) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47329119) q[1];
sx q[1];
rz(-1.4379825) q[1];
sx q[1];
rz(1.5584438) q[1];
rz(0.6289102) q[3];
sx q[3];
rz(-1.5987608) q[3];
sx q[3];
rz(-0.6599676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.21620096) q[2];
sx q[2];
rz(-0.4500469) q[2];
sx q[2];
rz(2.0897384) q[2];
rz(2.9205186) q[3];
sx q[3];
rz(-1.8408006) q[3];
sx q[3];
rz(-2.2033447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5939519) q[0];
sx q[0];
rz(-1.6407536) q[0];
sx q[0];
rz(-3.092691) q[0];
rz(-0.37638695) q[1];
sx q[1];
rz(-2.2793437) q[1];
sx q[1];
rz(1.9752165) q[1];
rz(-1.6025887) q[2];
sx q[2];
rz(-1.0009264) q[2];
sx q[2];
rz(1.4957503) q[2];
rz(0.54125124) q[3];
sx q[3];
rz(-0.79081906) q[3];
sx q[3];
rz(-2.9959903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
