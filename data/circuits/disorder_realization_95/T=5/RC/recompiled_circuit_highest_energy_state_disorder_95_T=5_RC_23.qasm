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
rz(1.3746102) q[0];
sx q[0];
rz(4.1144028) q[0];
sx q[0];
rz(10.65539) q[0];
rz(-1.542701) q[1];
sx q[1];
rz(-0.35294947) q[1];
sx q[1];
rz(0.43500873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6957851) q[0];
sx q[0];
rz(-1.2381993) q[0];
sx q[0];
rz(-0.10528721) q[0];
x q[1];
rz(0.22763197) q[2];
sx q[2];
rz(-1.442558) q[2];
sx q[2];
rz(-1.1007835) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5953345) q[1];
sx q[1];
rz(-1.8500684) q[1];
sx q[1];
rz(0.78650773) q[1];
rz(1.5391971) q[3];
sx q[3];
rz(-1.5571619) q[3];
sx q[3];
rz(-2.0054833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1970485) q[2];
sx q[2];
rz(-0.41434449) q[2];
sx q[2];
rz(2.22331) q[2];
rz(2.7783172) q[3];
sx q[3];
rz(-1.8922292) q[3];
sx q[3];
rz(-1.3242807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0927703) q[0];
sx q[0];
rz(-2.4019882) q[0];
sx q[0];
rz(1.1356461) q[0];
rz(-0.28226918) q[1];
sx q[1];
rz(-1.5156563) q[1];
sx q[1];
rz(2.938882) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36798635) q[0];
sx q[0];
rz(-1.752089) q[0];
sx q[0];
rz(1.7530935) q[0];
rz(-pi) q[1];
rz(1.0890929) q[2];
sx q[2];
rz(-1.559245) q[2];
sx q[2];
rz(-2.3318219) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.41676471) q[1];
sx q[1];
rz(-1.4484805) q[1];
sx q[1];
rz(-0.45682795) q[1];
x q[2];
rz(0.0052350076) q[3];
sx q[3];
rz(-1.5605032) q[3];
sx q[3];
rz(2.0106842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2477766) q[2];
sx q[2];
rz(-2.717369) q[2];
sx q[2];
rz(0.068537863) q[2];
rz(2.2872772) q[3];
sx q[3];
rz(-1.8770542) q[3];
sx q[3];
rz(-3.1356649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0900367) q[0];
sx q[0];
rz(-2.5904901) q[0];
sx q[0];
rz(1.0968444) q[0];
rz(-0.54496533) q[1];
sx q[1];
rz(-0.68383354) q[1];
sx q[1];
rz(0.1942689) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30667728) q[0];
sx q[0];
rz(-0.55473548) q[0];
sx q[0];
rz(0.86443211) q[0];
rz(-1.0610007) q[2];
sx q[2];
rz(-0.61697996) q[2];
sx q[2];
rz(-1.7327019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.771505) q[1];
sx q[1];
rz(-1.2125101) q[1];
sx q[1];
rz(2.6111228) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9241289) q[3];
sx q[3];
rz(-1.6969181) q[3];
sx q[3];
rz(2.452891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37476173) q[2];
sx q[2];
rz(-1.7685879) q[2];
sx q[2];
rz(0.72804803) q[2];
rz(-2.9211365) q[3];
sx q[3];
rz(-2.9900592) q[3];
sx q[3];
rz(2.5098586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8927638) q[0];
sx q[0];
rz(-1.4735824) q[0];
sx q[0];
rz(1.0695176) q[0];
rz(-2.2546774) q[1];
sx q[1];
rz(-0.506217) q[1];
sx q[1];
rz(-1.3217529) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66002125) q[0];
sx q[0];
rz(-1.1371326) q[0];
sx q[0];
rz(1.6138484) q[0];
rz(-1.5863877) q[2];
sx q[2];
rz(-1.8097477) q[2];
sx q[2];
rz(0.76214253) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9001213) q[1];
sx q[1];
rz(-1.6423499) q[1];
sx q[1];
rz(2.009719) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2668244) q[3];
sx q[3];
rz(-0.7255377) q[3];
sx q[3];
rz(2.1005132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6679823) q[2];
sx q[2];
rz(-2.0986291) q[2];
sx q[2];
rz(1.854151) q[2];
rz(-2.1227664) q[3];
sx q[3];
rz(-0.5580709) q[3];
sx q[3];
rz(-0.67600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92893112) q[0];
sx q[0];
rz(-2.9086845) q[0];
sx q[0];
rz(0.88053298) q[0];
rz(0.13678837) q[1];
sx q[1];
rz(-0.22577481) q[1];
sx q[1];
rz(-0.61019623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5870415) q[0];
sx q[0];
rz(-0.64233883) q[0];
sx q[0];
rz(1.5917529) q[0];
x q[1];
rz(-1.4177104) q[2];
sx q[2];
rz(-1.1436074) q[2];
sx q[2];
rz(1.9356021) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6366263) q[1];
sx q[1];
rz(-1.7646878) q[1];
sx q[1];
rz(-0.30047471) q[1];
rz(-1.4349157) q[3];
sx q[3];
rz(-1.1349832) q[3];
sx q[3];
rz(-2.6282981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5337164) q[2];
sx q[2];
rz(-0.72489649) q[2];
sx q[2];
rz(-1.9113212) q[2];
rz(2.2599334) q[3];
sx q[3];
rz(-1.6917546) q[3];
sx q[3];
rz(-1.3151883) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45516685) q[0];
sx q[0];
rz(-1.5979586) q[0];
sx q[0];
rz(2.0879188) q[0];
rz(-1.3764489) q[1];
sx q[1];
rz(-1.1076628) q[1];
sx q[1];
rz(-2.460316) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85515012) q[0];
sx q[0];
rz(-1.5871829) q[0];
sx q[0];
rz(-0.70006242) q[0];
rz(-pi) q[1];
rz(-2.3967016) q[2];
sx q[2];
rz(-1.3110263) q[2];
sx q[2];
rz(-1.7967173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14398814) q[1];
sx q[1];
rz(-1.0790842) q[1];
sx q[1];
rz(-2.8365447) q[1];
rz(1.525649) q[3];
sx q[3];
rz(-2.4490926) q[3];
sx q[3];
rz(0.11603234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.072558746) q[2];
sx q[2];
rz(-2.7636038) q[2];
sx q[2];
rz(0.47522137) q[2];
rz(-1.9184939) q[3];
sx q[3];
rz(-0.94448543) q[3];
sx q[3];
rz(2.9782915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2326736) q[0];
sx q[0];
rz(-2.1894046) q[0];
sx q[0];
rz(-2.9222144) q[0];
rz(-1.0064005) q[1];
sx q[1];
rz(-1.718113) q[1];
sx q[1];
rz(1.2054319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4426098) q[0];
sx q[0];
rz(-2.1694393) q[0];
sx q[0];
rz(-0.026431008) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57997616) q[2];
sx q[2];
rz(-2.6784228) q[2];
sx q[2];
rz(1.5214187) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6483375) q[1];
sx q[1];
rz(-1.770079) q[1];
sx q[1];
rz(-1.553345) q[1];
x q[2];
rz(2.0153644) q[3];
sx q[3];
rz(-1.6939112) q[3];
sx q[3];
rz(-3.0266704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.22126234) q[2];
sx q[2];
rz(-2.7483676) q[2];
sx q[2];
rz(2.4779251) q[2];
rz(0.71781939) q[3];
sx q[3];
rz(-0.65687537) q[3];
sx q[3];
rz(2.9653911) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8127036) q[0];
sx q[0];
rz(-0.95903522) q[0];
sx q[0];
rz(-1.2411728) q[0];
rz(-2.7136956) q[1];
sx q[1];
rz(-1.5623845) q[1];
sx q[1];
rz(1.4235628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9374838) q[0];
sx q[0];
rz(-2.9108469) q[0];
sx q[0];
rz(0.33270128) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85278089) q[2];
sx q[2];
rz(-1.4508392) q[2];
sx q[2];
rz(-0.44242417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1686624) q[1];
sx q[1];
rz(-2.2386547) q[1];
sx q[1];
rz(-3.1144823) q[1];
rz(-pi) q[2];
rz(-3.0377722) q[3];
sx q[3];
rz(-0.33664931) q[3];
sx q[3];
rz(0.81538425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9972035) q[2];
sx q[2];
rz(-2.453697) q[2];
sx q[2];
rz(2.8456861) q[2];
rz(-1.0812673) q[3];
sx q[3];
rz(-1.6176977) q[3];
sx q[3];
rz(0.14779873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57701552) q[0];
sx q[0];
rz(-2.5595589) q[0];
sx q[0];
rz(-2.4356595) q[0];
rz(2.9544746) q[1];
sx q[1];
rz(-0.83693224) q[1];
sx q[1];
rz(-0.44609889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69639403) q[0];
sx q[0];
rz(-0.15226752) q[0];
sx q[0];
rz(-0.012114004) q[0];
rz(-2.1438333) q[2];
sx q[2];
rz(-2.5964542) q[2];
sx q[2];
rz(-2.0968495) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6793788) q[1];
sx q[1];
rz(-2.129038) q[1];
sx q[1];
rz(-2.7427196) q[1];
rz(2.5650065) q[3];
sx q[3];
rz(-1.4268677) q[3];
sx q[3];
rz(-1.3670539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17911653) q[2];
sx q[2];
rz(-1.0255145) q[2];
sx q[2];
rz(1.7402488) q[2];
rz(-2.5380747) q[3];
sx q[3];
rz(-0.18134914) q[3];
sx q[3];
rz(-0.13874273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(3.1097581) q[0];
sx q[0];
rz(-1.4048445) q[0];
sx q[0];
rz(-0.11012125) q[0];
rz(-1.8203863) q[1];
sx q[1];
rz(-0.92480129) q[1];
sx q[1];
rz(-0.22470156) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49602717) q[0];
sx q[0];
rz(-1.1764045) q[0];
sx q[0];
rz(2.140122) q[0];
rz(-pi) q[1];
rz(0.16188936) q[2];
sx q[2];
rz(-1.8542552) q[2];
sx q[2];
rz(0.12609158) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1096829) q[1];
sx q[1];
rz(-2.2930298) q[1];
sx q[1];
rz(-1.7914045) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7240771) q[3];
sx q[3];
rz(-0.17440052) q[3];
sx q[3];
rz(-1.0518975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8297742) q[2];
sx q[2];
rz(-0.31121397) q[2];
sx q[2];
rz(-2.6149926) q[2];
rz(2.7569568) q[3];
sx q[3];
rz(-1.2472109) q[3];
sx q[3];
rz(2.9068936) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93152355) q[0];
sx q[0];
rz(-0.43545224) q[0];
sx q[0];
rz(-0.42238105) q[0];
rz(0.79413636) q[1];
sx q[1];
rz(-1.4310478) q[1];
sx q[1];
rz(-1.4337883) q[1];
rz(2.4985102) q[2];
sx q[2];
rz(-2.1529635) q[2];
sx q[2];
rz(-2.8827428) q[2];
rz(-0.724283) q[3];
sx q[3];
rz(-0.52762994) q[3];
sx q[3];
rz(2.4366196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
