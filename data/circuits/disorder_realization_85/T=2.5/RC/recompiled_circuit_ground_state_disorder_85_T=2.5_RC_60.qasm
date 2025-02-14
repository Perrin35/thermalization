OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6326132) q[0];
sx q[0];
rz(-2.296663) q[0];
sx q[0];
rz(-0.71339575) q[0];
rz(-0.20819918) q[1];
sx q[1];
rz(-2.9543076) q[1];
sx q[1];
rz(2.0050144) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5644585) q[0];
sx q[0];
rz(-0.73720471) q[0];
sx q[0];
rz(1.2080824) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4955253) q[2];
sx q[2];
rz(-2.4589361) q[2];
sx q[2];
rz(-0.45479506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7311095) q[1];
sx q[1];
rz(-1.4315194) q[1];
sx q[1];
rz(-0.21674214) q[1];
rz(0.30367476) q[3];
sx q[3];
rz(-1.4436385) q[3];
sx q[3];
rz(0.87005471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.95181695) q[2];
sx q[2];
rz(-1.2729898) q[2];
sx q[2];
rz(-0.58958685) q[2];
rz(-0.014852614) q[3];
sx q[3];
rz(-1.868052) q[3];
sx q[3];
rz(2.5428037) q[3];
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
rz(0.0066765229) q[0];
sx q[0];
rz(-2.0966457) q[0];
sx q[0];
rz(-0.45231229) q[0];
rz(0.85330883) q[1];
sx q[1];
rz(-2.6494458) q[1];
sx q[1];
rz(-2.7795627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872386) q[0];
sx q[0];
rz(-1.2487129) q[0];
sx q[0];
rz(-2.5117666) q[0];
rz(-0.30899991) q[2];
sx q[2];
rz(-0.51181343) q[2];
sx q[2];
rz(1.9644033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.39581842) q[1];
sx q[1];
rz(-0.89268301) q[1];
sx q[1];
rz(-0.78435244) q[1];
rz(2.9647847) q[3];
sx q[3];
rz(-0.71030092) q[3];
sx q[3];
rz(-1.7506404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8028458) q[2];
sx q[2];
rz(-2.5669528) q[2];
sx q[2];
rz(0.90633264) q[2];
rz(1.4551506) q[3];
sx q[3];
rz(-0.95821277) q[3];
sx q[3];
rz(1.5549392) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872658) q[0];
sx q[0];
rz(-1.9540906) q[0];
sx q[0];
rz(0.97386709) q[0];
rz(-2.5663238) q[1];
sx q[1];
rz(-1.560805) q[1];
sx q[1];
rz(-1.5633993) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78249507) q[0];
sx q[0];
rz(-1.5486915) q[0];
sx q[0];
rz(1.5074411) q[0];
x q[1];
rz(0.44610666) q[2];
sx q[2];
rz(-2.2488454) q[2];
sx q[2];
rz(-2.7715832) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.089019589) q[1];
sx q[1];
rz(-1.1088437) q[1];
sx q[1];
rz(-2.3506052) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3493912) q[3];
sx q[3];
rz(-0.81238895) q[3];
sx q[3];
rz(1.0426857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8156208) q[2];
sx q[2];
rz(-1.1824563) q[2];
sx q[2];
rz(1.6131442) q[2];
rz(3.1260887) q[3];
sx q[3];
rz(-1.9663234) q[3];
sx q[3];
rz(-1.6541245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7148606) q[0];
sx q[0];
rz(-0.69789129) q[0];
sx q[0];
rz(0.062407169) q[0];
rz(-1.9852253) q[1];
sx q[1];
rz(-0.9461177) q[1];
sx q[1];
rz(-0.83414331) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4164364) q[0];
sx q[0];
rz(-1.2761937) q[0];
sx q[0];
rz(-2.0935154) q[0];
rz(-pi) q[1];
rz(-1.9454003) q[2];
sx q[2];
rz(-3.0060177) q[2];
sx q[2];
rz(-1.6983528) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8319611) q[1];
sx q[1];
rz(-0.97339645) q[1];
sx q[1];
rz(-3.076665) q[1];
rz(-pi) q[2];
rz(0.20872727) q[3];
sx q[3];
rz(-2.0794915) q[3];
sx q[3];
rz(1.824348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9741029) q[2];
sx q[2];
rz(-0.17593273) q[2];
sx q[2];
rz(-1.2220194) q[2];
rz(-0.40056285) q[3];
sx q[3];
rz(-1.0223072) q[3];
sx q[3];
rz(2.9266973) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062531384) q[0];
sx q[0];
rz(-0.65199861) q[0];
sx q[0];
rz(-0.25667152) q[0];
rz(-1.7968862) q[1];
sx q[1];
rz(-1.2213629) q[1];
sx q[1];
rz(1.409097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2891319) q[0];
sx q[0];
rz(-2.2268804) q[0];
sx q[0];
rz(0.38087861) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1878384) q[2];
sx q[2];
rz(-2.8013392) q[2];
sx q[2];
rz(-0.34274064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6761634) q[1];
sx q[1];
rz(-1.6106642) q[1];
sx q[1];
rz(2.7207147) q[1];
x q[2];
rz(-1.872329) q[3];
sx q[3];
rz(-2.216194) q[3];
sx q[3];
rz(1.0848349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1008489) q[2];
sx q[2];
rz(-0.70793968) q[2];
sx q[2];
rz(2.9948998) q[2];
rz(-1.3692859) q[3];
sx q[3];
rz(-2.7761603) q[3];
sx q[3];
rz(2.9036314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9323102) q[0];
sx q[0];
rz(-2.3754061) q[0];
sx q[0];
rz(0.76195088) q[0];
rz(-0.85488287) q[1];
sx q[1];
rz(-1.8652752) q[1];
sx q[1];
rz(0.68861047) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5999187) q[0];
sx q[0];
rz(-1.1038938) q[0];
sx q[0];
rz(-0.093976973) q[0];
rz(1.5499849) q[2];
sx q[2];
rz(-1.6841128) q[2];
sx q[2];
rz(1.3993539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2999867) q[1];
sx q[1];
rz(-3.0191023) q[1];
sx q[1];
rz(-0.56614205) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9290885) q[3];
sx q[3];
rz(-0.52364615) q[3];
sx q[3];
rz(1.3646082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7410572) q[2];
sx q[2];
rz(-2.1964938) q[2];
sx q[2];
rz(-0.98089027) q[2];
rz(0.26997057) q[3];
sx q[3];
rz(-1.910784) q[3];
sx q[3];
rz(0.009416906) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16801676) q[0];
sx q[0];
rz(-1.4485757) q[0];
sx q[0];
rz(0.5434522) q[0];
rz(-1.5914397) q[1];
sx q[1];
rz(-2.2563069) q[1];
sx q[1];
rz(-0.98522225) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5352683) q[0];
sx q[0];
rz(-2.599975) q[0];
sx q[0];
rz(2.6846425) q[0];
rz(-pi) q[1];
rz(2.1357762) q[2];
sx q[2];
rz(-1.8881637) q[2];
sx q[2];
rz(-1.3591426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4812659) q[1];
sx q[1];
rz(-0.72176343) q[1];
sx q[1];
rz(-0.97246171) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6292218) q[3];
sx q[3];
rz(-1.2768043) q[3];
sx q[3];
rz(1.1666573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33124179) q[2];
sx q[2];
rz(-2.6369429) q[2];
sx q[2];
rz(-1.3387559) q[2];
rz(1.4618368) q[3];
sx q[3];
rz(-1.5906518) q[3];
sx q[3];
rz(2.4707879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93002334) q[0];
sx q[0];
rz(-1.1299364) q[0];
sx q[0];
rz(1.9071963) q[0];
rz(-1.3537539) q[1];
sx q[1];
rz(-0.45629558) q[1];
sx q[1];
rz(3.0455468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1559469) q[0];
sx q[0];
rz(-2.0237158) q[0];
sx q[0];
rz(2.3326796) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0158218) q[2];
sx q[2];
rz(-2.5979804) q[2];
sx q[2];
rz(2.8961492) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7754566) q[1];
sx q[1];
rz(-1.6153496) q[1];
sx q[1];
rz(2.7612727) q[1];
rz(-1.6173564) q[3];
sx q[3];
rz(-0.73107204) q[3];
sx q[3];
rz(-0.94207596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8058406) q[2];
sx q[2];
rz(-2.7554607) q[2];
sx q[2];
rz(-1.1248355) q[2];
rz(-2.2624894) q[3];
sx q[3];
rz(-1.9156888) q[3];
sx q[3];
rz(0.44900289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6833471) q[0];
sx q[0];
rz(-1.6793716) q[0];
sx q[0];
rz(-0.71518389) q[0];
rz(2.8257651) q[1];
sx q[1];
rz(-0.68604699) q[1];
sx q[1];
rz(0.16206965) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.049555) q[0];
sx q[0];
rz(-2.3020235) q[0];
sx q[0];
rz(0.58734307) q[0];
x q[1];
rz(0.58312441) q[2];
sx q[2];
rz(-1.7411425) q[2];
sx q[2];
rz(1.2254465) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9980364) q[1];
sx q[1];
rz(-2.2846235) q[1];
sx q[1];
rz(-1.803323) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8211179) q[3];
sx q[3];
rz(-1.2487914) q[3];
sx q[3];
rz(-2.1706949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7519303) q[2];
sx q[2];
rz(-1.1952362) q[2];
sx q[2];
rz(1.2644838) q[2];
rz(-1.3324995) q[3];
sx q[3];
rz(-1.8615362) q[3];
sx q[3];
rz(0.31681028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463117) q[0];
sx q[0];
rz(-2.3667211) q[0];
sx q[0];
rz(-2.0922022) q[0];
rz(0.28114444) q[1];
sx q[1];
rz(-2.099359) q[1];
sx q[1];
rz(-1.3220471) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78681386) q[0];
sx q[0];
rz(-0.6118868) q[0];
sx q[0];
rz(-2.6953681) q[0];
rz(0.63788267) q[2];
sx q[2];
rz(-1.97458) q[2];
sx q[2];
rz(-1.0255726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2019649) q[1];
sx q[1];
rz(-1.2117821) q[1];
sx q[1];
rz(-0.16965349) q[1];
rz(0.23909823) q[3];
sx q[3];
rz(-2.7452668) q[3];
sx q[3];
rz(-2.9843085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7634924) q[2];
sx q[2];
rz(-1.8201479) q[2];
sx q[2];
rz(0.05149252) q[2];
rz(-1.3242877) q[3];
sx q[3];
rz(-1.4824425) q[3];
sx q[3];
rz(-1.9342559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0180203) q[0];
sx q[0];
rz(-1.8089031) q[0];
sx q[0];
rz(-1.7261169) q[0];
rz(0.82072683) q[1];
sx q[1];
rz(-1.4438933) q[1];
sx q[1];
rz(-1.874598) q[1];
rz(-1.8078697) q[2];
sx q[2];
rz(-0.96482926) q[2];
sx q[2];
rz(0.0083487645) q[2];
rz(1.2352712) q[3];
sx q[3];
rz(-1.0809837) q[3];
sx q[3];
rz(-3.0726846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
