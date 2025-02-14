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
rz(2.9945381) q[0];
sx q[0];
rz(-1.7400063) q[0];
sx q[0];
rz(-0.92010486) q[0];
rz(-0.080634557) q[1];
sx q[1];
rz(3.7205003) q[1];
sx q[1];
rz(10.401934) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29229627) q[0];
sx q[0];
rz(-1.9417282) q[0];
sx q[0];
rz(-1.7975259) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5844272) q[2];
sx q[2];
rz(-1.4385828) q[2];
sx q[2];
rz(-1.7280098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0339573) q[1];
sx q[1];
rz(-2.1672241) q[1];
sx q[1];
rz(1.4075085) q[1];
rz(-pi) q[2];
rz(2.4529934) q[3];
sx q[3];
rz(-1.8103726) q[3];
sx q[3];
rz(1.857615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5620293) q[2];
sx q[2];
rz(-1.2144438) q[2];
sx q[2];
rz(2.5207632) q[2];
rz(0.51554716) q[3];
sx q[3];
rz(-1.7190869) q[3];
sx q[3];
rz(-2.2990885) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2466549) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(0.57902336) q[0];
rz(0.82194263) q[1];
sx q[1];
rz(-1.6252981) q[1];
sx q[1];
rz(-2.6878405) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4849783) q[0];
sx q[0];
rz(-1.3458283) q[0];
sx q[0];
rz(-2.5091835) q[0];
rz(2.6751509) q[2];
sx q[2];
rz(-1.4701519) q[2];
sx q[2];
rz(-2.538344) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.380341) q[1];
sx q[1];
rz(-0.56219343) q[1];
sx q[1];
rz(-2.9060049) q[1];
x q[2];
rz(2.2367495) q[3];
sx q[3];
rz(-0.65757759) q[3];
sx q[3];
rz(-0.50298728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73432505) q[2];
sx q[2];
rz(-3.0027323) q[2];
sx q[2];
rz(-0.56488758) q[2];
rz(-2.7867553) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(1.423665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1693717) q[0];
sx q[0];
rz(-0.2114978) q[0];
sx q[0];
rz(-0.55150223) q[0];
rz(-2.7768199) q[1];
sx q[1];
rz(-2.8021937) q[1];
sx q[1];
rz(0.50484467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4992139) q[0];
sx q[0];
rz(-1.9758245) q[0];
sx q[0];
rz(-1.3971055) q[0];
rz(0.41270035) q[2];
sx q[2];
rz(-1.6566212) q[2];
sx q[2];
rz(2.7049261) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80124679) q[1];
sx q[1];
rz(-0.97724229) q[1];
sx q[1];
rz(0.98538633) q[1];
rz(-pi) q[2];
rz(-0.64845131) q[3];
sx q[3];
rz(-2.5621427) q[3];
sx q[3];
rz(3.0045829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43368936) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(-3.0943387) q[2];
rz(-2.3047678) q[3];
sx q[3];
rz(-2.4182726) q[3];
sx q[3];
rz(-1.0251454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4243917) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(1.8151872) q[0];
rz(1.4855509) q[1];
sx q[1];
rz(-1.0842208) q[1];
sx q[1];
rz(-0.63492376) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7039873) q[0];
sx q[0];
rz(-0.018089596) q[0];
sx q[0];
rz(1.6419069) q[0];
x q[1];
rz(2.0697753) q[2];
sx q[2];
rz(-1.8868251) q[2];
sx q[2];
rz(-1.8763148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8532456) q[1];
sx q[1];
rz(-1.6430055) q[1];
sx q[1];
rz(1.9068524) q[1];
rz(0.67948273) q[3];
sx q[3];
rz(-2.3906997) q[3];
sx q[3];
rz(2.7487432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2403468) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(1.7690313) q[2];
rz(1.9339804) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(2.8670368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32886252) q[0];
sx q[0];
rz(-2.6565318) q[0];
sx q[0];
rz(3.0364756) q[0];
rz(-0.33991995) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(-1.0629268) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30468291) q[0];
sx q[0];
rz(-1.5567008) q[0];
sx q[0];
rz(2.430116) q[0];
rz(-pi) q[1];
rz(-1.8196443) q[2];
sx q[2];
rz(-2.3169059) q[2];
sx q[2];
rz(0.4363554) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.868456) q[1];
sx q[1];
rz(-1.5850164) q[1];
sx q[1];
rz(-0.058560024) q[1];
x q[2];
rz(-1.5751198) q[3];
sx q[3];
rz(-0.49884847) q[3];
sx q[3];
rz(0.81567848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0792599) q[2];
sx q[2];
rz(-1.7534813) q[2];
sx q[2];
rz(-1.3366535) q[2];
rz(3.0873599) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(-2.5469053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9688251) q[0];
sx q[0];
rz(-0.69197881) q[0];
sx q[0];
rz(1.2139976) q[0];
rz(-1.0460151) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(2.2878343) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3873685) q[0];
sx q[0];
rz(-1.4367391) q[0];
sx q[0];
rz(-0.06019528) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3860903) q[2];
sx q[2];
rz(-2.0167588) q[2];
sx q[2];
rz(1.310629) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5723107) q[1];
sx q[1];
rz(-2.1220868) q[1];
sx q[1];
rz(-2.2612919) q[1];
rz(0.34900174) q[3];
sx q[3];
rz(-1.4274538) q[3];
sx q[3];
rz(-0.81763148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13009109) q[2];
sx q[2];
rz(-1.4902196) q[2];
sx q[2];
rz(2.3940274) q[2];
rz(-2.2085564) q[3];
sx q[3];
rz(-1.7739762) q[3];
sx q[3];
rz(2.8396377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0316684) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(2.5001496) q[0];
rz(-0.96915069) q[1];
sx q[1];
rz(-2.0717924) q[1];
sx q[1];
rz(1.1136805) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9819558) q[0];
sx q[0];
rz(-0.41712077) q[0];
sx q[0];
rz(1.3288767) q[0];
rz(-pi) q[1];
rz(-2.2458057) q[2];
sx q[2];
rz(-0.75715827) q[2];
sx q[2];
rz(2.1799708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8889474) q[1];
sx q[1];
rz(-1.749649) q[1];
sx q[1];
rz(-1.0429383) q[1];
rz(-0.52033958) q[3];
sx q[3];
rz(-1.6005777) q[3];
sx q[3];
rz(-2.7958718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.574719) q[2];
sx q[2];
rz(-0.69872624) q[2];
sx q[2];
rz(0.75801545) q[2];
rz(2.8356683) q[3];
sx q[3];
rz(-2.7829792) q[3];
sx q[3];
rz(0.44736403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7253983) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(2.388227) q[0];
rz(2.2196409) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(3.0070378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050722402) q[0];
sx q[0];
rz(-2.7742552) q[0];
sx q[0];
rz(1.814117) q[0];
x q[1];
rz(1.8250663) q[2];
sx q[2];
rz(-0.94506028) q[2];
sx q[2];
rz(2.3601687) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1052195) q[1];
sx q[1];
rz(-1.4810307) q[1];
sx q[1];
rz(1.5516993) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55338545) q[3];
sx q[3];
rz(-2.3883005) q[3];
sx q[3];
rz(0.91915059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77384633) q[2];
sx q[2];
rz(-1.4464804) q[2];
sx q[2];
rz(-1.9473677) q[2];
rz(-1.1456683) q[3];
sx q[3];
rz(-2.001389) q[3];
sx q[3];
rz(1.0660508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0399748) q[0];
sx q[0];
rz(-2.6823253) q[0];
sx q[0];
rz(0.38247821) q[0];
rz(-0.76599145) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(-1.8006178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95582286) q[0];
sx q[0];
rz(-0.18322028) q[0];
sx q[0];
rz(-2.630156) q[0];
rz(-pi) q[1];
rz(-0.75720301) q[2];
sx q[2];
rz(-1.6348038) q[2];
sx q[2];
rz(-2.7724977) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.080440259) q[1];
sx q[1];
rz(-0.78394475) q[1];
sx q[1];
rz(0.32439167) q[1];
rz(1.789647) q[3];
sx q[3];
rz(-0.36389458) q[3];
sx q[3];
rz(-0.6536676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0535023) q[2];
sx q[2];
rz(-0.92038766) q[2];
sx q[2];
rz(-3.0255393) q[2];
rz(1.2608438) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(1.6837696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.5905404) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(-2.9569448) q[0];
rz(-2.2231936) q[1];
sx q[1];
rz(-1.5946439) q[1];
sx q[1];
rz(1.7041697) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6708095) q[0];
sx q[0];
rz(-0.52007404) q[0];
sx q[0];
rz(1.1519679) q[0];
x q[1];
rz(-0.080731656) q[2];
sx q[2];
rz(-0.7237607) q[2];
sx q[2];
rz(2.2249976) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5823707) q[1];
sx q[1];
rz(-2.5032024) q[1];
sx q[1];
rz(1.8965782) q[1];
x q[2];
rz(-2.4653696) q[3];
sx q[3];
rz(-1.8340602) q[3];
sx q[3];
rz(-0.77419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.664428) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(2.5028382) q[2];
rz(-0.81513682) q[3];
sx q[3];
rz(-0.4709979) q[3];
sx q[3];
rz(-0.42069978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4169793) q[0];
sx q[0];
rz(-1.8920349) q[0];
sx q[0];
rz(-2.98988) q[0];
rz(-2.3607415) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(-1.9736171) q[2];
sx q[2];
rz(-0.30277534) q[2];
sx q[2];
rz(1.4390611) q[2];
rz(-1.609997) q[3];
sx q[3];
rz(-1.8146252) q[3];
sx q[3];
rz(-1.6011325) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
