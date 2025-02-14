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
rz(-1.9953097) q[0];
sx q[0];
rz(-0.028385552) q[0];
sx q[0];
rz(-2.1838768) q[0];
rz(0.43342844) q[1];
sx q[1];
rz(-1.537701) q[1];
sx q[1];
rz(2.5133207) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0032438) q[0];
sx q[0];
rz(-1.4491557) q[0];
sx q[0];
rz(0.37183372) q[0];
rz(-pi) q[1];
rz(1.0721016) q[2];
sx q[2];
rz(-1.9737502) q[2];
sx q[2];
rz(0.13355959) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1294259) q[1];
sx q[1];
rz(-2.2291227) q[1];
sx q[1];
rz(0.40947394) q[1];
x q[2];
rz(3.0423156) q[3];
sx q[3];
rz(-1.9813207) q[3];
sx q[3];
rz(-1.6308189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4665224) q[2];
sx q[2];
rz(-1.0202946) q[2];
sx q[2];
rz(2.7543219) q[2];
rz(-1.9668503) q[3];
sx q[3];
rz(-1.4459123) q[3];
sx q[3];
rz(-2.2430879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3627477) q[0];
sx q[0];
rz(-1.400482) q[0];
sx q[0];
rz(-1.8722906) q[0];
rz(1.6954039) q[1];
sx q[1];
rz(-1.8480999) q[1];
sx q[1];
rz(2.2206025) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22494754) q[0];
sx q[0];
rz(-0.9592255) q[0];
sx q[0];
rz(2.7782604) q[0];
rz(1.3878421) q[2];
sx q[2];
rz(-0.91380807) q[2];
sx q[2];
rz(-0.25568257) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2650406) q[1];
sx q[1];
rz(-1.7608133) q[1];
sx q[1];
rz(-1.9273619) q[1];
x q[2];
rz(-0.030618592) q[3];
sx q[3];
rz(-2.9155882) q[3];
sx q[3];
rz(-1.5196368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9903119) q[2];
sx q[2];
rz(-1.9419779) q[2];
sx q[2];
rz(-2.3178103) q[2];
rz(-1.6173877) q[3];
sx q[3];
rz(-1.5882322) q[3];
sx q[3];
rz(-1.2113781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23564944) q[0];
sx q[0];
rz(-2.2610569) q[0];
sx q[0];
rz(-2.549951) q[0];
rz(0.94332424) q[1];
sx q[1];
rz(-2.0843518) q[1];
sx q[1];
rz(-2.1181769) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34020326) q[0];
sx q[0];
rz(-1.4638299) q[0];
sx q[0];
rz(1.4530327) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1259427) q[2];
sx q[2];
rz(-1.890939) q[2];
sx q[2];
rz(-2.5261836) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0630546) q[1];
sx q[1];
rz(-1.5397433) q[1];
sx q[1];
rz(-1.478315) q[1];
rz(-0.39822762) q[3];
sx q[3];
rz(-2.4169528) q[3];
sx q[3];
rz(-0.46996024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25980514) q[2];
sx q[2];
rz(-1.9116348) q[2];
sx q[2];
rz(-1.7477431) q[2];
rz(-1.4272089) q[3];
sx q[3];
rz(-2.2673159) q[3];
sx q[3];
rz(2.8425596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77089906) q[0];
sx q[0];
rz(-2.735266) q[0];
sx q[0];
rz(2.4942177) q[0];
rz(-2.8992843) q[1];
sx q[1];
rz(-1.7055885) q[1];
sx q[1];
rz(-1.6829596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577135) q[0];
sx q[0];
rz(-1.4365968) q[0];
sx q[0];
rz(-2.8154897) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5299923) q[2];
sx q[2];
rz(-1.0047874) q[2];
sx q[2];
rz(-1.6141487) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9414166) q[1];
sx q[1];
rz(-1.6072909) q[1];
sx q[1];
rz(2.0564276) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69497739) q[3];
sx q[3];
rz(-2.2314592) q[3];
sx q[3];
rz(1.7373178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.709939) q[2];
sx q[2];
rz(-2.2971051) q[2];
sx q[2];
rz(-1.9169774) q[2];
rz(2.8905408) q[3];
sx q[3];
rz(-2.4243088) q[3];
sx q[3];
rz(-0.93799463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5994023) q[0];
sx q[0];
rz(-1.1325862) q[0];
sx q[0];
rz(0.19790025) q[0];
rz(1.9056162) q[1];
sx q[1];
rz(-2.7695152) q[1];
sx q[1];
rz(0.0045675357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2777956) q[0];
sx q[0];
rz(-0.63710538) q[0];
sx q[0];
rz(2.7180014) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.058706) q[2];
sx q[2];
rz(-1.0382639) q[2];
sx q[2];
rz(-0.96585912) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98036042) q[1];
sx q[1];
rz(-1.8440108) q[1];
sx q[1];
rz(1.7262162) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7813579) q[3];
sx q[3];
rz(-1.2270842) q[3];
sx q[3];
rz(-0.76219073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68134754) q[2];
sx q[2];
rz(-0.494151) q[2];
sx q[2];
rz(0.34026185) q[2];
rz(-1.5334689) q[3];
sx q[3];
rz(-1.7483277) q[3];
sx q[3];
rz(-1.5197915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.915864) q[0];
sx q[0];
rz(-0.69750834) q[0];
sx q[0];
rz(-1.9541784) q[0];
rz(1.4215218) q[1];
sx q[1];
rz(-0.41821304) q[1];
sx q[1];
rz(-0.13030599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7447259) q[0];
sx q[0];
rz(-1.5123741) q[0];
sx q[0];
rz(-1.6944042) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2537986) q[2];
sx q[2];
rz(-1.0947795) q[2];
sx q[2];
rz(-2.7107875) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7728945) q[1];
sx q[1];
rz(-0.16783585) q[1];
sx q[1];
rz(-0.0040199587) q[1];
x q[2];
rz(3.1316787) q[3];
sx q[3];
rz(-1.8335206) q[3];
sx q[3];
rz(2.3309121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3336031) q[2];
sx q[2];
rz(-1.5895867) q[2];
sx q[2];
rz(2.1017334) q[2];
rz(-0.17397675) q[3];
sx q[3];
rz(-2.0128553) q[3];
sx q[3];
rz(2.4719293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.533605) q[0];
sx q[0];
rz(-1.0796115) q[0];
sx q[0];
rz(-2.3002891) q[0];
rz(-1.816642) q[1];
sx q[1];
rz(-0.94580301) q[1];
sx q[1];
rz(-1.4680877) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6852035) q[0];
sx q[0];
rz(-2.6273871) q[0];
sx q[0];
rz(-0.77001621) q[0];
rz(-pi) q[1];
rz(1.5505474) q[2];
sx q[2];
rz(-0.51873461) q[2];
sx q[2];
rz(0.27257365) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74875301) q[1];
sx q[1];
rz(-0.76828271) q[1];
sx q[1];
rz(2.4359951) q[1];
x q[2];
rz(0.87196799) q[3];
sx q[3];
rz(-2.6879394) q[3];
sx q[3];
rz(3.0986971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2738721) q[2];
sx q[2];
rz(-0.90110675) q[2];
sx q[2];
rz(-0.64485288) q[2];
rz(2.0598748) q[3];
sx q[3];
rz(-0.41926256) q[3];
sx q[3];
rz(-2.4206415) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4475187) q[0];
sx q[0];
rz(-2.9677291) q[0];
sx q[0];
rz(2.7287667) q[0];
rz(-1.6089571) q[1];
sx q[1];
rz(-0.91333476) q[1];
sx q[1];
rz(2.434381) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0577257) q[0];
sx q[0];
rz(-0.8208771) q[0];
sx q[0];
rz(2.7146401) q[0];
rz(1.7272495) q[2];
sx q[2];
rz(-2.3067637) q[2];
sx q[2];
rz(3.0184658) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8726996) q[1];
sx q[1];
rz(-1.9063236) q[1];
sx q[1];
rz(-0.56443946) q[1];
x q[2];
rz(1.8052638) q[3];
sx q[3];
rz(-2.5904561) q[3];
sx q[3];
rz(2.7296327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1327208) q[2];
sx q[2];
rz(-2.8421695) q[2];
sx q[2];
rz(0.51187619) q[2];
rz(-0.68073186) q[3];
sx q[3];
rz(-1.3635037) q[3];
sx q[3];
rz(2.9752922) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39636382) q[0];
sx q[0];
rz(-1.5928716) q[0];
sx q[0];
rz(-0.45528278) q[0];
rz(1.0633172) q[1];
sx q[1];
rz(-2.2472007) q[1];
sx q[1];
rz(3.0696226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6980879) q[0];
sx q[0];
rz(-0.035497276) q[0];
sx q[0];
rz(-0.3548236) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5939868) q[2];
sx q[2];
rz(-1.1299777) q[2];
sx q[2];
rz(-3.0435516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.2769952) q[1];
sx q[1];
rz(-2.2114843) q[1];
sx q[1];
rz(2.6990141) q[1];
x q[2];
rz(-0.25983475) q[3];
sx q[3];
rz(-1.8913811) q[3];
sx q[3];
rz(-1.1449006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7887743) q[2];
sx q[2];
rz(-0.41663751) q[2];
sx q[2];
rz(-2.519506) q[2];
rz(-2.3173053) q[3];
sx q[3];
rz(-1.0930748) q[3];
sx q[3];
rz(0.85339671) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033009919) q[0];
sx q[0];
rz(-2.0350631) q[0];
sx q[0];
rz(-0.95091096) q[0];
rz(-1.7225601) q[1];
sx q[1];
rz(-2.6585237) q[1];
sx q[1];
rz(-0.61666617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71910673) q[0];
sx q[0];
rz(-0.93084836) q[0];
sx q[0];
rz(-1.6011222) q[0];
rz(1.9497112) q[2];
sx q[2];
rz(-2.1466563) q[2];
sx q[2];
rz(2.4873231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.93393713) q[1];
sx q[1];
rz(-2.9300234) q[1];
sx q[1];
rz(3.1283698) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70722001) q[3];
sx q[3];
rz(-1.116718) q[3];
sx q[3];
rz(1.2421158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0509433) q[2];
sx q[2];
rz(-1.9885149) q[2];
sx q[2];
rz(-2.0743745) q[2];
rz(-2.0459335) q[3];
sx q[3];
rz(-1.2376384) q[3];
sx q[3];
rz(-2.1777976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6296366) q[0];
sx q[0];
rz(-0.98573276) q[0];
sx q[0];
rz(-1.0712256) q[0];
rz(-1.4818954) q[1];
sx q[1];
rz(-2.0493458) q[1];
sx q[1];
rz(0.58072166) q[1];
rz(-1.5956249) q[2];
sx q[2];
rz(-1.8332743) q[2];
sx q[2];
rz(-1.0854032) q[2];
rz(-2.7392503) q[3];
sx q[3];
rz(-0.33807031) q[3];
sx q[3];
rz(-2.8344179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
