OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91035834) q[0];
sx q[0];
rz(-2.2725821) q[0];
sx q[0];
rz(-1.0847217) q[0];
rz(1.9864858) q[1];
sx q[1];
rz(-2.3218563) q[1];
sx q[1];
rz(0.91135946) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.203043) q[0];
sx q[0];
rz(-1.8126376) q[0];
sx q[0];
rz(-0.77745243) q[0];
x q[1];
rz(-1.4248104) q[2];
sx q[2];
rz(-2.2772191) q[2];
sx q[2];
rz(-1.3024769) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.54420097) q[1];
sx q[1];
rz(-2.4185838) q[1];
sx q[1];
rz(-0.71031481) q[1];
x q[2];
rz(0.064685589) q[3];
sx q[3];
rz(-1.1380592) q[3];
sx q[3];
rz(-2.329934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.800941) q[2];
sx q[2];
rz(-1.1110577) q[2];
sx q[2];
rz(0.079455376) q[2];
rz(-2.5331412) q[3];
sx q[3];
rz(-1.918101) q[3];
sx q[3];
rz(2.2505545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65421739) q[0];
sx q[0];
rz(-0.86367622) q[0];
sx q[0];
rz(-1.3551711) q[0];
rz(1.0379418) q[1];
sx q[1];
rz(-2.4619921) q[1];
sx q[1];
rz(-2.3341446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45189598) q[0];
sx q[0];
rz(-0.97366714) q[0];
sx q[0];
rz(-2.7886765) q[0];
rz(-0.37952559) q[2];
sx q[2];
rz(-2.1132601) q[2];
sx q[2];
rz(1.2778953) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8958853) q[1];
sx q[1];
rz(-1.6370607) q[1];
sx q[1];
rz(-2.9016728) q[1];
x q[2];
rz(0.1065527) q[3];
sx q[3];
rz(-1.7037183) q[3];
sx q[3];
rz(2.7850341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9281533) q[2];
sx q[2];
rz(-1.5636874) q[2];
sx q[2];
rz(-1.6820924) q[2];
rz(-0.45808074) q[3];
sx q[3];
rz(-0.57772485) q[3];
sx q[3];
rz(-0.95988449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9816575) q[0];
sx q[0];
rz(-2.0913048) q[0];
sx q[0];
rz(-0.67498573) q[0];
rz(-2.5534897) q[1];
sx q[1];
rz(-0.85398483) q[1];
sx q[1];
rz(1.3311707) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0597525) q[0];
sx q[0];
rz(-2.4911103) q[0];
sx q[0];
rz(-2.4516134) q[0];
x q[1];
rz(2.476112) q[2];
sx q[2];
rz(-1.988171) q[2];
sx q[2];
rz(2.5145384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.71215502) q[1];
sx q[1];
rz(-0.78199103) q[1];
sx q[1];
rz(-1.9351134) q[1];
rz(-1.0717594) q[3];
sx q[3];
rz(-1.3924358) q[3];
sx q[3];
rz(-0.87770977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.62446928) q[2];
sx q[2];
rz(-0.34999592) q[2];
sx q[2];
rz(2.9208753) q[2];
rz(-2.3366426) q[3];
sx q[3];
rz(-1.5419818) q[3];
sx q[3];
rz(1.8055003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6014366) q[0];
sx q[0];
rz(-0.24208459) q[0];
sx q[0];
rz(2.5228187) q[0];
rz(-0.082911804) q[1];
sx q[1];
rz(-0.34268788) q[1];
sx q[1];
rz(-1.6212911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0205958) q[0];
sx q[0];
rz(-1.474664) q[0];
sx q[0];
rz(2.0970048) q[0];
rz(2.7882468) q[2];
sx q[2];
rz(-2.5328703) q[2];
sx q[2];
rz(-0.48369394) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6449086) q[1];
sx q[1];
rz(-2.0808947) q[1];
sx q[1];
rz(-3.0472075) q[1];
rz(-pi) q[2];
rz(1.7509364) q[3];
sx q[3];
rz(-1.2015752) q[3];
sx q[3];
rz(2.957291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3380022) q[2];
sx q[2];
rz(-0.10927304) q[2];
sx q[2];
rz(-2.9962311) q[2];
rz(0.27211443) q[3];
sx q[3];
rz(-1.4067255) q[3];
sx q[3];
rz(-1.1468148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014932545) q[0];
sx q[0];
rz(-1.3260051) q[0];
sx q[0];
rz(3.0354011) q[0];
rz(-0.95942489) q[1];
sx q[1];
rz(-0.54490772) q[1];
sx q[1];
rz(0.19439654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9656669) q[0];
sx q[0];
rz(-2.3011294) q[0];
sx q[0];
rz(1.4552119) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.528605) q[2];
sx q[2];
rz(-0.11494177) q[2];
sx q[2];
rz(-1.7670222) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55479706) q[1];
sx q[1];
rz(-1.5298843) q[1];
sx q[1];
rz(2.552383) q[1];
rz(-1.7484457) q[3];
sx q[3];
rz(-2.7294949) q[3];
sx q[3];
rz(0.92890152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.83547366) q[2];
sx q[2];
rz(-1.6872311) q[2];
sx q[2];
rz(2.5731738) q[2];
rz(0.052637188) q[3];
sx q[3];
rz(-1.6391552) q[3];
sx q[3];
rz(-0.13681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0282054) q[0];
sx q[0];
rz(-2.5882692) q[0];
sx q[0];
rz(-2.9521039) q[0];
rz(-2.1151309) q[1];
sx q[1];
rz(-0.84954134) q[1];
sx q[1];
rz(-0.75526563) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8403977) q[0];
sx q[0];
rz(-0.9009255) q[0];
sx q[0];
rz(2.2706881) q[0];
rz(0.83024518) q[2];
sx q[2];
rz(-2.4355222) q[2];
sx q[2];
rz(2.490807) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4800075) q[1];
sx q[1];
rz(-2.9094041) q[1];
sx q[1];
rz(-0.8187553) q[1];
x q[2];
rz(-0.80249287) q[3];
sx q[3];
rz(-2.6518648) q[3];
sx q[3];
rz(-2.6306613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75817627) q[2];
sx q[2];
rz(-1.5864317) q[2];
sx q[2];
rz(1.4413393) q[2];
rz(-1.3759184) q[3];
sx q[3];
rz(-2.223189) q[3];
sx q[3];
rz(-0.70022303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.7021084) q[0];
sx q[0];
rz(-1.0478042) q[0];
sx q[0];
rz(-1.1442319) q[0];
rz(2.8817835) q[1];
sx q[1];
rz(-2.3016498) q[1];
sx q[1];
rz(-0.98794404) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081962498) q[0];
sx q[0];
rz(-0.10082997) q[0];
sx q[0];
rz(2.9743845) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6778474) q[2];
sx q[2];
rz(-0.79163359) q[2];
sx q[2];
rz(2.3196057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.64558951) q[1];
sx q[1];
rz(-2.2119446) q[1];
sx q[1];
rz(1.8899797) q[1];
rz(-0.99616237) q[3];
sx q[3];
rz(-1.8271433) q[3];
sx q[3];
rz(-2.4202895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40954956) q[2];
sx q[2];
rz(-1.6403551) q[2];
sx q[2];
rz(-0.6130971) q[2];
rz(0.71550718) q[3];
sx q[3];
rz(-0.77176538) q[3];
sx q[3];
rz(2.7806921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.20329423) q[0];
sx q[0];
rz(-0.84380117) q[0];
sx q[0];
rz(2.8734558) q[0];
rz(-2.9109491) q[1];
sx q[1];
rz(-1.5430887) q[1];
sx q[1];
rz(3.1275829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5174311) q[0];
sx q[0];
rz(-0.78963477) q[0];
sx q[0];
rz(-0.29615088) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1300541) q[2];
sx q[2];
rz(-1.4258372) q[2];
sx q[2];
rz(3.0359389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6790633) q[1];
sx q[1];
rz(-2.6025297) q[1];
sx q[1];
rz(-0.087058914) q[1];
rz(-pi) q[2];
rz(-3.1118891) q[3];
sx q[3];
rz(-1.4280768) q[3];
sx q[3];
rz(2.5102455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9739428) q[2];
sx q[2];
rz(-1.3385945) q[2];
sx q[2];
rz(-0.30926427) q[2];
rz(-0.15657982) q[3];
sx q[3];
rz(-1.4045818) q[3];
sx q[3];
rz(1.0369161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88368791) q[0];
sx q[0];
rz(-2.4990999) q[0];
sx q[0];
rz(0.25074348) q[0];
rz(2.5550487) q[1];
sx q[1];
rz(-1.5637014) q[1];
sx q[1];
rz(-2.3768545) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.667812) q[0];
sx q[0];
rz(-0.43536738) q[0];
sx q[0];
rz(0.82447585) q[0];
rz(2.370442) q[2];
sx q[2];
rz(-1.7131107) q[2];
sx q[2];
rz(-2.3754295) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2231317) q[1];
sx q[1];
rz(-1.4479965) q[1];
sx q[1];
rz(1.5614913) q[1];
rz(1.5161726) q[3];
sx q[3];
rz(-1.9491073) q[3];
sx q[3];
rz(2.9448933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5137198) q[2];
sx q[2];
rz(-0.79078117) q[2];
sx q[2];
rz(-2.9676843) q[2];
rz(0.74226132) q[3];
sx q[3];
rz(-2.3895013) q[3];
sx q[3];
rz(-0.043005634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1327508) q[0];
sx q[0];
rz(-2.371377) q[0];
sx q[0];
rz(2.8736864) q[0];
rz(-1.8065037) q[1];
sx q[1];
rz(-2.2309512) q[1];
sx q[1];
rz(-1.3722027) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55937476) q[0];
sx q[0];
rz(-1.1211044) q[0];
sx q[0];
rz(0.26828464) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4698974) q[2];
sx q[2];
rz(-1.865662) q[2];
sx q[2];
rz(3.1050668) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47169301) q[1];
sx q[1];
rz(-2.0834384) q[1];
sx q[1];
rz(1.0197082) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2413512) q[3];
sx q[3];
rz(-1.2678896) q[3];
sx q[3];
rz(-1.1166513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.043896) q[2];
sx q[2];
rz(-1.0363657) q[2];
sx q[2];
rz(-1.3807266) q[2];
rz(-2.0128287) q[3];
sx q[3];
rz(-0.94996101) q[3];
sx q[3];
rz(1.5560163) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783405) q[0];
sx q[0];
rz(-1.6077519) q[0];
sx q[0];
rz(2.0335249) q[0];
rz(-0.68645984) q[1];
sx q[1];
rz(-1.3931128) q[1];
sx q[1];
rz(-1.211094) q[1];
rz(0.49439597) q[2];
sx q[2];
rz(-1.313268) q[2];
sx q[2];
rz(3.1368844) q[2];
rz(1.0330647) q[3];
sx q[3];
rz(-1.929639) q[3];
sx q[3];
rz(0.7076984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
