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
rz(2.702873) q[0];
sx q[0];
rz(-0.60957849) q[0];
sx q[0];
rz(2.4513945) q[0];
rz(1.26735) q[1];
sx q[1];
rz(-0.48935088) q[1];
sx q[1];
rz(2.7971921) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9448954) q[0];
sx q[0];
rz(-1.718194) q[0];
sx q[0];
rz(-1.3582673) q[0];
x q[1];
rz(1.3252668) q[2];
sx q[2];
rz(-1.9787162) q[2];
sx q[2];
rz(-1.7657042) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8616383) q[1];
sx q[1];
rz(-0.32740232) q[1];
sx q[1];
rz(-2.2285217) q[1];
rz(-pi) q[2];
rz(-2.0627414) q[3];
sx q[3];
rz(-1.4189093) q[3];
sx q[3];
rz(-1.2152176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0979746) q[2];
sx q[2];
rz(-2.2100885) q[2];
sx q[2];
rz(1.0574868) q[2];
rz(-0.99772325) q[3];
sx q[3];
rz(-1.7505373) q[3];
sx q[3];
rz(-1.9690431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69720307) q[0];
sx q[0];
rz(-2.3219705) q[0];
sx q[0];
rz(-0.75253734) q[0];
rz(-1.2731816) q[1];
sx q[1];
rz(-1.9877142) q[1];
sx q[1];
rz(2.2862327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0146926) q[0];
sx q[0];
rz(-1.784526) q[0];
sx q[0];
rz(-3.041211) q[0];
rz(-pi) q[1];
rz(1.2307421) q[2];
sx q[2];
rz(-2.4916561) q[2];
sx q[2];
rz(0.75071819) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2084045) q[1];
sx q[1];
rz(-1.734184) q[1];
sx q[1];
rz(-2.7524314) q[1];
rz(-pi) q[2];
rz(0.12491183) q[3];
sx q[3];
rz(-0.20044573) q[3];
sx q[3];
rz(-0.87504609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7684043) q[2];
sx q[2];
rz(-1.5311925) q[2];
sx q[2];
rz(1.1676403) q[2];
rz(-3.043637) q[3];
sx q[3];
rz(-2.8601213) q[3];
sx q[3];
rz(-1.9907985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34404594) q[0];
sx q[0];
rz(-2.5757289) q[0];
sx q[0];
rz(-1.4669482) q[0];
rz(-3.103718) q[1];
sx q[1];
rz(-1.9297618) q[1];
sx q[1];
rz(1.1143335) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1087462) q[0];
sx q[0];
rz(-1.683273) q[0];
sx q[0];
rz(1.1576018) q[0];
x q[1];
rz(0.75598209) q[2];
sx q[2];
rz(-2.4376483) q[2];
sx q[2];
rz(2.9005425) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0689286) q[1];
sx q[1];
rz(-1.1411828) q[1];
sx q[1];
rz(1.0547897) q[1];
rz(-2.8381078) q[3];
sx q[3];
rz(-1.0110098) q[3];
sx q[3];
rz(1.289285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36180878) q[2];
sx q[2];
rz(-2.5461758) q[2];
sx q[2];
rz(-0.53100604) q[2];
rz(2.2011254) q[3];
sx q[3];
rz(-2.2898424) q[3];
sx q[3];
rz(-1.4550335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4130212) q[0];
sx q[0];
rz(-2.72609) q[0];
sx q[0];
rz(-0.79237932) q[0];
rz(0.13262311) q[1];
sx q[1];
rz(-0.68999973) q[1];
sx q[1];
rz(0.92540583) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6974119) q[0];
sx q[0];
rz(-1.385293) q[0];
sx q[0];
rz(0.15794396) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6907352) q[2];
sx q[2];
rz(-0.45443568) q[2];
sx q[2];
rz(1.1392405) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5645743) q[1];
sx q[1];
rz(-0.98688236) q[1];
sx q[1];
rz(2.9008032) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6033591) q[3];
sx q[3];
rz(-2.5550014) q[3];
sx q[3];
rz(-0.80249062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7856019) q[2];
sx q[2];
rz(-1.8428558) q[2];
sx q[2];
rz(1.7388657) q[2];
rz(-1.8897024) q[3];
sx q[3];
rz(-1.8391049) q[3];
sx q[3];
rz(-0.29786626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4414325) q[0];
sx q[0];
rz(-2.1692363) q[0];
sx q[0];
rz(-2.1527619) q[0];
rz(-1.5160457) q[1];
sx q[1];
rz(-1.6700309) q[1];
sx q[1];
rz(2.5942047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.103391) q[0];
sx q[0];
rz(-0.48990881) q[0];
sx q[0];
rz(-3.0730877) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21018942) q[2];
sx q[2];
rz(-2.7700305) q[2];
sx q[2];
rz(-0.7084058) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0802287) q[1];
sx q[1];
rz(-1.974611) q[1];
sx q[1];
rz(-2.3263127) q[1];
x q[2];
rz(1.6689273) q[3];
sx q[3];
rz(-1.8456015) q[3];
sx q[3];
rz(2.0637728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0543412) q[2];
sx q[2];
rz(-2.6332899) q[2];
sx q[2];
rz(-2.9906719) q[2];
rz(-1.5058676) q[3];
sx q[3];
rz(-1.3003636) q[3];
sx q[3];
rz(-1.6747564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.48581377) q[0];
sx q[0];
rz(-1.9447615) q[0];
sx q[0];
rz(2.6659513) q[0];
rz(-2.7173243) q[1];
sx q[1];
rz(-1.905966) q[1];
sx q[1];
rz(-1.7686527) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9638393) q[0];
sx q[0];
rz(-2.8256157) q[0];
sx q[0];
rz(0.58167268) q[0];
x q[1];
rz(0.59117909) q[2];
sx q[2];
rz(-1.9918963) q[2];
sx q[2];
rz(2.3003701) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81022963) q[1];
sx q[1];
rz(-0.64953564) q[1];
sx q[1];
rz(-1.002921) q[1];
rz(1.3282534) q[3];
sx q[3];
rz(-2.6879592) q[3];
sx q[3];
rz(-0.88874528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0515392) q[2];
sx q[2];
rz(-2.1681483) q[2];
sx q[2];
rz(3.0494087) q[2];
rz(-1.1719545) q[3];
sx q[3];
rz(-2.6086174) q[3];
sx q[3];
rz(2.2251825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6072657) q[0];
sx q[0];
rz(-0.83616513) q[0];
sx q[0];
rz(2.109206) q[0];
rz(-0.1296002) q[1];
sx q[1];
rz(-1.7584636) q[1];
sx q[1];
rz(-1.786877) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.860703) q[0];
sx q[0];
rz(-2.1797545) q[0];
sx q[0];
rz(0.97228284) q[0];
rz(-pi) q[1];
rz(-0.35667901) q[2];
sx q[2];
rz(-1.891815) q[2];
sx q[2];
rz(2.4090846) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9216283) q[1];
sx q[1];
rz(-0.79571834) q[1];
sx q[1];
rz(-0.74337767) q[1];
x q[2];
rz(1.7223825) q[3];
sx q[3];
rz(-2.6612284) q[3];
sx q[3];
rz(-2.9472873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9356392) q[2];
sx q[2];
rz(-0.26198584) q[2];
sx q[2];
rz(-1.4865173) q[2];
rz(2.4691811) q[3];
sx q[3];
rz(-1.0786062) q[3];
sx q[3];
rz(-0.068517223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9392149) q[0];
sx q[0];
rz(-3.0164533) q[0];
sx q[0];
rz(2.8676721) q[0];
rz(2.9778453) q[1];
sx q[1];
rz(-1.3122908) q[1];
sx q[1];
rz(-2.7408677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3794601) q[0];
sx q[0];
rz(-0.5665938) q[0];
sx q[0];
rz(-2.6521126) q[0];
rz(-0.12376229) q[2];
sx q[2];
rz(-1.1699737) q[2];
sx q[2];
rz(-1.4126029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54998794) q[1];
sx q[1];
rz(-2.2428136) q[1];
sx q[1];
rz(-3.1179291) q[1];
rz(0.89426453) q[3];
sx q[3];
rz(-0.52519631) q[3];
sx q[3];
rz(-2.0390455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.87216941) q[2];
sx q[2];
rz(-2.9672406) q[2];
sx q[2];
rz(-1.6677469) q[2];
rz(0.10146865) q[3];
sx q[3];
rz(-1.2227819) q[3];
sx q[3];
rz(1.7469223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2673016) q[0];
sx q[0];
rz(-0.75031459) q[0];
sx q[0];
rz(3.1399723) q[0];
rz(0.56456176) q[1];
sx q[1];
rz(-1.3357342) q[1];
sx q[1];
rz(-2.6394305) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.314437) q[0];
sx q[0];
rz(-1.0183882) q[0];
sx q[0];
rz(-2.3236426) q[0];
x q[1];
rz(-0.55353005) q[2];
sx q[2];
rz(-1.9268254) q[2];
sx q[2];
rz(-1.9218307) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50599033) q[1];
sx q[1];
rz(-0.5491623) q[1];
sx q[1];
rz(-0.378774) q[1];
rz(-pi) q[2];
rz(-1.2099482) q[3];
sx q[3];
rz(-0.54665414) q[3];
sx q[3];
rz(-2.9846559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.042171176) q[2];
sx q[2];
rz(-1.1974988) q[2];
sx q[2];
rz(-1.8966804) q[2];
rz(-0.20243195) q[3];
sx q[3];
rz(-1.3916241) q[3];
sx q[3];
rz(-2.1421053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56364432) q[0];
sx q[0];
rz(-2.7751594) q[0];
sx q[0];
rz(1.4165437) q[0];
rz(-1.1997403) q[1];
sx q[1];
rz(-1.9681294) q[1];
sx q[1];
rz(2.8660668) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4352132) q[0];
sx q[0];
rz(-1.2093028) q[0];
sx q[0];
rz(2.778591) q[0];
rz(-pi) q[1];
rz(1.4813384) q[2];
sx q[2];
rz(-1.2209792) q[2];
sx q[2];
rz(-1.1773156) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4135189) q[1];
sx q[1];
rz(-0.46183837) q[1];
sx q[1];
rz(1.6175265) q[1];
rz(-0.59861981) q[3];
sx q[3];
rz(-1.4685923) q[3];
sx q[3];
rz(-0.24598611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.766091) q[2];
sx q[2];
rz(-1.4644863) q[2];
sx q[2];
rz(-0.38140934) q[2];
rz(-2.8219847) q[3];
sx q[3];
rz(-0.8173129) q[3];
sx q[3];
rz(-0.6680502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850591) q[0];
sx q[0];
rz(-1.1548797) q[0];
sx q[0];
rz(-0.67697939) q[0];
rz(-2.1398687) q[1];
sx q[1];
rz(-1.4717419) q[1];
sx q[1];
rz(-0.91632661) q[1];
rz(-0.15663319) q[2];
sx q[2];
rz(-0.67121438) q[2];
sx q[2];
rz(0.12987655) q[2];
rz(0.14409625) q[3];
sx q[3];
rz(-2.477705) q[3];
sx q[3];
rz(-2.4841819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
