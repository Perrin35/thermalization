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
rz(-0.43871969) q[0];
sx q[0];
rz(-2.5320142) q[0];
sx q[0];
rz(0.69019812) q[0];
rz(-1.8742427) q[1];
sx q[1];
rz(-2.6522418) q[1];
sx q[1];
rz(0.34440053) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1699088) q[0];
sx q[0];
rz(-0.2580041) q[0];
sx q[0];
rz(0.95746104) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6292388) q[2];
sx q[2];
rz(-0.47253451) q[2];
sx q[2];
rz(-2.3290881) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92235451) q[1];
sx q[1];
rz(-1.7686756) q[1];
sx q[1];
rz(1.3082275) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1719443) q[3];
sx q[3];
rz(-2.056582) q[3];
sx q[3];
rz(-2.7051089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0979746) q[2];
sx q[2];
rz(-2.2100885) q[2];
sx q[2];
rz(-1.0574868) q[2];
rz(2.1438694) q[3];
sx q[3];
rz(-1.7505373) q[3];
sx q[3];
rz(1.1725496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69720307) q[0];
sx q[0];
rz(-2.3219705) q[0];
sx q[0];
rz(0.75253734) q[0];
rz(-1.2731816) q[1];
sx q[1];
rz(-1.1538785) q[1];
sx q[1];
rz(-2.2862327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68356791) q[0];
sx q[0];
rz(-0.23580256) q[0];
sx q[0];
rz(2.0033512) q[0];
x q[1];
rz(0.24829243) q[2];
sx q[2];
rz(-0.96370164) q[2];
sx q[2];
rz(2.8090629) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2084045) q[1];
sx q[1];
rz(-1.4074086) q[1];
sx q[1];
rz(-2.7524314) q[1];
x q[2];
rz(0.19892502) q[3];
sx q[3];
rz(-1.5459877) q[3];
sx q[3];
rz(-2.5682784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37318834) q[2];
sx q[2];
rz(-1.5311925) q[2];
sx q[2];
rz(1.1676403) q[2];
rz(-3.043637) q[3];
sx q[3];
rz(-0.28147134) q[3];
sx q[3];
rz(1.9907985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34404594) q[0];
sx q[0];
rz(-0.56586376) q[0];
sx q[0];
rz(1.6746445) q[0];
rz(3.103718) q[1];
sx q[1];
rz(-1.2118309) q[1];
sx q[1];
rz(-2.0272592) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5112202) q[0];
sx q[0];
rz(-1.9812222) q[0];
sx q[0];
rz(0.12271304) q[0];
rz(-pi) q[1];
rz(0.55338316) q[2];
sx q[2];
rz(-1.1107365) q[2];
sx q[2];
rz(-1.9529238) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13521299) q[1];
sx q[1];
rz(-2.4828382) q[1];
sx q[1];
rz(-0.82243311) q[1];
x q[2];
rz(0.98975269) q[3];
sx q[3];
rz(-1.8268181) q[3];
sx q[3];
rz(0.44629249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.36180878) q[2];
sx q[2];
rz(-0.59541687) q[2];
sx q[2];
rz(-0.53100604) q[2];
rz(0.94046721) q[3];
sx q[3];
rz(-2.2898424) q[3];
sx q[3];
rz(1.4550335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4130212) q[0];
sx q[0];
rz(-0.41550264) q[0];
sx q[0];
rz(-0.79237932) q[0];
rz(0.13262311) q[1];
sx q[1];
rz(-0.68999973) q[1];
sx q[1];
rz(-2.2161868) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9856095) q[0];
sx q[0];
rz(-1.7260084) q[0];
sx q[0];
rz(-1.3830091) q[0];
x q[1];
rz(1.4508574) q[2];
sx q[2];
rz(-2.687157) q[2];
sx q[2];
rz(1.1392405) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1455411) q[1];
sx q[1];
rz(-0.62623238) q[1];
sx q[1];
rz(1.9171417) q[1];
x q[2];
rz(-0.53823353) q[3];
sx q[3];
rz(-0.58659121) q[3];
sx q[3];
rz(0.80249062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35599071) q[2];
sx q[2];
rz(-1.2987368) q[2];
sx q[2];
rz(1.7388657) q[2];
rz(1.8897024) q[3];
sx q[3];
rz(-1.8391049) q[3];
sx q[3];
rz(0.29786626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-1.4715618) q[1];
sx q[1];
rz(0.54738799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.669466) q[0];
sx q[0];
rz(-1.5385813) q[0];
sx q[0];
rz(2.6526582) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6519189) q[2];
sx q[2];
rz(-1.2077959) q[2];
sx q[2];
rz(-2.2081019) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84505397) q[1];
sx q[1];
rz(-0.88857874) q[1];
sx q[1];
rz(-0.53081546) q[1];
rz(-pi) q[2];
rz(0.27606729) q[3];
sx q[3];
rz(-1.4763586) q[3];
sx q[3];
rz(-2.6219079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0543412) q[2];
sx q[2];
rz(-2.6332899) q[2];
sx q[2];
rz(0.15092078) q[2];
rz(-1.5058676) q[3];
sx q[3];
rz(-1.8412291) q[3];
sx q[3];
rz(1.6747564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48581377) q[0];
sx q[0];
rz(-1.9447615) q[0];
sx q[0];
rz(-0.47564137) q[0];
rz(0.42426839) q[1];
sx q[1];
rz(-1.905966) q[1];
sx q[1];
rz(-1.7686527) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.951648) q[0];
sx q[0];
rz(-1.7423672) q[0];
sx q[0];
rz(-2.8749332) q[0];
x q[1];
rz(-2.0654997) q[2];
sx q[2];
rz(-2.1044136) q[2];
sx q[2];
rz(0.99737203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13481278) q[1];
sx q[1];
rz(-1.0357417) q[1];
sx q[1];
rz(-0.3877918) q[1];
rz(0.11656363) q[3];
sx q[3];
rz(-1.1313843) q[3];
sx q[3];
rz(-1.9842465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0900535) q[2];
sx q[2];
rz(-2.1681483) q[2];
sx q[2];
rz(-3.0494087) q[2];
rz(-1.9696382) q[3];
sx q[3];
rz(-2.6086174) q[3];
sx q[3];
rz(0.91641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5343269) q[0];
sx q[0];
rz(-2.3054275) q[0];
sx q[0];
rz(-1.0323866) q[0];
rz(-3.0119925) q[1];
sx q[1];
rz(-1.7584636) q[1];
sx q[1];
rz(-1.3547156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537405) q[0];
sx q[0];
rz(-0.82621413) q[0];
sx q[0];
rz(0.67954845) q[0];
rz(-pi) q[1];
rz(1.9117891) q[2];
sx q[2];
rz(-1.908506) q[2];
sx q[2];
rz(-2.4203398) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4413905) q[1];
sx q[1];
rz(-1.017015) q[1];
sx q[1];
rz(0.9662083) q[1];
rz(-1.4192102) q[3];
sx q[3];
rz(-0.48036423) q[3];
sx q[3];
rz(2.9472873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9356392) q[2];
sx q[2];
rz(-2.8796068) q[2];
sx q[2];
rz(1.4865173) q[2];
rz(-2.4691811) q[3];
sx q[3];
rz(-1.0786062) q[3];
sx q[3];
rz(-3.0730754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2023778) q[0];
sx q[0];
rz(-3.0164533) q[0];
sx q[0];
rz(0.2739206) q[0];
rz(0.1637474) q[1];
sx q[1];
rz(-1.3122908) q[1];
sx q[1];
rz(-0.40072498) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3253097) q[0];
sx q[0];
rz(-1.0772711) q[0];
sx q[0];
rz(-1.8614344) q[0];
rz(-3.0178304) q[2];
sx q[2];
rz(-1.971619) q[2];
sx q[2];
rz(-1.4126029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5536062) q[1];
sx q[1];
rz(-0.67236908) q[1];
sx q[1];
rz(1.5410627) q[1];
rz(-pi) q[2];
rz(-0.89426453) q[3];
sx q[3];
rz(-2.6163963) q[3];
sx q[3];
rz(-2.0390455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87216941) q[2];
sx q[2];
rz(-2.9672406) q[2];
sx q[2];
rz(-1.6677469) q[2];
rz(-3.040124) q[3];
sx q[3];
rz(-1.2227819) q[3];
sx q[3];
rz(1.7469223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87429109) q[0];
sx q[0];
rz(-0.75031459) q[0];
sx q[0];
rz(3.1399723) q[0];
rz(2.5770309) q[1];
sx q[1];
rz(-1.8058585) q[1];
sx q[1];
rz(-2.6394305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9410132) q[0];
sx q[0];
rz(-0.94958011) q[0];
sx q[0];
rz(-2.4401779) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61567581) q[2];
sx q[2];
rz(-2.4936495) q[2];
sx q[2];
rz(0.16251646) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.069418176) q[1];
sx q[1];
rz(-2.0771306) q[1];
sx q[1];
rz(-1.7933374) q[1];
rz(-pi) q[2];
rz(-1.05324) q[3];
sx q[3];
rz(-1.386214) q[3];
sx q[3];
rz(-1.7257168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0994215) q[2];
sx q[2];
rz(-1.9440938) q[2];
sx q[2];
rz(1.8966804) q[2];
rz(2.9391607) q[3];
sx q[3];
rz(-1.3916241) q[3];
sx q[3];
rz(-2.1421053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779483) q[0];
sx q[0];
rz(-0.36643323) q[0];
sx q[0];
rz(-1.7250489) q[0];
rz(-1.1997403) q[1];
sx q[1];
rz(-1.1734633) q[1];
sx q[1];
rz(-2.8660668) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4107127) q[0];
sx q[0];
rz(-1.2322324) q[0];
sx q[0];
rz(1.186446) q[0];
rz(-pi) q[1];
rz(-2.7904835) q[2];
sx q[2];
rz(-1.654823) q[2];
sx q[2];
rz(2.71738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80088204) q[1];
sx q[1];
rz(-1.5499797) q[1];
sx q[1];
rz(2.0321991) q[1];
rz(0.59861981) q[3];
sx q[3];
rz(-1.4685923) q[3];
sx q[3];
rz(0.24598611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3755017) q[2];
sx q[2];
rz(-1.4644863) q[2];
sx q[2];
rz(-0.38140934) q[2];
rz(2.8219847) q[3];
sx q[3];
rz(-2.3242798) q[3];
sx q[3];
rz(-0.6680502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1565336) q[0];
sx q[0];
rz(-1.986713) q[0];
sx q[0];
rz(2.4646133) q[0];
rz(2.1398687) q[1];
sx q[1];
rz(-1.6698508) q[1];
sx q[1];
rz(2.225266) q[1];
rz(-1.4475293) q[2];
sx q[2];
rz(-0.90926778) q[2];
sx q[2];
rz(-2.8127083) q[2];
rz(-1.6826717) q[3];
sx q[3];
rz(-2.2266012) q[3];
sx q[3];
rz(0.83960017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
