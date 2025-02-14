OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.51337564) q[0];
sx q[0];
rz(-1.6802508) q[0];
sx q[0];
rz(2.2267377) q[0];
rz(2.272361) q[1];
sx q[1];
rz(-0.72328049) q[1];
sx q[1];
rz(-1.7887315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2268838) q[0];
sx q[0];
rz(-2.0830724) q[0];
sx q[0];
rz(-0.73988503) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74341821) q[2];
sx q[2];
rz(-1.8525436) q[2];
sx q[2];
rz(-2.8670058) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1754419) q[1];
sx q[1];
rz(-1.3805318) q[1];
sx q[1];
rz(-0.89118608) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9269591) q[3];
sx q[3];
rz(-2.1691536) q[3];
sx q[3];
rz(-0.90902381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4065518) q[2];
sx q[2];
rz(-0.81321365) q[2];
sx q[2];
rz(1.0068033) q[2];
rz(-1.462228) q[3];
sx q[3];
rz(-1.6961325) q[3];
sx q[3];
rz(1.538895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2891069) q[0];
sx q[0];
rz(-2.2302581) q[0];
sx q[0];
rz(3.015633) q[0];
rz(1.1360629) q[1];
sx q[1];
rz(-1.845153) q[1];
sx q[1];
rz(-2.0819285) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7505109) q[0];
sx q[0];
rz(-1.5727057) q[0];
sx q[0];
rz(-1.5715412) q[0];
rz(-pi) q[1];
rz(2.4074209) q[2];
sx q[2];
rz(-2.3245272) q[2];
sx q[2];
rz(1.2881607) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4821533) q[1];
sx q[1];
rz(-1.1350613) q[1];
sx q[1];
rz(-2.3547291) q[1];
rz(-0.35925389) q[3];
sx q[3];
rz(-1.6174966) q[3];
sx q[3];
rz(2.8253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.58553592) q[2];
sx q[2];
rz(-1.5289331) q[2];
sx q[2];
rz(-2.5118828) q[2];
rz(1.2398047) q[3];
sx q[3];
rz(-2.3972378) q[3];
sx q[3];
rz(2.0201717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12104163) q[0];
sx q[0];
rz(-0.23153767) q[0];
sx q[0];
rz(-1.9112021) q[0];
rz(-2.4379099) q[1];
sx q[1];
rz(-0.92862248) q[1];
sx q[1];
rz(1.5812662) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9890922) q[0];
sx q[0];
rz(-1.9643503) q[0];
sx q[0];
rz(3.1253184) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4034924) q[2];
sx q[2];
rz(-1.7154652) q[2];
sx q[2];
rz(-1.6044782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.81654233) q[1];
sx q[1];
rz(-2.5245164) q[1];
sx q[1];
rz(-0.80775921) q[1];
rz(-2.7625691) q[3];
sx q[3];
rz(-2.3190917) q[3];
sx q[3];
rz(1.8387972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33502093) q[2];
sx q[2];
rz(-2.0609914) q[2];
sx q[2];
rz(3.0873155) q[2];
rz(-2.055376) q[3];
sx q[3];
rz(-0.79038668) q[3];
sx q[3];
rz(-0.49200341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.7739968) q[0];
sx q[0];
rz(-0.01570276) q[0];
sx q[0];
rz(-0.98454654) q[0];
rz(1.7355851) q[1];
sx q[1];
rz(-1.5407298) q[1];
sx q[1];
rz(2.9071009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2024883) q[0];
sx q[0];
rz(-1.2959769) q[0];
sx q[0];
rz(2.2780134) q[0];
rz(-pi) q[1];
rz(2.8963887) q[2];
sx q[2];
rz(-1.8498932) q[2];
sx q[2];
rz(-2.9989105) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.6813184) q[1];
sx q[1];
rz(-1.360962) q[1];
sx q[1];
rz(2.5833681) q[1];
rz(-pi) q[2];
rz(-2.1470039) q[3];
sx q[3];
rz(-1.5458917) q[3];
sx q[3];
rz(-1.2664317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8977114) q[2];
sx q[2];
rz(-2.477024) q[2];
sx q[2];
rz(2.393874) q[2];
rz(-1.0866577) q[3];
sx q[3];
rz(-2.1472411) q[3];
sx q[3];
rz(2.7896816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20027593) q[0];
sx q[0];
rz(-0.87012297) q[0];
sx q[0];
rz(1.548832) q[0];
rz(3.0039655) q[1];
sx q[1];
rz(-0.96447861) q[1];
sx q[1];
rz(2.4639938) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9842868) q[0];
sx q[0];
rz(-1.5781227) q[0];
sx q[0];
rz(3.1276032) q[0];
x q[1];
rz(-1.4196012) q[2];
sx q[2];
rz(-1.7624904) q[2];
sx q[2];
rz(-0.58133948) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4487344) q[1];
sx q[1];
rz(-3.0046607) q[1];
sx q[1];
rz(-1.6137692) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8749488) q[3];
sx q[3];
rz(-2.3867749) q[3];
sx q[3];
rz(-2.9637869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4981093) q[2];
sx q[2];
rz(-2.6124239) q[2];
sx q[2];
rz(-2.7672178) q[2];
rz(-2.3937285) q[3];
sx q[3];
rz(-1.6467843) q[3];
sx q[3];
rz(-1.098746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13136524) q[0];
sx q[0];
rz(-1.7759198) q[0];
sx q[0];
rz(-0.22450547) q[0];
rz(0.35047105) q[1];
sx q[1];
rz(-1.1843362) q[1];
sx q[1];
rz(-1.9559466) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4476337) q[0];
sx q[0];
rz(-2.9835851) q[0];
sx q[0];
rz(-1.9290646) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4746488) q[2];
sx q[2];
rz(-1.1091057) q[2];
sx q[2];
rz(-2.8682414) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95690475) q[1];
sx q[1];
rz(-0.84057759) q[1];
sx q[1];
rz(-0.59112494) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91573651) q[3];
sx q[3];
rz(-1.011542) q[3];
sx q[3];
rz(2.9489234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6422358) q[2];
sx q[2];
rz(-0.88358742) q[2];
sx q[2];
rz(0.16481608) q[2];
rz(-0.76588255) q[3];
sx q[3];
rz(-0.82847786) q[3];
sx q[3];
rz(-0.59144056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.0931452) q[0];
sx q[0];
rz(-0.026938139) q[0];
sx q[0];
rz(2.7222166) q[0];
rz(-1.558149) q[1];
sx q[1];
rz(-2.7132468) q[1];
sx q[1];
rz(1.8289808) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93363491) q[0];
sx q[0];
rz(-2.5267753) q[0];
sx q[0];
rz(-1.8551679) q[0];
rz(-pi) q[1];
rz(-0.68676853) q[2];
sx q[2];
rz(-2.1109952) q[2];
sx q[2];
rz(1.5093925) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53322843) q[1];
sx q[1];
rz(-1.3238412) q[1];
sx q[1];
rz(0.19503959) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26793865) q[3];
sx q[3];
rz(-1.4456141) q[3];
sx q[3];
rz(2.1942558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7340362) q[2];
sx q[2];
rz(-2.4702256) q[2];
sx q[2];
rz(-2.7835795) q[2];
rz(2.9426306) q[3];
sx q[3];
rz(-2.3301061) q[3];
sx q[3];
rz(-0.097804047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99893779) q[0];
sx q[0];
rz(-1.0336579) q[0];
sx q[0];
rz(2.348483) q[0];
rz(-2.2456887) q[1];
sx q[1];
rz(-1.9205807) q[1];
sx q[1];
rz(-2.6633247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60837854) q[0];
sx q[0];
rz(-1.694931) q[0];
sx q[0];
rz(-0.45093765) q[0];
x q[1];
rz(1.2094648) q[2];
sx q[2];
rz(-1.1602957) q[2];
sx q[2];
rz(2.362006) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1898124) q[1];
sx q[1];
rz(-2.5972022) q[1];
sx q[1];
rz(-1.0448827) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0555154) q[3];
sx q[3];
rz(-1.4751504) q[3];
sx q[3];
rz(-0.0064004504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6200977) q[2];
sx q[2];
rz(-1.5767117) q[2];
sx q[2];
rz(3.133797) q[2];
rz(1.1826285) q[3];
sx q[3];
rz(-1.7595485) q[3];
sx q[3];
rz(-1.8462605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9919392) q[0];
sx q[0];
rz(-0.45435926) q[0];
sx q[0];
rz(-0.40865189) q[0];
rz(-2.0463792) q[1];
sx q[1];
rz(-1.6588255) q[1];
sx q[1];
rz(-2.2832787) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33108157) q[0];
sx q[0];
rz(-0.16180049) q[0];
sx q[0];
rz(-1.1463548) q[0];
x q[1];
rz(-1.3038396) q[2];
sx q[2];
rz(-0.46707312) q[2];
sx q[2];
rz(2.2427223) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6950001) q[1];
sx q[1];
rz(-1.6808676) q[1];
sx q[1];
rz(1.7677444) q[1];
rz(-2.8204042) q[3];
sx q[3];
rz(-0.62231718) q[3];
sx q[3];
rz(0.27143196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7909164) q[2];
sx q[2];
rz(-2.3886267) q[2];
sx q[2];
rz(-0.6616627) q[2];
rz(2.3412797) q[3];
sx q[3];
rz(-1.9789663) q[3];
sx q[3];
rz(0.77320981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8372339) q[0];
sx q[0];
rz(-0.14166129) q[0];
sx q[0];
rz(2.4270571) q[0];
rz(2.6658658) q[1];
sx q[1];
rz(-2.5560502) q[1];
sx q[1];
rz(0.48019662) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6960951) q[0];
sx q[0];
rz(-1.0871743) q[0];
sx q[0];
rz(0.66894834) q[0];
x q[1];
rz(-0.91861208) q[2];
sx q[2];
rz(-1.3038074) q[2];
sx q[2];
rz(0.3912386) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.46327766) q[1];
sx q[1];
rz(-2.2147605) q[1];
sx q[1];
rz(-1.2106845) q[1];
rz(2.1289292) q[3];
sx q[3];
rz(-0.50202657) q[3];
sx q[3];
rz(3.0691911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51018888) q[2];
sx q[2];
rz(-1.9437342) q[2];
sx q[2];
rz(-0.27604827) q[2];
rz(2.7035642) q[3];
sx q[3];
rz(-1.5166401) q[3];
sx q[3];
rz(2.5146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.940687) q[0];
sx q[0];
rz(-1.9872266) q[0];
sx q[0];
rz(-3.0845597) q[0];
rz(-0.036046473) q[1];
sx q[1];
rz(-1.6135975) q[1];
sx q[1];
rz(1.5560908) q[1];
rz(-2.8677058) q[2];
sx q[2];
rz(-2.0487493) q[2];
sx q[2];
rz(-0.4977939) q[2];
rz(-1.4344169) q[3];
sx q[3];
rz(-1.3188667) q[3];
sx q[3];
rz(1.8018166) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
