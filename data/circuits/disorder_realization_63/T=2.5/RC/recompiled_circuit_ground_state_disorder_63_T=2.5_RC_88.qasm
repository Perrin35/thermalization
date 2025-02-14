OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11290057) q[0];
sx q[0];
rz(2.8960462) q[0];
sx q[0];
rz(9.4774376) q[0];
rz(1.117299) q[1];
sx q[1];
rz(-2.8007562) q[1];
sx q[1];
rz(-2.053082) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65050012) q[0];
sx q[0];
rz(-0.40084565) q[0];
sx q[0];
rz(1.4087138) q[0];
x q[1];
rz(-0.91204466) q[2];
sx q[2];
rz(-0.17870644) q[2];
sx q[2];
rz(2.2026874) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6793078) q[1];
sx q[1];
rz(-2.4801697) q[1];
sx q[1];
rz(0.2703305) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5616712) q[3];
sx q[3];
rz(-0.41968981) q[3];
sx q[3];
rz(1.4008824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6344305) q[2];
sx q[2];
rz(-1.3619225) q[2];
sx q[2];
rz(0.6081028) q[2];
rz(-2.0173343) q[3];
sx q[3];
rz(-1.2048771) q[3];
sx q[3];
rz(2.2500136) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267777) q[0];
sx q[0];
rz(-1.5730653) q[0];
sx q[0];
rz(2.538105) q[0];
rz(-2.6314349) q[1];
sx q[1];
rz(-0.77168232) q[1];
sx q[1];
rz(-2.1147494) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6525243) q[0];
sx q[0];
rz(-1.7260255) q[0];
sx q[0];
rz(2.9672769) q[0];
x q[1];
rz(-2.2910094) q[2];
sx q[2];
rz(-0.89909485) q[2];
sx q[2];
rz(1.4631866) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3852147) q[1];
sx q[1];
rz(-2.2526388) q[1];
sx q[1];
rz(2.4663976) q[1];
rz(-2.6072826) q[3];
sx q[3];
rz(-2.1624193) q[3];
sx q[3];
rz(-0.1402157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.55871955) q[2];
sx q[2];
rz(-1.9931953) q[2];
sx q[2];
rz(0.64644512) q[2];
rz(-1.8630113) q[3];
sx q[3];
rz(-2.1918112) q[3];
sx q[3];
rz(-1.4518552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.932514) q[0];
sx q[0];
rz(-0.13020733) q[0];
sx q[0];
rz(1.5721488) q[0];
rz(-2.7945844) q[1];
sx q[1];
rz(-1.6667112) q[1];
sx q[1];
rz(-2.2775547) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0873252) q[0];
sx q[0];
rz(-1.010276) q[0];
sx q[0];
rz(1.818324) q[0];
x q[1];
rz(1.8225602) q[2];
sx q[2];
rz(-1.994311) q[2];
sx q[2];
rz(-0.26811312) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7622704) q[1];
sx q[1];
rz(-2.218609) q[1];
sx q[1];
rz(0.76384441) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7517463) q[3];
sx q[3];
rz(-2.6084508) q[3];
sx q[3];
rz(1.0764493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.049909264) q[2];
sx q[2];
rz(-0.20541643) q[2];
sx q[2];
rz(1.4570215) q[2];
rz(1.0157887) q[3];
sx q[3];
rz(-2.6088645) q[3];
sx q[3];
rz(-0.31639019) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006573) q[0];
sx q[0];
rz(-2.7641986) q[0];
sx q[0];
rz(1.5628016) q[0];
rz(-1.0904301) q[1];
sx q[1];
rz(-0.54160392) q[1];
sx q[1];
rz(-1.9900367) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3131696) q[0];
sx q[0];
rz(-2.1439664) q[0];
sx q[0];
rz(-1.8445119) q[0];
rz(-0.76049034) q[2];
sx q[2];
rz(-1.9273238) q[2];
sx q[2];
rz(-0.52887756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87428625) q[1];
sx q[1];
rz(-1.5360254) q[1];
sx q[1];
rz(1.0465996) q[1];
rz(1.5308446) q[3];
sx q[3];
rz(-0.60725799) q[3];
sx q[3];
rz(1.1049529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2942723) q[2];
sx q[2];
rz(-2.5966094) q[2];
sx q[2];
rz(-1.4997743) q[2];
rz(-3.0070987) q[3];
sx q[3];
rz(-1.8650863) q[3];
sx q[3];
rz(2.3597778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0624369) q[0];
sx q[0];
rz(-0.9117313) q[0];
sx q[0];
rz(-2.675918) q[0];
rz(-1.8648719) q[1];
sx q[1];
rz(-1.0036422) q[1];
sx q[1];
rz(1.0328971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55556923) q[0];
sx q[0];
rz(-2.9831671) q[0];
sx q[0];
rz(-1.6371631) q[0];
x q[1];
rz(-3.0014702) q[2];
sx q[2];
rz(-1.6701785) q[2];
sx q[2];
rz(1.6665719) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.031449854) q[1];
sx q[1];
rz(-1.3198766) q[1];
sx q[1];
rz(0.67278426) q[1];
rz(2.1139268) q[3];
sx q[3];
rz(-1.2333721) q[3];
sx q[3];
rz(1.8130482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59139645) q[2];
sx q[2];
rz(-1.1893716) q[2];
sx q[2];
rz(-2.2097394) q[2];
rz(3.0469117) q[3];
sx q[3];
rz(-2.2414312) q[3];
sx q[3];
rz(1.3100821) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4513627) q[0];
sx q[0];
rz(-2.9424423) q[0];
sx q[0];
rz(-3.0354101) q[0];
rz(-0.55446082) q[1];
sx q[1];
rz(-0.78840557) q[1];
sx q[1];
rz(-1.5415446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0672749) q[0];
sx q[0];
rz(-1.5467001) q[0];
sx q[0];
rz(2.7510608) q[0];
x q[1];
rz(-2.0279858) q[2];
sx q[2];
rz(-1.3234252) q[2];
sx q[2];
rz(0.24277011) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8733124) q[1];
sx q[1];
rz(-2.2619216) q[1];
sx q[1];
rz(-1.7296438) q[1];
rz(-pi) q[2];
rz(-2.8951449) q[3];
sx q[3];
rz(-1.6250027) q[3];
sx q[3];
rz(-1.1185631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23325486) q[2];
sx q[2];
rz(-0.7408064) q[2];
sx q[2];
rz(1.6439269) q[2];
rz(-0.46105591) q[3];
sx q[3];
rz(-2.4888829) q[3];
sx q[3];
rz(1.3155931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6719565) q[0];
sx q[0];
rz(-0.75103432) q[0];
sx q[0];
rz(-2.0627956) q[0];
rz(-0.59393334) q[1];
sx q[1];
rz(-2.1682231) q[1];
sx q[1];
rz(2.914391) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7861601) q[0];
sx q[0];
rz(-0.18635145) q[0];
sx q[0];
rz(0.92599286) q[0];
rz(-1.960177) q[2];
sx q[2];
rz(-2.7837278) q[2];
sx q[2];
rz(0.030454446) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46887423) q[1];
sx q[1];
rz(-1.0212082) q[1];
sx q[1];
rz(-0.36512111) q[1];
rz(-pi) q[2];
x q[2];
rz(3.041275) q[3];
sx q[3];
rz(-1.9411818) q[3];
sx q[3];
rz(-2.5357192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5220773) q[2];
sx q[2];
rz(-0.41831133) q[2];
sx q[2];
rz(2.412001) q[2];
rz(1.6884165) q[3];
sx q[3];
rz(-1.6875234) q[3];
sx q[3];
rz(0.61354536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59298092) q[0];
sx q[0];
rz(-2.225003) q[0];
sx q[0];
rz(2.7899637) q[0];
rz(2.8278606) q[1];
sx q[1];
rz(-1.2064563) q[1];
sx q[1];
rz(1.3765913) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54395318) q[0];
sx q[0];
rz(-2.7386129) q[0];
sx q[0];
rz(3.0209994) q[0];
rz(-pi) q[1];
rz(-0.11560924) q[2];
sx q[2];
rz(-1.7358276) q[2];
sx q[2];
rz(-0.045265667) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2089858) q[1];
sx q[1];
rz(-2.5521005) q[1];
sx q[1];
rz(-2.8578651) q[1];
rz(-1.2010716) q[3];
sx q[3];
rz(-1.4784357) q[3];
sx q[3];
rz(-0.28932387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54479638) q[2];
sx q[2];
rz(-0.72303855) q[2];
sx q[2];
rz(-1.5210305) q[2];
rz(2.9936786) q[3];
sx q[3];
rz(-1.3444129) q[3];
sx q[3];
rz(1.3676876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3677597) q[0];
sx q[0];
rz(-1.2667043) q[0];
sx q[0];
rz(1.2441147) q[0];
rz(-1.7077839) q[1];
sx q[1];
rz(-1.5807187) q[1];
sx q[1];
rz(2.921385) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6235636) q[0];
sx q[0];
rz(-1.9401778) q[0];
sx q[0];
rz(0.9250888) q[0];
x q[1];
rz(2.037165) q[2];
sx q[2];
rz(-1.4254071) q[2];
sx q[2];
rz(0.25823247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9088194) q[1];
sx q[1];
rz(-0.72130433) q[1];
sx q[1];
rz(0.90296332) q[1];
rz(-pi) q[2];
rz(0.098252849) q[3];
sx q[3];
rz(-1.6643401) q[3];
sx q[3];
rz(-0.52631179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9571017) q[2];
sx q[2];
rz(-1.9944921) q[2];
sx q[2];
rz(-1.5000878) q[2];
rz(-3.0575276) q[3];
sx q[3];
rz(-1.2596687) q[3];
sx q[3];
rz(3.0715004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641007) q[0];
sx q[0];
rz(-2.3832432) q[0];
sx q[0];
rz(0.41859928) q[0];
rz(-0.94340008) q[1];
sx q[1];
rz(-1.6523596) q[1];
sx q[1];
rz(-2.8387866) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3269791) q[0];
sx q[0];
rz(-1.7075065) q[0];
sx q[0];
rz(1.7390521) q[0];
rz(2.3934028) q[2];
sx q[2];
rz(-2.8837969) q[2];
sx q[2];
rz(-0.92743826) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9040519) q[1];
sx q[1];
rz(-1.6142323) q[1];
sx q[1];
rz(-3.0588845) q[1];
x q[2];
rz(-1.2262086) q[3];
sx q[3];
rz(-1.4610963) q[3];
sx q[3];
rz(-1.4668224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4960949) q[2];
sx q[2];
rz(-0.87254137) q[2];
sx q[2];
rz(0.75054753) q[2];
rz(-1.6802855) q[3];
sx q[3];
rz(-0.76992005) q[3];
sx q[3];
rz(-0.51699483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208329) q[0];
sx q[0];
rz(-1.5626386) q[0];
sx q[0];
rz(1.563969) q[0];
rz(1.9137406) q[1];
sx q[1];
rz(-0.49929437) q[1];
sx q[1];
rz(1.1566537) q[1];
rz(0.56671178) q[2];
sx q[2];
rz(-1.7253582) q[2];
sx q[2];
rz(2.2127989) q[2];
rz(1.0629366) q[3];
sx q[3];
rz(-2.776317) q[3];
sx q[3];
rz(-1.1767514) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
