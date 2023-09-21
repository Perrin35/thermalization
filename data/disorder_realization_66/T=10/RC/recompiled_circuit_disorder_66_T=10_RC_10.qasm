OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(4.1132676) q[0];
sx q[0];
rz(10.905807) q[0];
rz(3.0193168) q[1];
sx q[1];
rz(-3.0552157) q[1];
sx q[1];
rz(3.1187305) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1296596) q[0];
sx q[0];
rz(-1.4225905) q[0];
sx q[0];
rz(0.11797842) q[0];
rz(-pi) q[1];
rz(2.4542698) q[2];
sx q[2];
rz(-0.63185531) q[2];
sx q[2];
rz(-0.99422115) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.80427985) q[1];
sx q[1];
rz(-2.7645281) q[1];
sx q[1];
rz(-1.1537329) q[1];
rz(2.888527) q[3];
sx q[3];
rz(-0.96732891) q[3];
sx q[3];
rz(-0.30482182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7227398) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(2.9700759) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(-0.086659327) q[0];
rz(0.86241972) q[1];
sx q[1];
rz(-0.54310596) q[1];
sx q[1];
rz(-0.08509732) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9239685) q[0];
sx q[0];
rz(-0.94857615) q[0];
sx q[0];
rz(-1.5505962) q[0];
x q[1];
rz(-0.30479635) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(-0.09588974) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7460691) q[1];
sx q[1];
rz(-1.4847401) q[1];
sx q[1];
rz(1.6834016) q[1];
rz(-pi) q[2];
rz(-2.0410791) q[3];
sx q[3];
rz(-1.2840052) q[3];
sx q[3];
rz(1.4834529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.4651728) q[2];
rz(1.9654467) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(2.235967) q[0];
rz(2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.3844301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6838283) q[0];
sx q[0];
rz(-1.0965075) q[0];
sx q[0];
rz(-1.6636687) q[0];
rz(-pi) q[1];
x q[1];
rz(2.076782) q[2];
sx q[2];
rz(-1.9188606) q[2];
sx q[2];
rz(-2.6315174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.78450173) q[1];
sx q[1];
rz(-0.75575268) q[1];
sx q[1];
rz(3.109039) q[1];
rz(0.82377394) q[3];
sx q[3];
rz(-0.71411055) q[3];
sx q[3];
rz(-2.2658474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(-1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06552799) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(-1.7975851) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(2.9342594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5642693) q[0];
sx q[0];
rz(-3.0992357) q[0];
sx q[0];
rz(-1.7839412) q[0];
rz(-3.1111654) q[2];
sx q[2];
rz(-1.7446339) q[2];
sx q[2];
rz(-2.0174842) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.281665) q[1];
sx q[1];
rz(-1.7237701) q[1];
sx q[1];
rz(-1.8106736) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5563577) q[3];
sx q[3];
rz(-2.0824021) q[3];
sx q[3];
rz(0.59359854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8308782) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-2.3332398) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(-0.19367735) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(2.2461058) q[0];
rz(-2.8764309) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-0.53363824) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14347178) q[0];
sx q[0];
rz(-0.40092418) q[0];
sx q[0];
rz(0.801416) q[0];
rz(-pi) q[1];
rz(2.9075378) q[2];
sx q[2];
rz(-2.3599527) q[2];
sx q[2];
rz(-0.78125886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67709778) q[1];
sx q[1];
rz(-2.8765196) q[1];
sx q[1];
rz(-2.5527843) q[1];
rz(2.93612) q[3];
sx q[3];
rz(-0.95147248) q[3];
sx q[3];
rz(-1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-2.5732102) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50826532) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(1.5559224) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(1.0046545) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6509103) q[0];
sx q[0];
rz(-2.6255529) q[0];
sx q[0];
rz(2.3923621) q[0];
rz(1.651628) q[2];
sx q[2];
rz(-1.9664552) q[2];
sx q[2];
rz(0.8880907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1130044) q[1];
sx q[1];
rz(-2.5341946) q[1];
sx q[1];
rz(0.44087704) q[1];
rz(-pi) q[2];
rz(1.1569571) q[3];
sx q[3];
rz(-1.5598179) q[3];
sx q[3];
rz(-1.699284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4040318) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(-0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(-1.8572042) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(2.4694494) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7710966) q[0];
sx q[0];
rz(-1.5633977) q[0];
sx q[0];
rz(-1.7367944) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0124112) q[2];
sx q[2];
rz(-1.1814983) q[2];
sx q[2];
rz(-2.0814975) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0828515) q[1];
sx q[1];
rz(-2.5141659) q[1];
sx q[1];
rz(-0.65545603) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70127212) q[3];
sx q[3];
rz(-1.6810732) q[3];
sx q[3];
rz(0.14227223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(0.82459015) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(-0.59857541) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44052112) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(-0.80668443) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0643443) q[0];
sx q[0];
rz(-1.4532538) q[0];
sx q[0];
rz(2.006152) q[0];
x q[1];
rz(1.6786472) q[2];
sx q[2];
rz(-1.3662405) q[2];
sx q[2];
rz(-2.4754935) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.403277) q[1];
sx q[1];
rz(-2.3936831) q[1];
sx q[1];
rz(0.43055375) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8716378) q[3];
sx q[3];
rz(-1.6013718) q[3];
sx q[3];
rz(1.9105063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5640101) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(-3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-2.0172393) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(0.28655562) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(2.8009169) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17250241) q[0];
sx q[0];
rz(-1.889125) q[0];
sx q[0];
rz(-0.83842917) q[0];
x q[1];
rz(-1.9129487) q[2];
sx q[2];
rz(-2.5492382) q[2];
sx q[2];
rz(2.5387788) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77010158) q[1];
sx q[1];
rz(-1.0189459) q[1];
sx q[1];
rz(2.9279207) q[1];
rz(-2.4386028) q[3];
sx q[3];
rz(-0.25745108) q[3];
sx q[3];
rz(-1.6182181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75491536) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(2.6436451) q[2];
rz(2.8399816) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(-0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122445) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(0.27858946) q[0];
rz(2.5623698) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(-0.07671193) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0522703) q[0];
sx q[0];
rz(-0.33903402) q[0];
sx q[0];
rz(-1.9498755) q[0];
rz(-0.17148359) q[2];
sx q[2];
rz(-2.137261) q[2];
sx q[2];
rz(-0.4085853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0127276) q[1];
sx q[1];
rz(-2.7755133) q[1];
sx q[1];
rz(3.0831343) q[1];
rz(1.65927) q[3];
sx q[3];
rz(-1.7469329) q[3];
sx q[3];
rz(2.4000771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-2.6860766) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(-0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-2.5554399) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(-1.0473245) q[2];
sx q[2];
rz(-2.6101255) q[2];
sx q[2];
rz(2.5892467) q[2];
rz(-2.7145731) q[3];
sx q[3];
rz(-1.9058766) q[3];
sx q[3];
rz(1.8591892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];