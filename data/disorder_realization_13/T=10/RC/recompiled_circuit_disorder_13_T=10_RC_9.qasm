OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9632602) q[0];
sx q[0];
rz(-1.652521) q[0];
sx q[0];
rz(0.89515495) q[0];
rz(-0.31495467) q[1];
sx q[1];
rz(-2.0575674) q[1];
sx q[1];
rz(-1.6853583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47269868) q[0];
sx q[0];
rz(-1.5905252) q[0];
sx q[0];
rz(1.701645) q[0];
rz(-2.1158754) q[2];
sx q[2];
rz(-1.8946049) q[2];
sx q[2];
rz(-0.88534249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97988843) q[1];
sx q[1];
rz(-0.6941422) q[1];
sx q[1];
rz(2.90467) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5159357) q[3];
sx q[3];
rz(-1.8901955) q[3];
sx q[3];
rz(-1.1064135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3866117) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(0.45271978) q[2];
rz(-0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2125856) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(-1.1520977) q[0];
rz(1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(2.6706085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5841056) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(-1.6886061) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2330301) q[2];
sx q[2];
rz(-2.517189) q[2];
sx q[2];
rz(2.117346) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2742548) q[1];
sx q[1];
rz(-1.5161637) q[1];
sx q[1];
rz(1.0824624) q[1];
rz(3.0456411) q[3];
sx q[3];
rz(-2.4437685) q[3];
sx q[3];
rz(0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(-0.2557959) q[2];
rz(-1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(-0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8783022) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(-0.27134744) q[0];
rz(0.73633206) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(0.43930611) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96734756) q[0];
sx q[0];
rz(-0.085701533) q[0];
sx q[0];
rz(1.218319) q[0];
x q[1];
rz(1.5191684) q[2];
sx q[2];
rz(-1.2106967) q[2];
sx q[2];
rz(-2.1215631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9648793) q[1];
sx q[1];
rz(-0.63767725) q[1];
sx q[1];
rz(1.8646851) q[1];
x q[2];
rz(2.1268232) q[3];
sx q[3];
rz(-1.1960256) q[3];
sx q[3];
rz(2.8127363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1674041) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(0.75508562) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077483594) q[0];
sx q[0];
rz(-1.5492726) q[0];
sx q[0];
rz(-1.2444929) q[0];
rz(0.81047932) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(-2.2669852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65507209) q[0];
sx q[0];
rz(-0.97210303) q[0];
sx q[0];
rz(-0.18874164) q[0];
rz(2.7565016) q[2];
sx q[2];
rz(-1.4355806) q[2];
sx q[2];
rz(-1.9584292) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8352656) q[1];
sx q[1];
rz(-1.7976465) q[1];
sx q[1];
rz(-0.38362417) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10505418) q[3];
sx q[3];
rz(-1.4964364) q[3];
sx q[3];
rz(-0.33945938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0041634) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(-1.7144263) q[2];
rz(-3.0754722) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(-1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(-2.5710035) q[0];
rz(2.5866306) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(2.3983009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63440454) q[0];
sx q[0];
rz(-2.2162063) q[0];
sx q[0];
rz(-2.8651644) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18851738) q[2];
sx q[2];
rz(-0.50893116) q[2];
sx q[2];
rz(2.237042) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1764644) q[1];
sx q[1];
rz(-1.7904518) q[1];
sx q[1];
rz(-1.3419282) q[1];
x q[2];
rz(2.8652142) q[3];
sx q[3];
rz(-0.82377269) q[3];
sx q[3];
rz(-0.36229047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9479998) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(1.0466446) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(-2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.627581) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(-2.6348689) q[0];
rz(2.8880033) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(-2.0862897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2282288) q[0];
sx q[0];
rz(-0.95634395) q[0];
sx q[0];
rz(-2.9173304) q[0];
rz(0.52168092) q[2];
sx q[2];
rz(-0.98419596) q[2];
sx q[2];
rz(2.4751543) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9195337) q[1];
sx q[1];
rz(-2.2026081) q[1];
sx q[1];
rz(2.6421089) q[1];
rz(-2.3985732) q[3];
sx q[3];
rz(-1.7457186) q[3];
sx q[3];
rz(-2.7878441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6598597) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(-0.8824904) q[2];
rz(-2.4957538) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(-1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(2.3690467) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(0.57377446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086512452) q[0];
sx q[0];
rz(-1.7614613) q[0];
sx q[0];
rz(3.029682) q[0];
rz(3.119486) q[2];
sx q[2];
rz(-1.7709641) q[2];
sx q[2];
rz(-2.9526763) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37807235) q[1];
sx q[1];
rz(-1.0135279) q[1];
sx q[1];
rz(-2.4850363) q[1];
rz(-pi) q[2];
rz(1.8049559) q[3];
sx q[3];
rz(-1.7171211) q[3];
sx q[3];
rz(-0.22406604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1022169) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8687826) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(1.1707206) q[0];
rz(-2.6314578) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(-1.8458813) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8745236) q[0];
sx q[0];
rz(-1.1133615) q[0];
sx q[0];
rz(1.9382856) q[0];
rz(-pi) q[1];
rz(-1.2175351) q[2];
sx q[2];
rz(-1.3152939) q[2];
sx q[2];
rz(-3.0614292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1295373) q[1];
sx q[1];
rz(-0.66857282) q[1];
sx q[1];
rz(1.936391) q[1];
rz(-pi) q[2];
rz(1.6298953) q[3];
sx q[3];
rz(-2.700138) q[3];
sx q[3];
rz(0.7837226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43508139) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(-0.56813017) q[2];
rz(2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-2.9476416) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(-2.4639159) q[0];
rz(2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(-2.303404) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4351589) q[0];
sx q[0];
rz(-2.7088232) q[0];
sx q[0];
rz(-0.5612527) q[0];
rz(-0.81460641) q[2];
sx q[2];
rz(-1.6199154) q[2];
sx q[2];
rz(-2.5533822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2182525) q[1];
sx q[1];
rz(-2.372101) q[1];
sx q[1];
rz(1.8583276) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30672726) q[3];
sx q[3];
rz(-2.8711651) q[3];
sx q[3];
rz(-0.21817792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80205408) q[2];
sx q[2];
rz(-2.4513117) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(0.96261111) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(-0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72120136) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(-2.904073) q[0];
rz(2.1233842) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(2.9097897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7983539) q[0];
sx q[0];
rz(-2.624247) q[0];
sx q[0];
rz(1.377064) q[0];
rz(-pi) q[1];
rz(-0.78328697) q[2];
sx q[2];
rz(-2.8056393) q[2];
sx q[2];
rz(-0.17620262) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3540436) q[1];
sx q[1];
rz(-1.7122867) q[1];
sx q[1];
rz(-2.8741807) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2139123) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(-3.1090581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(-2.7534289) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(-2.3378519) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5768455) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(0.3248365) q[2];
sx q[2];
rz(-1.3424716) q[2];
sx q[2];
rz(-2.9906103) q[2];
rz(0.10235056) q[3];
sx q[3];
rz(-2.9526738) q[3];
sx q[3];
rz(-1.916193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
