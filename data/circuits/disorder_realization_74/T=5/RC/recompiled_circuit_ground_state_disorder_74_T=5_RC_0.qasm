OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57269078) q[0];
sx q[0];
rz(-2.1532018) q[0];
sx q[0];
rz(0.44901499) q[0];
rz(2.9721337) q[1];
sx q[1];
rz(-0.13154498) q[1];
sx q[1];
rz(2.0102672) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9207391) q[0];
sx q[0];
rz(-1.8628696) q[0];
sx q[0];
rz(0.9432442) q[0];
x q[1];
rz(-0.65931084) q[2];
sx q[2];
rz(-0.74987312) q[2];
sx q[2];
rz(1.4091968) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4428569) q[1];
sx q[1];
rz(-1.3971824) q[1];
sx q[1];
rz(0.97877494) q[1];
x q[2];
rz(-1.851397) q[3];
sx q[3];
rz(-2.5803356) q[3];
sx q[3];
rz(3.0292976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1863056) q[2];
sx q[2];
rz(-0.39262843) q[2];
sx q[2];
rz(0.10360959) q[2];
rz(-2.7747532) q[3];
sx q[3];
rz(-1.6678526) q[3];
sx q[3];
rz(-1.8240671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9570626) q[0];
sx q[0];
rz(-2.6657031) q[0];
sx q[0];
rz(2.6153508) q[0];
rz(-1.2435675) q[1];
sx q[1];
rz(-1.7310111) q[1];
sx q[1];
rz(-2.9454561) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.552452) q[0];
sx q[0];
rz(-1.6375443) q[0];
sx q[0];
rz(1.1132973) q[0];
rz(-2.4320658) q[2];
sx q[2];
rz(-1.6975743) q[2];
sx q[2];
rz(-1.4122054) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6013697) q[1];
sx q[1];
rz(-1.1246846) q[1];
sx q[1];
rz(-0.15769728) q[1];
rz(-3.1297471) q[3];
sx q[3];
rz(-1.5812173) q[3];
sx q[3];
rz(1.1529779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5474995) q[2];
sx q[2];
rz(-0.27442351) q[2];
sx q[2];
rz(-0.37626949) q[2];
rz(-2.2606692) q[3];
sx q[3];
rz(-1.8062402) q[3];
sx q[3];
rz(0.81203619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45102099) q[0];
sx q[0];
rz(-2.8132827) q[0];
sx q[0];
rz(-2.1300533) q[0];
rz(-0.66863376) q[1];
sx q[1];
rz(-1.1625544) q[1];
sx q[1];
rz(0.54723251) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035346383) q[0];
sx q[0];
rz(-1.5900471) q[0];
sx q[0];
rz(-2.9401642) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7906495) q[2];
sx q[2];
rz(-1.8230652) q[2];
sx q[2];
rz(2.6961435) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4302252) q[1];
sx q[1];
rz(-2.2177296) q[1];
sx q[1];
rz(2.837846) q[1];
rz(-0.8143592) q[3];
sx q[3];
rz(-1.4297155) q[3];
sx q[3];
rz(0.94179487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68613595) q[2];
sx q[2];
rz(-2.8563209) q[2];
sx q[2];
rz(1.6395462) q[2];
rz(2.5655668) q[3];
sx q[3];
rz(-1.6582158) q[3];
sx q[3];
rz(-2.7801133) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6308924) q[0];
sx q[0];
rz(-2.3402813) q[0];
sx q[0];
rz(-2.8019688) q[0];
rz(-2.6009808) q[1];
sx q[1];
rz(-2.4410591) q[1];
sx q[1];
rz(-2.2115754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7226573) q[0];
sx q[0];
rz(-1.5003073) q[0];
sx q[0];
rz(-1.7138636) q[0];
rz(-pi) q[1];
rz(-2.5901718) q[2];
sx q[2];
rz(-2.4822771) q[2];
sx q[2];
rz(-2.2267411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2423289) q[1];
sx q[1];
rz(-1.5643969) q[1];
sx q[1];
rz(2.0501627) q[1];
rz(-pi) q[2];
rz(-2.2100979) q[3];
sx q[3];
rz(-2.3911797) q[3];
sx q[3];
rz(-1.9281333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0129619) q[2];
sx q[2];
rz(-1.7065115) q[2];
sx q[2];
rz(1.5677412) q[2];
rz(-2.3980127) q[3];
sx q[3];
rz(-0.79135528) q[3];
sx q[3];
rz(-0.96405205) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6745233) q[0];
sx q[0];
rz(-0.85551298) q[0];
sx q[0];
rz(-3.0849482) q[0];
rz(1.665834) q[1];
sx q[1];
rz(-1.1328127) q[1];
sx q[1];
rz(-1.3585565) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18902738) q[0];
sx q[0];
rz(-1.0882821) q[0];
sx q[0];
rz(0.58084647) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3681202) q[2];
sx q[2];
rz(-2.2182825) q[2];
sx q[2];
rz(1.7062239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92734071) q[1];
sx q[1];
rz(-2.3433629) q[1];
sx q[1];
rz(2.4386474) q[1];
x q[2];
rz(-0.18601619) q[3];
sx q[3];
rz(-2.5721492) q[3];
sx q[3];
rz(1.525477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.27611) q[2];
sx q[2];
rz(-1.7296187) q[2];
sx q[2];
rz(1.126368) q[2];
rz(-1.3364835) q[3];
sx q[3];
rz(-1.0765321) q[3];
sx q[3];
rz(-2.0520463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4514076) q[0];
sx q[0];
rz(-1.0253588) q[0];
sx q[0];
rz(3.0080646) q[0];
rz(-2.1639157) q[1];
sx q[1];
rz(-0.97594273) q[1];
sx q[1];
rz(-0.28688637) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7909951) q[0];
sx q[0];
rz(-2.4387601) q[0];
sx q[0];
rz(-1.4969456) q[0];
x q[1];
rz(1.8989766) q[2];
sx q[2];
rz(-1.6115453) q[2];
sx q[2];
rz(-0.894067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62090767) q[1];
sx q[1];
rz(-1.0093099) q[1];
sx q[1];
rz(-0.26446277) q[1];
rz(-pi) q[2];
x q[2];
rz(1.828015) q[3];
sx q[3];
rz(-2.1495499) q[3];
sx q[3];
rz(-0.46427765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7202683) q[2];
sx q[2];
rz(-1.4807533) q[2];
sx q[2];
rz(-1.486091) q[2];
rz(0.36763516) q[3];
sx q[3];
rz(-2.1143819) q[3];
sx q[3];
rz(-2.0462842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7169645) q[0];
sx q[0];
rz(-2.3724738) q[0];
sx q[0];
rz(0.47384438) q[0];
rz(-0.66420707) q[1];
sx q[1];
rz(-1.9240446) q[1];
sx q[1];
rz(-1.3833822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5055588) q[0];
sx q[0];
rz(-2.2657388) q[0];
sx q[0];
rz(2.4916584) q[0];
rz(-pi) q[1];
x q[1];
rz(1.152731) q[2];
sx q[2];
rz(-1.5454614) q[2];
sx q[2];
rz(-2.096614) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6205299) q[1];
sx q[1];
rz(-2.8179666) q[1];
sx q[1];
rz(0.77038295) q[1];
rz(-1.7960127) q[3];
sx q[3];
rz(-1.3966421) q[3];
sx q[3];
rz(0.98291493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.15709269) q[2];
sx q[2];
rz(-1.7326771) q[2];
sx q[2];
rz(-0.3248997) q[2];
rz(0.38069185) q[3];
sx q[3];
rz(-0.84563962) q[3];
sx q[3];
rz(1.0977753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0710881) q[0];
sx q[0];
rz(-0.66362137) q[0];
sx q[0];
rz(-0.93884236) q[0];
rz(-0.70010575) q[1];
sx q[1];
rz(-0.63799262) q[1];
sx q[1];
rz(-0.28900388) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5225713) q[0];
sx q[0];
rz(-0.24174015) q[0];
sx q[0];
rz(1.2416583) q[0];
rz(-pi) q[1];
rz(2.2174066) q[2];
sx q[2];
rz(-1.1614387) q[2];
sx q[2];
rz(-1.5034624) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6470312) q[1];
sx q[1];
rz(-0.53003487) q[1];
sx q[1];
rz(1.7575592) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6427755) q[3];
sx q[3];
rz(-1.9508617) q[3];
sx q[3];
rz(0.38004181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2603904) q[2];
sx q[2];
rz(-1.5217047) q[2];
sx q[2];
rz(2.728906) q[2];
rz(0.9203426) q[3];
sx q[3];
rz(-2.6350382) q[3];
sx q[3];
rz(-1.9119561) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5145787) q[0];
sx q[0];
rz(-2.161442) q[0];
sx q[0];
rz(-2.8421616) q[0];
rz(-2.0191655) q[1];
sx q[1];
rz(-1.9405245) q[1];
sx q[1];
rz(1.0312414) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64410463) q[0];
sx q[0];
rz(-1.9369164) q[0];
sx q[0];
rz(-3.1300504) q[0];
rz(-pi) q[1];
rz(2.6793807) q[2];
sx q[2];
rz(-1.3920857) q[2];
sx q[2];
rz(-1.4437087) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9523176) q[1];
sx q[1];
rz(-0.77042033) q[1];
sx q[1];
rz(-1.1483869) q[1];
rz(-0.93017756) q[3];
sx q[3];
rz(-2.0559466) q[3];
sx q[3];
rz(1.5826011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9791744) q[2];
sx q[2];
rz(-1.3434429) q[2];
sx q[2];
rz(0.24270414) q[2];
rz(1.8414712) q[3];
sx q[3];
rz(-1.0588812) q[3];
sx q[3];
rz(0.21568957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31371394) q[0];
sx q[0];
rz(-2.9254881) q[0];
sx q[0];
rz(1.7207654) q[0];
rz(2.8315262) q[1];
sx q[1];
rz(-2.2646751) q[1];
sx q[1];
rz(1.5975331) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9102073) q[0];
sx q[0];
rz(-0.69150309) q[0];
sx q[0];
rz(0.67759902) q[0];
rz(-pi) q[1];
rz(1.8158967) q[2];
sx q[2];
rz(-2.7900006) q[2];
sx q[2];
rz(1.0838255) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.032736445) q[1];
sx q[1];
rz(-1.8305155) q[1];
sx q[1];
rz(1.3459567) q[1];
x q[2];
rz(0.15940729) q[3];
sx q[3];
rz(-1.4324354) q[3];
sx q[3];
rz(2.0970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0073504) q[2];
sx q[2];
rz(-1.728936) q[2];
sx q[2];
rz(-1.7751815) q[2];
rz(0.8775231) q[3];
sx q[3];
rz(-1.4656504) q[3];
sx q[3];
rz(-2.2519978) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14961814) q[0];
sx q[0];
rz(-1.5852954) q[0];
sx q[0];
rz(-1.4854767) q[0];
rz(-0.86434518) q[1];
sx q[1];
rz(-2.1280011) q[1];
sx q[1];
rz(-0.43307532) q[1];
rz(2.6411459) q[2];
sx q[2];
rz(-1.3896349) q[2];
sx q[2];
rz(-1.096772) q[2];
rz(-1.6771562) q[3];
sx q[3];
rz(-0.48135664) q[3];
sx q[3];
rz(2.4873747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
