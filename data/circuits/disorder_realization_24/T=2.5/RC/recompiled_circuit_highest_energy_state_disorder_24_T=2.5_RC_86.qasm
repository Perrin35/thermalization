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
rz(2.6296122) q[0];
sx q[0];
rz(2.5706302) q[0];
sx q[0];
rz(9.6124967) q[0];
rz(0.81387782) q[1];
sx q[1];
rz(-1.6698281) q[1];
sx q[1];
rz(-1.8522813) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16351249) q[0];
sx q[0];
rz(-2.8239282) q[0];
sx q[0];
rz(-2.4207522) q[0];
x q[1];
rz(-2.5709349) q[2];
sx q[2];
rz(-1.5458115) q[2];
sx q[2];
rz(-2.5236584) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6481406) q[1];
sx q[1];
rz(-1.3228184) q[1];
sx q[1];
rz(0.47391717) q[1];
rz(-0.23203316) q[3];
sx q[3];
rz(-1.9153144) q[3];
sx q[3];
rz(-1.0330878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2671555) q[2];
sx q[2];
rz(-2.0384553) q[2];
sx q[2];
rz(-0.77679408) q[2];
rz(1.3022425) q[3];
sx q[3];
rz(-1.4812508) q[3];
sx q[3];
rz(-1.9384025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(2.5341107) q[0];
sx q[0];
rz(-2.6597839) q[0];
sx q[0];
rz(-0.99824655) q[0];
rz(2.7718995) q[1];
sx q[1];
rz(-1.0414711) q[1];
sx q[1];
rz(-2.0236156) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6551483) q[0];
sx q[0];
rz(-2.110743) q[0];
sx q[0];
rz(2.6686431) q[0];
x q[1];
rz(-2.4547365) q[2];
sx q[2];
rz(-1.1026898) q[2];
sx q[2];
rz(1.3324225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0675903) q[1];
sx q[1];
rz(-2.6704881) q[1];
sx q[1];
rz(-1.2420688) q[1];
rz(1.9309773) q[3];
sx q[3];
rz(-0.66347117) q[3];
sx q[3];
rz(-1.0960032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2758241) q[2];
sx q[2];
rz(-1.3752702) q[2];
sx q[2];
rz(2.7562874) q[2];
rz(3.1028808) q[3];
sx q[3];
rz(-2.7888515) q[3];
sx q[3];
rz(-0.64046162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5019219) q[0];
sx q[0];
rz(-2.6646035) q[0];
sx q[0];
rz(-1.3013526) q[0];
rz(3.0838857) q[1];
sx q[1];
rz(-0.4464018) q[1];
sx q[1];
rz(-1.3267964) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0327137) q[0];
sx q[0];
rz(-1.098806) q[0];
sx q[0];
rz(2.1817939) q[0];
rz(0.16043255) q[2];
sx q[2];
rz(-1.4331054) q[2];
sx q[2];
rz(-1.4169803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88736278) q[1];
sx q[1];
rz(-1.0471724) q[1];
sx q[1];
rz(-1.9033405) q[1];
x q[2];
rz(-0.84872021) q[3];
sx q[3];
rz(-1.8406788) q[3];
sx q[3];
rz(3.0218642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79746276) q[2];
sx q[2];
rz(-2.3064488) q[2];
sx q[2];
rz(2.3878035) q[2];
rz(-1.7720743) q[3];
sx q[3];
rz(-1.6794208) q[3];
sx q[3];
rz(2.4863825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264483) q[0];
sx q[0];
rz(-2.0560052) q[0];
sx q[0];
rz(-1.7468859) q[0];
rz(-1.9649547) q[1];
sx q[1];
rz(-1.6329012) q[1];
sx q[1];
rz(-2.7746157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1934051) q[0];
sx q[0];
rz(-0.087554878) q[0];
sx q[0];
rz(2.2746536) q[0];
rz(-2.8072678) q[2];
sx q[2];
rz(-1.3204638) q[2];
sx q[2];
rz(1.3368397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9865468) q[1];
sx q[1];
rz(-1.9268225) q[1];
sx q[1];
rz(1.3638391) q[1];
rz(2.9212679) q[3];
sx q[3];
rz(-0.37543618) q[3];
sx q[3];
rz(-2.8130949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41161141) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(1.8709987) q[2];
rz(2.1805084) q[3];
sx q[3];
rz(-1.4189439) q[3];
sx q[3];
rz(3.0911176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5807895) q[0];
sx q[0];
rz(-2.2152948) q[0];
sx q[0];
rz(1.8286937) q[0];
rz(1.987223) q[1];
sx q[1];
rz(-1.1416953) q[1];
sx q[1];
rz(-0.55783522) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2247705) q[0];
sx q[0];
rz(-2.0799412) q[0];
sx q[0];
rz(0.96352412) q[0];
rz(-pi) q[1];
rz(-1.4935605) q[2];
sx q[2];
rz(-1.7309411) q[2];
sx q[2];
rz(2.3004265) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6184885) q[1];
sx q[1];
rz(-2.2918211) q[1];
sx q[1];
rz(-0.40722653) q[1];
rz(-2.6799503) q[3];
sx q[3];
rz(-1.4351234) q[3];
sx q[3];
rz(1.5380579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4981726) q[2];
sx q[2];
rz(-1.8489445) q[2];
sx q[2];
rz(-0.84929973) q[2];
rz(2.9863206) q[3];
sx q[3];
rz(-0.97584358) q[3];
sx q[3];
rz(2.7707905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7840541) q[0];
sx q[0];
rz(-0.58670601) q[0];
sx q[0];
rz(2.6724755) q[0];
rz(-2.6672089) q[1];
sx q[1];
rz(-2.3553039) q[1];
sx q[1];
rz(-0.75538409) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82925382) q[0];
sx q[0];
rz(-2.8702998) q[0];
sx q[0];
rz(-1.5872699) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9017436) q[2];
sx q[2];
rz(-2.1734383) q[2];
sx q[2];
rz(1.861426) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1052221) q[1];
sx q[1];
rz(-1.3579988) q[1];
sx q[1];
rz(0.12963055) q[1];
x q[2];
rz(0.42879019) q[3];
sx q[3];
rz(-0.8800104) q[3];
sx q[3];
rz(0.59861983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0581806) q[2];
sx q[2];
rz(-1.8184793) q[2];
sx q[2];
rz(-0.18784909) q[2];
rz(1.8958873) q[3];
sx q[3];
rz(-1.3789504) q[3];
sx q[3];
rz(-2.0511621) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.93194) q[0];
sx q[0];
rz(-1.8745475) q[0];
sx q[0];
rz(-2.7506822) q[0];
rz(2.2115425) q[1];
sx q[1];
rz(-1.4314194) q[1];
sx q[1];
rz(1.8375058) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025607312) q[0];
sx q[0];
rz(-1.5918333) q[0];
sx q[0];
rz(1.9217092) q[0];
x q[1];
rz(-1.6565336) q[2];
sx q[2];
rz(-1.9112327) q[2];
sx q[2];
rz(2.1742333) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3871349) q[1];
sx q[1];
rz(-1.4031193) q[1];
sx q[1];
rz(-1.1134336) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.061873925) q[3];
sx q[3];
rz(-0.29682595) q[3];
sx q[3];
rz(0.59943953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8730674) q[2];
sx q[2];
rz(-2.8748942) q[2];
sx q[2];
rz(0.11338691) q[2];
rz(-3.125627) q[3];
sx q[3];
rz(-1.6458052) q[3];
sx q[3];
rz(2.9769843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575344) q[0];
sx q[0];
rz(-1.6679732) q[0];
sx q[0];
rz(3.0175324) q[0];
rz(-3.0198174) q[1];
sx q[1];
rz(-0.75141326) q[1];
sx q[1];
rz(1.6990936) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.34262) q[0];
sx q[0];
rz(-0.39726394) q[0];
sx q[0];
rz(1.308206) q[0];
rz(-pi) q[1];
rz(0.29192544) q[2];
sx q[2];
rz(-1.5291844) q[2];
sx q[2];
rz(-2.6440563) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5296764) q[1];
sx q[1];
rz(-1.9135336) q[1];
sx q[1];
rz(-0.14118282) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4115218) q[3];
sx q[3];
rz(-1.3833117) q[3];
sx q[3];
rz(1.9188251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.96559912) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(2.3967801) q[2];
rz(0.92721573) q[3];
sx q[3];
rz(-1.6330481) q[3];
sx q[3];
rz(1.0923227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6532779) q[0];
sx q[0];
rz(-1.7019685) q[0];
sx q[0];
rz(2.9162245) q[0];
rz(1.1124181) q[1];
sx q[1];
rz(-2.367159) q[1];
sx q[1];
rz(2.2427799) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8161801) q[0];
sx q[0];
rz(-1.0409779) q[0];
sx q[0];
rz(-2.2384032) q[0];
x q[1];
rz(0.65943879) q[2];
sx q[2];
rz(-2.5213679) q[2];
sx q[2];
rz(-2.8218215) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.73279335) q[1];
sx q[1];
rz(-2.2299754) q[1];
sx q[1];
rz(2.885347) q[1];
rz(1.8636673) q[3];
sx q[3];
rz(-1.5808269) q[3];
sx q[3];
rz(0.18126479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0972458) q[2];
sx q[2];
rz(-1.4298507) q[2];
sx q[2];
rz(-0.35935768) q[2];
rz(-2.8906631) q[3];
sx q[3];
rz(-1.072262) q[3];
sx q[3];
rz(2.3738764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74442416) q[0];
sx q[0];
rz(-3.011062) q[0];
sx q[0];
rz(0.8771483) q[0];
rz(1.5646704) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(0.78561479) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6495384) q[0];
sx q[0];
rz(-1.7888594) q[0];
sx q[0];
rz(1.8137003) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.142611) q[2];
sx q[2];
rz(-1.0773563) q[2];
sx q[2];
rz(-2.7827415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4686376) q[1];
sx q[1];
rz(-1.5002999) q[1];
sx q[1];
rz(-2.7967909) q[1];
rz(1.4203868) q[3];
sx q[3];
rz(-2.0724943) q[3];
sx q[3];
rz(3.0034163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2239573) q[2];
sx q[2];
rz(-0.79006299) q[2];
sx q[2];
rz(-2.4033974) q[2];
rz(-2.2186642) q[3];
sx q[3];
rz(-1.8332558) q[3];
sx q[3];
rz(-0.2963399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34687635) q[0];
sx q[0];
rz(-1.3608169) q[0];
sx q[0];
rz(0.97186744) q[0];
rz(-2.234266) q[1];
sx q[1];
rz(-1.0846039) q[1];
sx q[1];
rz(-0.32122282) q[1];
rz(-2.1512866) q[2];
sx q[2];
rz(-1.9038426) q[2];
sx q[2];
rz(2.7833084) q[2];
rz(-2.3774556) q[3];
sx q[3];
rz(-0.69307477) q[3];
sx q[3];
rz(0.57144036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
