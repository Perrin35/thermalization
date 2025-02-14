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
rz(-0.1470546) q[0];
sx q[0];
rz(-1.4015863) q[0];
sx q[0];
rz(0.92010486) q[0];
rz(-0.080634557) q[1];
sx q[1];
rz(-2.562685) q[1];
sx q[1];
rz(0.97715598) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.282519) q[0];
sx q[0];
rz(-0.4319829) q[0];
sx q[0];
rz(0.52406128) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0394568) q[2];
sx q[2];
rz(-0.13291026) q[2];
sx q[2];
rz(-1.5166211) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0339573) q[1];
sx q[1];
rz(-2.1672241) q[1];
sx q[1];
rz(1.7340842) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68859921) q[3];
sx q[3];
rz(-1.3312201) q[3];
sx q[3];
rz(1.857615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(-2.5207632) q[2];
rz(-2.6260455) q[3];
sx q[3];
rz(-1.7190869) q[3];
sx q[3];
rz(-2.2990885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2466549) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(2.5625693) q[0];
rz(-0.82194263) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(0.45375219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61868405) q[0];
sx q[0];
rz(-2.4755619) q[0];
sx q[0];
rz(0.36938195) q[0];
rz(-1.4582107) q[2];
sx q[2];
rz(-2.0346918) q[2];
sx q[2];
rz(0.91700208) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7612517) q[1];
sx q[1];
rz(-0.56219343) q[1];
sx q[1];
rz(2.9060049) q[1];
rz(1.0250799) q[3];
sx q[3];
rz(-1.9580055) q[3];
sx q[3];
rz(-2.6301165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4072676) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(-2.5767051) q[2];
rz(0.35483739) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(-1.7179276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9722209) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(2.5900904) q[0];
rz(-2.7768199) q[1];
sx q[1];
rz(-2.8021937) q[1];
sx q[1];
rz(0.50484467) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0802949) q[0];
sx q[0];
rz(-0.43879959) q[0];
sx q[0];
rz(0.3831692) q[0];
rz(1.4771492) q[2];
sx q[2];
rz(-1.159707) q[2];
sx q[2];
rz(-1.1716441) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6712196) q[1];
sx q[1];
rz(-0.80802901) q[1];
sx q[1];
rz(-0.68617448) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1944207) q[3];
sx q[3];
rz(-1.1191812) q[3];
sx q[3];
rz(-2.268689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43368936) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(-0.047253963) q[2];
rz(-0.83682483) q[3];
sx q[3];
rz(-2.4182726) q[3];
sx q[3];
rz(1.0251454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4243917) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(1.8151872) q[0];
rz(1.4855509) q[1];
sx q[1];
rz(-1.0842208) q[1];
sx q[1];
rz(2.5066689) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0795006) q[0];
sx q[0];
rz(-1.5720815) q[0];
sx q[0];
rz(1.5527524) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35648326) q[2];
sx q[2];
rz(-2.0429789) q[2];
sx q[2];
rz(3.0038578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0628218) q[1];
sx q[1];
rz(-0.34343849) q[1];
sx q[1];
rz(1.354864) q[1];
x q[2];
rz(-0.62795378) q[3];
sx q[3];
rz(-2.0139004) q[3];
sx q[3];
rz(-1.7114044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9012458) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(-1.3725613) q[2];
rz(-1.9339804) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(-2.8670368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32886252) q[0];
sx q[0];
rz(-0.48506081) q[0];
sx q[0];
rz(-0.10511705) q[0];
rz(2.8016727) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(2.0786659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8369097) q[0];
sx q[0];
rz(-1.5848918) q[0];
sx q[0];
rz(-2.430116) q[0];
x q[1];
rz(-1.8196443) q[2];
sx q[2];
rz(-0.82468678) q[2];
sx q[2];
rz(-0.4363554) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.868456) q[1];
sx q[1];
rz(-1.5850164) q[1];
sx q[1];
rz(-0.058560024) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5664728) q[3];
sx q[3];
rz(-2.6427442) q[3];
sx q[3];
rz(2.3259142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0792599) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(1.3366535) q[2];
rz(0.054232728) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(-0.59468734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.9688251) q[0];
sx q[0];
rz(-0.69197881) q[0];
sx q[0];
rz(-1.927595) q[0];
rz(2.0955775) q[1];
sx q[1];
rz(-1.0449301) q[1];
sx q[1];
rz(-2.2878343) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1778422) q[0];
sx q[0];
rz(-0.14687777) q[0];
sx q[0];
rz(-1.9903723) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5327024) q[2];
sx q[2];
rz(-0.85431803) q[2];
sx q[2];
rz(2.4520055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5484838) q[1];
sx q[1];
rz(-0.9973155) q[1];
sx q[1];
rz(-2.4683263) q[1];
rz(-1.4183956) q[3];
sx q[3];
rz(-1.9160685) q[3];
sx q[3];
rz(-2.4403646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13009109) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(-2.3940274) q[2];
rz(0.93303624) q[3];
sx q[3];
rz(-1.3676164) q[3];
sx q[3];
rz(-2.8396377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.1099243) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(2.5001496) q[0];
rz(2.172442) q[1];
sx q[1];
rz(-2.0717924) q[1];
sx q[1];
rz(-2.0279121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4232491) q[0];
sx q[0];
rz(-1.166543) q[0];
sx q[0];
rz(0.1057616) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53345726) q[2];
sx q[2];
rz(-2.1367624) q[2];
sx q[2];
rz(0.12803687) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42150527) q[1];
sx q[1];
rz(-2.0893789) q[1];
sx q[1];
rz(2.9353013) q[1];
rz(-pi) q[2];
rz(1.6051172) q[3];
sx q[3];
rz(-2.0908818) q[3];
sx q[3];
rz(1.8994562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.574719) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(-2.3835772) q[2];
rz(-0.30592439) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(-0.44736403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4161943) q[0];
sx q[0];
rz(-2.5328126) q[0];
sx q[0];
rz(-0.7533657) q[0];
rz(-2.2196409) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(0.13455483) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20924231) q[0];
sx q[0];
rz(-1.9268231) q[0];
sx q[0];
rz(0.092445157) q[0];
rz(-pi) q[1];
rz(-0.64132787) q[2];
sx q[2];
rz(-1.7761163) q[2];
sx q[2];
rz(-0.63831282) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8953029) q[1];
sx q[1];
rz(-3.0498235) q[1];
sx q[1];
rz(0.20905881) q[1];
x q[2];
rz(-0.67340019) q[3];
sx q[3];
rz(-1.9385425) q[3];
sx q[3];
rz(-1.0750225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3677463) q[2];
sx q[2];
rz(-1.4464804) q[2];
sx q[2];
rz(1.9473677) q[2];
rz(1.9959244) q[3];
sx q[3];
rz(-1.1402036) q[3];
sx q[3];
rz(2.0755419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1016178) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(-0.38247821) q[0];
rz(0.76599145) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(-1.3409748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4744953) q[0];
sx q[0];
rz(-1.7303559) q[0];
sx q[0];
rz(-1.480353) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4828311) q[2];
sx q[2];
rz(-0.81552699) q[2];
sx q[2];
rz(1.2620827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7786918) q[1];
sx q[1];
rz(-2.3039989) q[1];
sx q[1];
rz(1.2630839) q[1];
rz(-0.082499113) q[3];
sx q[3];
rz(-1.9256221) q[3];
sx q[3];
rz(0.88731836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0880903) q[2];
sx q[2];
rz(-0.92038766) q[2];
sx q[2];
rz(0.11605334) q[2];
rz(1.8807489) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(1.457823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5905404) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(0.1846479) q[0];
rz(0.91839904) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(1.437423) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6708095) q[0];
sx q[0];
rz(-0.52007404) q[0];
sx q[0];
rz(1.9896247) q[0];
rz(-pi) q[1];
rz(-1.499648) q[2];
sx q[2];
rz(-0.84991036) q[2];
sx q[2];
rz(-0.80903731) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9574163) q[1];
sx q[1];
rz(-0.97089689) q[1];
sx q[1];
rz(0.23317144) q[1];
rz(0.67622306) q[3];
sx q[3];
rz(-1.8340602) q[3];
sx q[3];
rz(-0.77419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.47716466) q[2];
sx q[2];
rz(-1.2457341) q[2];
sx q[2];
rz(2.5028382) q[2];
rz(-0.81513682) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(-2.7208929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.7246134) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(0.78085113) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(3.0197418) q[2];
sx q[2];
rz(-1.8486628) q[2];
sx q[2];
rz(1.8589414) q[2];
rz(-1.609997) q[3];
sx q[3];
rz(-1.8146252) q[3];
sx q[3];
rz(-1.6011325) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
