OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(-2.3119976) q[0];
sx q[0];
rz(-0.15396804) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(-2.1492465) q[1];
sx q[1];
rz(-0.33831236) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1725537) q[0];
sx q[0];
rz(-1.8882897) q[0];
sx q[0];
rz(-2.88455) q[0];
rz(-pi) q[1];
rz(-0.89262427) q[2];
sx q[2];
rz(-1.8632338) q[2];
sx q[2];
rz(2.6543648) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4813862) q[1];
sx q[1];
rz(-1.8567137) q[1];
sx q[1];
rz(1.6260765) q[1];
rz(3.1258718) q[3];
sx q[3];
rz(-1.058488) q[3];
sx q[3];
rz(-0.44954625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9989495) q[2];
sx q[2];
rz(-0.34037408) q[2];
sx q[2];
rz(-1.1738698) q[2];
rz(-0.075803444) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.8006111) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(3.0766292) q[0];
rz(0.57463542) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(1.2423135) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64105469) q[0];
sx q[0];
rz(-0.28028742) q[0];
sx q[0];
rz(2.6602402) q[0];
rz(-1.3556446) q[2];
sx q[2];
rz(-0.64385507) q[2];
sx q[2];
rz(1.7807963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6537135) q[1];
sx q[1];
rz(-2.5839845) q[1];
sx q[1];
rz(0.30456581) q[1];
x q[2];
rz(0.55450704) q[3];
sx q[3];
rz(-0.79075659) q[3];
sx q[3];
rz(-0.33695541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3339281) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(2.5675473) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(3.0103502) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42049256) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(0.72845355) q[0];
rz(1.4942253) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(-1.0167936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3166312) q[0];
sx q[0];
rz(-1.5350071) q[0];
sx q[0];
rz(0.081598452) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74480199) q[2];
sx q[2];
rz(-0.6664657) q[2];
sx q[2];
rz(1.9772066) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6779855) q[1];
sx q[1];
rz(-1.4336168) q[1];
sx q[1];
rz(-0.82171085) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4025027) q[3];
sx q[3];
rz(-0.28545359) q[3];
sx q[3];
rz(-1.1807549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3399405) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(-1.0495079) q[2];
rz(-0.63878757) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8124354) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(1.6695492) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-2.3627294) q[1];
sx q[1];
rz(0.24681117) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4287764) q[0];
sx q[0];
rz(-1.4220337) q[0];
sx q[0];
rz(2.2674198) q[0];
rz(-pi) q[1];
rz(-0.39101379) q[2];
sx q[2];
rz(-2.6303929) q[2];
sx q[2];
rz(2.9875987) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81891638) q[1];
sx q[1];
rz(-0.37441844) q[1];
sx q[1];
rz(-0.14426343) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28624268) q[3];
sx q[3];
rz(-2.1677368) q[3];
sx q[3];
rz(-2.3846574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6639158) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(-1.5412615) q[2];
rz(2.4345496) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33070579) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(-1.3274308) q[0];
rz(-1.56303) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(-0.24838233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7029593) q[0];
sx q[0];
rz(-2.6002433) q[0];
sx q[0];
rz(-0.066141733) q[0];
rz(2.9551198) q[2];
sx q[2];
rz(-1.7261793) q[2];
sx q[2];
rz(-1.1764256) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3779113) q[1];
sx q[1];
rz(-1.3994201) q[1];
sx q[1];
rz(-2.548449) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31110839) q[3];
sx q[3];
rz(-2.4690383) q[3];
sx q[3];
rz(-2.6274519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7053232) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(-2.752839) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(-0.50271547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0080863) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(1.3457993) q[0];
rz(-2.0603518) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(3.016901) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9799177) q[0];
sx q[0];
rz(-1.7470164) q[0];
sx q[0];
rz(-1.6582703) q[0];
rz(-pi) q[1];
rz(-0.6090392) q[2];
sx q[2];
rz(-1.3281203) q[2];
sx q[2];
rz(0.86953029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5387419) q[1];
sx q[1];
rz(-0.83487836) q[1];
sx q[1];
rz(1.49453) q[1];
rz(-pi) q[2];
rz(0.74495875) q[3];
sx q[3];
rz(-0.31674851) q[3];
sx q[3];
rz(1.484364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51320118) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(0.61895269) q[2];
rz(1.0533054) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(-0.9296023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597647) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(-1.0429617) q[0];
rz(-0.46328059) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(-1.0707062) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8541504) q[0];
sx q[0];
rz(-1.5702015) q[0];
sx q[0];
rz(2.6575412) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3473862) q[2];
sx q[2];
rz(-1.6958478) q[2];
sx q[2];
rz(-1.8206247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6003905) q[1];
sx q[1];
rz(-1.4711079) q[1];
sx q[1];
rz(0.77460918) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.058810874) q[3];
sx q[3];
rz(-2.8088514) q[3];
sx q[3];
rz(-2.9174093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36859194) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-2.8179742) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(-2.0543082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535646) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(1.9288829) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(0.94747296) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5652126) q[0];
sx q[0];
rz(-2.1064261) q[0];
sx q[0];
rz(3.0239848) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1768441) q[2];
sx q[2];
rz(-1.7272005) q[2];
sx q[2];
rz(-2.1795189) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.84711134) q[1];
sx q[1];
rz(-1.0080907) q[1];
sx q[1];
rz(-0.28114762) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34668563) q[3];
sx q[3];
rz(-2.3035435) q[3];
sx q[3];
rz(0.29688641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56132135) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(-2.6718111) q[2];
rz(-1.8404768) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(-0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8895421) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(-1.8126194) q[0];
rz(2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-2.9387617) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751574) q[0];
sx q[0];
rz(-0.046910722) q[0];
sx q[0];
rz(-0.58880083) q[0];
x q[1];
rz(-1.780026) q[2];
sx q[2];
rz(-1.3522569) q[2];
sx q[2];
rz(0.92313672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.030414) q[1];
sx q[1];
rz(-1.5215538) q[1];
sx q[1];
rz(-1.4700252) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6796954) q[3];
sx q[3];
rz(-1.4326722) q[3];
sx q[3];
rz(-1.4702597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41708502) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(-0.99651304) q[2];
rz(2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(-0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080169454) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(-2.9598575) q[0];
rz(0.043047992) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(2.8607686) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8107287) q[0];
sx q[0];
rz(-0.64635902) q[0];
sx q[0];
rz(-2.3848563) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7634723) q[2];
sx q[2];
rz(-1.780605) q[2];
sx q[2];
rz(3.0378621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3481969) q[1];
sx q[1];
rz(-1.6284202) q[1];
sx q[1];
rz(0.036199526) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9032352) q[3];
sx q[3];
rz(-1.7566163) q[3];
sx q[3];
rz(2.9591054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3165555) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(-2.5184856) q[2];
rz(2.1394219) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(-0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(2.2676246) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-0.95332425) q[2];
sx q[2];
rz(-0.8669903) q[2];
sx q[2];
rz(0.16624761) q[2];
rz(1.9813886) q[3];
sx q[3];
rz(-1.2413597) q[3];
sx q[3];
rz(1.6551457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
