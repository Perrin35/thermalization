OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(-0.46407035) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(-0.067117604) q[1];
sx q[1];
rz(1.2844515) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0425134) q[0];
sx q[0];
rz(-0.49603841) q[0];
sx q[0];
rz(-0.92143671) q[0];
rz(-pi) q[1];
rz(-1.60381) q[2];
sx q[2];
rz(-2.0602977) q[2];
sx q[2];
rz(-2.2061493) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9570934) q[1];
sx q[1];
rz(-0.59809369) q[1];
sx q[1];
rz(2.4615272) q[1];
rz(-pi) q[2];
rz(-1.5873317) q[3];
sx q[3];
rz(-1.4201122) q[3];
sx q[3];
rz(2.2076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7444732) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(2.822067) q[2];
rz(-0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-0.67392504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4085061) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(0.52655667) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(0.79663509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986886) q[0];
sx q[0];
rz(-1.8968549) q[0];
sx q[0];
rz(-2.23404) q[0];
rz(-2.8393306) q[2];
sx q[2];
rz(-0.59603359) q[2];
sx q[2];
rz(0.1421393) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14256515) q[1];
sx q[1];
rz(-1.5579281) q[1];
sx q[1];
rz(0.60728118) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2315352) q[3];
sx q[3];
rz(-0.72477341) q[3];
sx q[3];
rz(-1.4409325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9600296) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(0.78655085) q[2];
rz(-2.6484047) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86984533) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.7012117) q[0];
rz(-2.4213743) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(-2.4386141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7754606) q[0];
sx q[0];
rz(-1.4538987) q[0];
sx q[0];
rz(-1.2890105) q[0];
x q[1];
rz(-0.59962745) q[2];
sx q[2];
rz(-2.5587974) q[2];
sx q[2];
rz(-1.8011013) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2762201) q[1];
sx q[1];
rz(-3.0225388) q[1];
sx q[1];
rz(-0.40465506) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3316783) q[3];
sx q[3];
rz(-2.5044887) q[3];
sx q[3];
rz(-1.0106196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26677033) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(-0.91397816) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.5660969) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(1.874118) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(-0.56328303) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7960498) q[0];
sx q[0];
rz(-2.4420218) q[0];
sx q[0];
rz(-2.0027341) q[0];
x q[1];
rz(-1.2345384) q[2];
sx q[2];
rz(-2.2328937) q[2];
sx q[2];
rz(-1.7788356) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56983419) q[1];
sx q[1];
rz(-0.27184871) q[1];
sx q[1];
rz(-3.1175201) q[1];
rz(-0.33470811) q[3];
sx q[3];
rz(-0.74581205) q[3];
sx q[3];
rz(-1.0583744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48878601) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-2.9768067) q[2];
rz(0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.16211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441372) q[0];
sx q[0];
rz(-1.9636969) q[0];
sx q[0];
rz(-0.59209728) q[0];
rz(-2.2010872) q[2];
sx q[2];
rz(-2.5224707) q[2];
sx q[2];
rz(1.3319912) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1197549) q[1];
sx q[1];
rz(-2.2720085) q[1];
sx q[1];
rz(2.5874058) q[1];
rz(-1.5858298) q[3];
sx q[3];
rz(-0.84926499) q[3];
sx q[3];
rz(-2.815747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.44624415) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(1.2472786) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(0.22918992) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(-1.3060588) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7040946) q[0];
sx q[0];
rz(-1.9537582) q[0];
sx q[0];
rz(-1.0136481) q[0];
rz(-pi) q[1];
rz(-1.3681709) q[2];
sx q[2];
rz(-2.0261923) q[2];
sx q[2];
rz(0.95603285) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29014153) q[1];
sx q[1];
rz(-2.4212004) q[1];
sx q[1];
rz(-0.011522567) q[1];
rz(-1.7197051) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7541472) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(2.5409017) q[2];
rz(-2.1389652) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(2.7887662) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3951185) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(1.0466928) q[0];
rz(1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(2.7244862) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313275) q[0];
sx q[0];
rz(-1.7908887) q[0];
sx q[0];
rz(-0.025028153) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8646556) q[2];
sx q[2];
rz(-1.8655348) q[2];
sx q[2];
rz(2.7851832) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6228094) q[1];
sx q[1];
rz(-2.2768524) q[1];
sx q[1];
rz(0.6219567) q[1];
rz(-2.7583371) q[3];
sx q[3];
rz(-0.68985046) q[3];
sx q[3];
rz(-2.3243429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0059011857) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(1.3674412) q[2];
rz(0.55316365) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(-1.9518071) q[0];
rz(1.4272383) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(-3.0292125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0880786) q[0];
sx q[0];
rz(-1.3980165) q[0];
sx q[0];
rz(-1.0900351) q[0];
rz(-pi) q[1];
rz(2.8078812) q[2];
sx q[2];
rz(-1.6534272) q[2];
sx q[2];
rz(2.4090648) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8466859) q[1];
sx q[1];
rz(-1.2290188) q[1];
sx q[1];
rz(0.81975598) q[1];
rz(-pi) q[2];
rz(1.3120679) q[3];
sx q[3];
rz(-1.6885763) q[3];
sx q[3];
rz(3.0974914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0743951) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.3486264) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(0.034974139) q[0];
rz(-0.84683013) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(2.2299178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.864894) q[0];
sx q[0];
rz(-1.8407341) q[0];
sx q[0];
rz(-0.2322659) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5981204) q[2];
sx q[2];
rz(-1.4620355) q[2];
sx q[2];
rz(-0.47765884) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.94757838) q[1];
sx q[1];
rz(-1.7742426) q[1];
sx q[1];
rz(2.6737763) q[1];
x q[2];
rz(0.53578844) q[3];
sx q[3];
rz(-1.1402604) q[3];
sx q[3];
rz(-2.146194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.70790616) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(0.76672673) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(-0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837759) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(-2.7695079) q[0];
rz(-2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.7262329) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2863662) q[0];
sx q[0];
rz(-0.33272538) q[0];
sx q[0];
rz(-0.65441982) q[0];
x q[1];
rz(-2.4731589) q[2];
sx q[2];
rz(-2.4032776) q[2];
sx q[2];
rz(0.43503209) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33503867) q[1];
sx q[1];
rz(-1.6400669) q[1];
sx q[1];
rz(0.075242234) q[1];
rz(-pi) q[2];
rz(-1.0269208) q[3];
sx q[3];
rz(-1.2380935) q[3];
sx q[3];
rz(-1.1247016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8156478) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(0.20467219) q[2];
rz(1.7278016) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(-2.0457101) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6407912) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(1.3265058) q[2];
sx q[2];
rz(-2.6752224) q[2];
sx q[2];
rz(0.35717076) q[2];
rz(1.6987883) q[3];
sx q[3];
rz(-0.57590579) q[3];
sx q[3];
rz(-0.87344575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
