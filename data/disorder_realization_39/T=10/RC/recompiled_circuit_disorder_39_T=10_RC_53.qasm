OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(-2.0733641) q[0];
sx q[0];
rz(0.46407035) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0811631) q[0];
sx q[0];
rz(-1.2788749) q[0];
sx q[0];
rz(-1.1638327) q[0];
rz(-pi) q[1];
rz(0.48972763) q[2];
sx q[2];
rz(-1.5999319) q[2];
sx q[2];
rz(2.4907128) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7317036) q[1];
sx q[1];
rz(-1.1176425) q[1];
sx q[1];
rz(1.1660006) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5873317) q[3];
sx q[3];
rz(-1.7214805) q[3];
sx q[3];
rz(-2.2076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7444732) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(2.822067) q[2];
rz(-2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330866) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(0.52655667) q[0];
rz(-2.5780442) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(-2.3449576) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8684727) q[0];
sx q[0];
rz(-2.1935049) q[0];
sx q[0];
rz(0.40533439) q[0];
x q[1];
rz(1.37155) q[2];
sx q[2];
rz(-1.0052048) q[2];
sx q[2];
rz(2.9233962) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9990275) q[1];
sx q[1];
rz(-1.5579281) q[1];
sx q[1];
rz(-0.60728118) q[1];
rz(-pi) q[2];
rz(-2.2315352) q[3];
sx q[3];
rz(-0.72477341) q[3];
sx q[3];
rz(1.4409325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18156302) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(-0.44737679) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(1.440381) q[0];
rz(0.72021833) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(-2.4386141) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1709135) q[0];
sx q[0];
rz(-1.850607) q[0];
sx q[0];
rz(-3.0199416) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5419652) q[2];
sx q[2];
rz(-0.58279524) q[2];
sx q[2];
rz(-1.3404913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8382593) q[1];
sx q[1];
rz(-1.6175744) q[1];
sx q[1];
rz(0.10951885) q[1];
rz(0.94743518) q[3];
sx q[3];
rz(-1.4294335) q[3];
sx q[3];
rz(0.75368222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(0.91153574) q[2];
rz(-0.91397816) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(-2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754958) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(2.3655868) q[0];
rz(-1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-2.5783096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7960498) q[0];
sx q[0];
rz(-2.4420218) q[0];
sx q[0];
rz(1.1388586) q[0];
rz(-pi) q[1];
rz(-1.9070542) q[2];
sx q[2];
rz(-2.2328937) q[2];
sx q[2];
rz(1.7788356) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56983419) q[1];
sx q[1];
rz(-2.8697439) q[1];
sx q[1];
rz(0.024072577) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8654278) q[3];
sx q[3];
rz(-0.8751103) q[3];
sx q[3];
rz(0.61616117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(-2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.9794827) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441372) q[0];
sx q[0];
rz(-1.1778957) q[0];
sx q[0];
rz(-2.5494954) q[0];
rz(1.048462) q[2];
sx q[2];
rz(-1.9198717) q[2];
sx q[2];
rz(-0.77490865) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4967242) q[1];
sx q[1];
rz(-0.86360303) q[1];
sx q[1];
rz(2.1281388) q[1];
rz(2.4200053) q[3];
sx q[3];
rz(-1.5595094) q[3];
sx q[3];
rz(1.23502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.44624415) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(-1.8943141) q[2];
rz(0.013109664) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(-0.22918992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86876774) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(-0.75063467) q[0];
rz(2.0320832) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40697843) q[0];
sx q[0];
rz(-2.4771871) q[0];
sx q[0];
rz(-2.2218496) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3681709) q[2];
sx q[2];
rz(-2.0261923) q[2];
sx q[2];
rz(2.1855598) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8667824) q[1];
sx q[1];
rz(-0.85046235) q[1];
sx q[1];
rz(-1.5606828) q[1];
rz(0.2209729) q[3];
sx q[3];
rz(-2.5413725) q[3];
sx q[3];
rz(-0.61629399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3874454) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(-2.5409017) q[2];
rz(-2.1389652) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(-0.35282648) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464741) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(1.5294317) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(-0.41710645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1244303) q[0];
sx q[0];
rz(-2.9201047) q[0];
sx q[0];
rz(1.4593967) q[0];
x q[1];
rz(2.8344526) q[2];
sx q[2];
rz(-1.2899613) q[2];
sx q[2];
rz(-1.3020696) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.48737803) q[1];
sx q[1];
rz(-2.0301135) q[1];
sx q[1];
rz(2.3801801) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2715289) q[3];
sx q[3];
rz(-0.93942681) q[3];
sx q[3];
rz(-2.805998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1356915) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(0.55316365) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.816514) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(1.9518071) q[0];
rz(-1.4272383) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(3.0292125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9427129) q[0];
sx q[0];
rz(-2.6330224) q[0];
sx q[0];
rz(-1.93165) q[0];
x q[1];
rz(-2.8078812) q[2];
sx q[2];
rz(-1.4881655) q[2];
sx q[2];
rz(-0.73252788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29490678) q[1];
sx q[1];
rz(-1.9125738) q[1];
sx q[1];
rz(2.3218367) q[1];
rz(-1.8295248) q[3];
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
rz(2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.3486264) q[2];
rz(1.2049234) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(-0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(3.1066185) q[0];
rz(2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(2.2299178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9104886) q[0];
sx q[0];
rz(-1.7945053) q[0];
sx q[0];
rz(-1.847812) q[0];
x q[1];
rz(2.5981204) q[2];
sx q[2];
rz(-1.6795571) q[2];
sx q[2];
rz(-2.6639338) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8988077) q[1];
sx q[1];
rz(-2.6344732) q[1];
sx q[1];
rz(-2.7125263) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73260143) q[3];
sx q[3];
rz(-2.4676975) q[3];
sx q[3];
rz(3.1042838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(2.8923477) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(-2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.3578167) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(-0.37208474) q[0];
rz(-0.58139873) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(-1.4153597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6045195) q[0];
sx q[0];
rz(-1.8329289) q[0];
sx q[0];
rz(-1.3634691) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62016983) q[2];
sx q[2];
rz(-2.0010741) q[2];
sx q[2];
rz(2.534453) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.806554) q[1];
sx q[1];
rz(-1.5015258) q[1];
sx q[1];
rz(3.0663504) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1595702) q[3];
sx q[3];
rz(-0.62871274) q[3];
sx q[3];
rz(0.94129896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(-2.9369205) q[2];
rz(-1.4137911) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.50080147) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(1.8150868) q[2];
sx q[2];
rz(-0.46637022) q[2];
sx q[2];
rz(-2.7844219) q[2];
rz(-3.0588991) q[3];
sx q[3];
rz(-2.1413998) q[3];
sx q[3];
rz(-1.0257046) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
