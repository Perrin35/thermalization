OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(-1.6337448) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(3.9197796) q[1];
sx q[1];
rz(9.9240886) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0190174) q[0];
sx q[0];
rz(-2.8589773) q[0];
sx q[0];
rz(-0.22960381) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7668031) q[2];
sx q[2];
rz(-0.73050806) q[2];
sx q[2];
rz(1.0175878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62568134) q[1];
sx q[1];
rz(-1.7209007) q[1];
sx q[1];
rz(2.4019269) q[1];
rz(1.4673759) q[3];
sx q[3];
rz(-1.5973583) q[3];
sx q[3];
rz(2.6288222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(2.9075918) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4734128) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(-2.1372674) q[0];
rz(-1.6218119) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(-1.0027592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5402055) q[0];
sx q[0];
rz(-1.1367072) q[0];
sx q[0];
rz(-2.4853404) q[0];
rz(-pi) q[1];
rz(-0.82680306) q[2];
sx q[2];
rz(-2.8892592) q[2];
sx q[2];
rz(-1.2834872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7748389) q[1];
sx q[1];
rz(-2.5336821) q[1];
sx q[1];
rz(1.2971836) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59870436) q[3];
sx q[3];
rz(-1.8667392) q[3];
sx q[3];
rz(0.89950365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(-0.15094748) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(-2.7422854) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-2.2580106) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8944091) q[0];
sx q[0];
rz(-2.7054225) q[0];
sx q[0];
rz(2.7437074) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1951139) q[2];
sx q[2];
rz(-0.4779856) q[2];
sx q[2];
rz(-1.1849272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98595847) q[1];
sx q[1];
rz(-1.2859207) q[1];
sx q[1];
rz(3.1180711) q[1];
rz(1.4530573) q[3];
sx q[3];
rz(-0.91500926) q[3];
sx q[3];
rz(1.9829139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3016004) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(-2.1515576) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(-2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7097968) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(2.0344095) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(-2.591419) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.195939) q[0];
sx q[0];
rz(-1.6473856) q[0];
sx q[0];
rz(0.73777317) q[0];
rz(-pi) q[1];
rz(2.8407211) q[2];
sx q[2];
rz(-2.8340325) q[2];
sx q[2];
rz(-3.0645264) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6039227) q[1];
sx q[1];
rz(-2.1426755) q[1];
sx q[1];
rz(-3.0043383) q[1];
rz(1.5189734) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(-2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6198373) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6977285) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(0.91526389) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6633776) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(-2.1992654) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5259597) q[2];
sx q[2];
rz(-1.0460639) q[2];
sx q[2];
rz(0.63477883) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.535714) q[1];
sx q[1];
rz(-1.5389581) q[1];
sx q[1];
rz(-1.4153626) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4004585) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(2.4853064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(-2.6780224) q[2];
rz(-2.5772337) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(-0.92818964) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(1.0790496) q[0];
rz(-0.59533978) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(0.39658305) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88543372) q[0];
sx q[0];
rz(-2.0246756) q[0];
sx q[0];
rz(-2.8594349) q[0];
rz(-pi) q[1];
rz(-1.2601389) q[2];
sx q[2];
rz(-1.1144181) q[2];
sx q[2];
rz(0.059547101) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5844333) q[1];
sx q[1];
rz(-1.6646619) q[1];
sx q[1];
rz(0.90969109) q[1];
rz(-pi) q[2];
rz(-0.58432213) q[3];
sx q[3];
rz(-1.366368) q[3];
sx q[3];
rz(-2.1401329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1534999) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(-1.5863824) q[2];
rz(-1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-2.3144408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.39847386) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(-1.8709042) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(-2.0152337) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2369909) q[0];
sx q[0];
rz(-0.31123529) q[0];
sx q[0];
rz(-0.29445946) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0790714) q[2];
sx q[2];
rz(-1.695343) q[2];
sx q[2];
rz(1.6342271) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78325242) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(-1.7794442) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.095110006) q[3];
sx q[3];
rz(-1.2246119) q[3];
sx q[3];
rz(-0.1971052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9324947) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(-3.0287108) q[0];
rz(2.1408634) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(3.1088366) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.207064) q[0];
sx q[0];
rz(-1.3158187) q[0];
sx q[0];
rz(-1.2364788) q[0];
x q[1];
rz(1.007349) q[2];
sx q[2];
rz(-1.9512366) q[2];
sx q[2];
rz(3.0055339) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9084839) q[1];
sx q[1];
rz(-1.9875803) q[1];
sx q[1];
rz(-3.0770739) q[1];
rz(-pi) q[2];
rz(1.7912346) q[3];
sx q[3];
rz(-2.7535097) q[3];
sx q[3];
rz(-2.2409852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(0.62210554) q[2];
rz(-1.9744251) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.099982925) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(1.6149678) q[0];
rz(-2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(0.48535767) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55886666) q[0];
sx q[0];
rz(-0.14478806) q[0];
sx q[0];
rz(2.7607714) q[0];
rz(-pi) q[1];
rz(1.7301241) q[2];
sx q[2];
rz(-2.1630641) q[2];
sx q[2];
rz(1.2194022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7945054) q[1];
sx q[1];
rz(-2.3350888) q[1];
sx q[1];
rz(2.0054714) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0210578) q[3];
sx q[3];
rz(-2.1053208) q[3];
sx q[3];
rz(0.083732097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1099403) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(2.9821441) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(1.2592978) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(-0.5982582) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3529417) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
x q[1];
rz(-1.7190476) q[2];
sx q[2];
rz(-1.7808) q[2];
sx q[2];
rz(-1.6809747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64893374) q[1];
sx q[1];
rz(-1.4760541) q[1];
sx q[1];
rz(-1.2050864) q[1];
rz(-pi) q[2];
rz(0.80501276) q[3];
sx q[3];
rz(-0.22634889) q[3];
sx q[3];
rz(-1.5387907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.082211994) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(0.65199488) q[2];
rz(2.5478798) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(2.0098856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3020637) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(-1.9051753) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(2.9700206) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(-0.89647722) q[3];
sx q[3];
rz(-1.258068) q[3];
sx q[3];
rz(0.51545959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
