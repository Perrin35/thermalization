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
rz(0.23586805) q[0];
sx q[0];
rz(-1.1455102) q[0];
sx q[0];
rz(-3.0930162) q[0];
rz(1.5161169) q[1];
sx q[1];
rz(-2.0487509) q[1];
sx q[1];
rz(-0.71576524) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5465821) q[0];
sx q[0];
rz(-1.6002185) q[0];
sx q[0];
rz(1.3999697) q[0];
x q[1];
rz(-0.20697941) q[2];
sx q[2];
rz(-1.845115) q[2];
sx q[2];
rz(1.5766608) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6977344) q[1];
sx q[1];
rz(-0.83699534) q[1];
sx q[1];
rz(-0.55056449) q[1];
rz(-pi) q[2];
rz(-0.64438246) q[3];
sx q[3];
rz(-2.5231859) q[3];
sx q[3];
rz(1.9559086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6526661) q[2];
sx q[2];
rz(-0.55997866) q[2];
sx q[2];
rz(2.0860591) q[2];
rz(-0.40600285) q[3];
sx q[3];
rz(-1.8279653) q[3];
sx q[3];
rz(0.82160151) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15776289) q[0];
sx q[0];
rz(-1.0394179) q[0];
sx q[0];
rz(-2.5854172) q[0];
rz(1.3293386) q[1];
sx q[1];
rz(-2.5228597) q[1];
sx q[1];
rz(-3.0583196) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0389685) q[0];
sx q[0];
rz(-0.73154615) q[0];
sx q[0];
rz(0.49517531) q[0];
x q[1];
rz(-2.1883972) q[2];
sx q[2];
rz(-2.0591271) q[2];
sx q[2];
rz(-0.40541475) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.00035380414) q[1];
sx q[1];
rz(-1.7906931) q[1];
sx q[1];
rz(-2.6882437) q[1];
rz(-pi) q[2];
rz(-1.0038297) q[3];
sx q[3];
rz(-1.1440082) q[3];
sx q[3];
rz(0.7184283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.1978153) q[2];
sx q[2];
rz(-1.3894206) q[2];
sx q[2];
rz(0.49276349) q[2];
rz(2.1150151) q[3];
sx q[3];
rz(-2.8387098) q[3];
sx q[3];
rz(-1.4125642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3781085) q[0];
sx q[0];
rz(-0.56270993) q[0];
sx q[0];
rz(2.7028131) q[0];
rz(1.6890866) q[1];
sx q[1];
rz(-0.32183281) q[1];
sx q[1];
rz(0.39295331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38848061) q[0];
sx q[0];
rz(-1.5428758) q[0];
sx q[0];
rz(2.0620538) q[0];
rz(2.0991481) q[2];
sx q[2];
rz(-1.5425041) q[2];
sx q[2];
rz(-1.229778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4247228) q[1];
sx q[1];
rz(-0.42415127) q[1];
sx q[1];
rz(-1.7440026) q[1];
rz(-pi) q[2];
rz(-3.0370029) q[3];
sx q[3];
rz(-2.4896693) q[3];
sx q[3];
rz(-0.014217941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4517333) q[2];
sx q[2];
rz(-2.2971575) q[2];
sx q[2];
rz(2.8680657) q[2];
rz(-2.5290329) q[3];
sx q[3];
rz(-1.4140244) q[3];
sx q[3];
rz(-0.85649049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81371152) q[0];
sx q[0];
rz(-2.5484945) q[0];
sx q[0];
rz(2.3408422) q[0];
rz(0.71209359) q[1];
sx q[1];
rz(-1.1199896) q[1];
sx q[1];
rz(2.6584279) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8022435) q[0];
sx q[0];
rz(-2.3988945) q[0];
sx q[0];
rz(-0.00089641103) q[0];
rz(-2.0676548) q[2];
sx q[2];
rz(-2.0570847) q[2];
sx q[2];
rz(-0.72975791) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3222683) q[1];
sx q[1];
rz(-1.5996278) q[1];
sx q[1];
rz(-2.4954456) q[1];
rz(0.40164803) q[3];
sx q[3];
rz(-0.73026087) q[3];
sx q[3];
rz(0.83056565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.157865) q[2];
sx q[2];
rz(-2.4309776) q[2];
sx q[2];
rz(1.857081) q[2];
rz(-0.48786783) q[3];
sx q[3];
rz(-1.7043461) q[3];
sx q[3];
rz(-1.6871066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81216413) q[0];
sx q[0];
rz(-0.85947961) q[0];
sx q[0];
rz(-0.46464768) q[0];
rz(-2.1931785) q[1];
sx q[1];
rz(-0.83810884) q[1];
sx q[1];
rz(-1.4533739) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3188948) q[0];
sx q[0];
rz(-1.5359211) q[0];
sx q[0];
rz(-0.025392763) q[0];
x q[1];
rz(-2.862013) q[2];
sx q[2];
rz(-0.76020066) q[2];
sx q[2];
rz(0.81499824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0665474) q[1];
sx q[1];
rz(-2.5270224) q[1];
sx q[1];
rz(3.0250508) q[1];
x q[2];
rz(-0.57027633) q[3];
sx q[3];
rz(-0.81109427) q[3];
sx q[3];
rz(0.18660422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98443085) q[2];
sx q[2];
rz(-1.0474019) q[2];
sx q[2];
rz(2.7471527) q[2];
rz(1.1557584) q[3];
sx q[3];
rz(-2.211536) q[3];
sx q[3];
rz(-1.3084779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406469) q[0];
sx q[0];
rz(-2.7257958) q[0];
sx q[0];
rz(-2.0796602) q[0];
rz(0.95147079) q[1];
sx q[1];
rz(-0.87363344) q[1];
sx q[1];
rz(1.5207312) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5439107) q[0];
sx q[0];
rz(-1.8476227) q[0];
sx q[0];
rz(-0.055276543) q[0];
x q[1];
rz(-0.64487793) q[2];
sx q[2];
rz(-1.0928003) q[2];
sx q[2];
rz(2.1410112) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1870059) q[1];
sx q[1];
rz(-2.329064) q[1];
sx q[1];
rz(1.7962667) q[1];
rz(2.6470408) q[3];
sx q[3];
rz(-0.44122094) q[3];
sx q[3];
rz(2.0989204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.025853) q[2];
sx q[2];
rz(-1.9151177) q[2];
sx q[2];
rz(-0.24570492) q[2];
rz(1.687441) q[3];
sx q[3];
rz(-1.629963) q[3];
sx q[3];
rz(-0.83034849) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65569735) q[0];
sx q[0];
rz(-2.4035154) q[0];
sx q[0];
rz(0.6947211) q[0];
rz(3.0409536) q[1];
sx q[1];
rz(-2.1058319) q[1];
sx q[1];
rz(-0.86520854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9964906) q[0];
sx q[0];
rz(-1.9887513) q[0];
sx q[0];
rz(-1.6161902) q[0];
rz(-pi) q[1];
rz(2.2312282) q[2];
sx q[2];
rz(-1.5558793) q[2];
sx q[2];
rz(0.93176022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2336967) q[1];
sx q[1];
rz(-1.421126) q[1];
sx q[1];
rz(2.187633) q[1];
x q[2];
rz(2.4669924) q[3];
sx q[3];
rz(-1.1230506) q[3];
sx q[3];
rz(0.49481667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9575017) q[2];
sx q[2];
rz(-2.3881193) q[2];
sx q[2];
rz(-0.91020477) q[2];
rz(-1.4522067) q[3];
sx q[3];
rz(-1.9096749) q[3];
sx q[3];
rz(0.46428251) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8978187) q[0];
sx q[0];
rz(-2.0938566) q[0];
sx q[0];
rz(-2.5183103) q[0];
rz(-2.9858164) q[1];
sx q[1];
rz(-0.77729762) q[1];
sx q[1];
rz(2.8782841) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52180852) q[0];
sx q[0];
rz(-1.5814879) q[0];
sx q[0];
rz(1.7858206) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66923334) q[2];
sx q[2];
rz(-0.097429052) q[2];
sx q[2];
rz(2.214746) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1340848) q[1];
sx q[1];
rz(-0.99564394) q[1];
sx q[1];
rz(-1.47912) q[1];
rz(-2.1144763) q[3];
sx q[3];
rz(-0.86302084) q[3];
sx q[3];
rz(-2.2563427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89173633) q[2];
sx q[2];
rz(-2.8920434) q[2];
sx q[2];
rz(0.88805324) q[2];
rz(-0.92787162) q[3];
sx q[3];
rz(-1.1279305) q[3];
sx q[3];
rz(-2.1685062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.4230098) q[0];
sx q[0];
rz(-2.6319478) q[0];
sx q[0];
rz(2.6408559) q[0];
rz(-2.0210733) q[1];
sx q[1];
rz(-0.78103939) q[1];
sx q[1];
rz(-0.80148554) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.00068110739) q[0];
sx q[0];
rz(-1.4998224) q[0];
sx q[0];
rz(1.5727726) q[0];
x q[1];
rz(-0.74023418) q[2];
sx q[2];
rz(-2.4317305) q[2];
sx q[2];
rz(-0.20698337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3858955) q[1];
sx q[1];
rz(-0.79078005) q[1];
sx q[1];
rz(-2.5189931) q[1];
x q[2];
rz(-0.19480199) q[3];
sx q[3];
rz(-1.8240989) q[3];
sx q[3];
rz(1.4089438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72006172) q[2];
sx q[2];
rz(-1.3214107) q[2];
sx q[2];
rz(0.23019543) q[2];
rz(3.1318393) q[3];
sx q[3];
rz(-1.6246656) q[3];
sx q[3];
rz(2.9856288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941876) q[0];
sx q[0];
rz(-1.7740086) q[0];
sx q[0];
rz(2.3719924) q[0];
rz(-0.96254483) q[1];
sx q[1];
rz(-1.3178408) q[1];
sx q[1];
rz(2.1544971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55961032) q[0];
sx q[0];
rz(-1.6405237) q[0];
sx q[0];
rz(-2.0980673) q[0];
x q[1];
rz(-1.3279506) q[2];
sx q[2];
rz(-1.7195903) q[2];
sx q[2];
rz(0.75243581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35008123) q[1];
sx q[1];
rz(-1.529585) q[1];
sx q[1];
rz(-3.1084395) q[1];
rz(-pi) q[2];
rz(2.6392699) q[3];
sx q[3];
rz(-2.2993616) q[3];
sx q[3];
rz(0.80020927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7139682) q[2];
sx q[2];
rz(-2.2043113) q[2];
sx q[2];
rz(3.0734708) q[2];
rz(-0.44975975) q[3];
sx q[3];
rz(-2.7263548) q[3];
sx q[3];
rz(1.4382039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.149067) q[0];
sx q[0];
rz(-1.5328007) q[0];
sx q[0];
rz(-1.3971064) q[0];
rz(2.8515011) q[1];
sx q[1];
rz(-2.7128704) q[1];
sx q[1];
rz(-2.6147978) q[1];
rz(-2.5505603) q[2];
sx q[2];
rz(-2.5606511) q[2];
sx q[2];
rz(-2.64369) q[2];
rz(0.3984821) q[3];
sx q[3];
rz(-1.5881817) q[3];
sx q[3];
rz(2.6276799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
