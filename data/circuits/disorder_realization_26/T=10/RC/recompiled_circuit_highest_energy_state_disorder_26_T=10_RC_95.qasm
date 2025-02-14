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
rz(2.8925722) q[0];
sx q[0];
rz(-2.1828716) q[0];
sx q[0];
rz(-2.9160461) q[0];
rz(-1.072999) q[1];
sx q[1];
rz(3.5993242) q[1];
sx q[1];
rz(9.2158894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5208369) q[0];
sx q[0];
rz(-1.136047) q[0];
sx q[0];
rz(1.8692383) q[0];
rz(-pi) q[1];
rz(-1.289008) q[2];
sx q[2];
rz(-0.95880858) q[2];
sx q[2];
rz(-0.91152292) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7876405) q[1];
sx q[1];
rz(-2.5293416) q[1];
sx q[1];
rz(2.5933215) q[1];
x q[2];
rz(-2.8257583) q[3];
sx q[3];
rz(-1.2441934) q[3];
sx q[3];
rz(-2.0706319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0245584) q[2];
sx q[2];
rz(-2.890675) q[2];
sx q[2];
rz(2.2158902) q[2];
rz(2.3671345) q[3];
sx q[3];
rz(-1.2239933) q[3];
sx q[3];
rz(-1.8584049) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9722612) q[0];
sx q[0];
rz(-2.4915578) q[0];
sx q[0];
rz(1.6768804) q[0];
rz(0.88893923) q[1];
sx q[1];
rz(-2.2496532) q[1];
sx q[1];
rz(-1.3533786) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.618093) q[0];
sx q[0];
rz(-0.67636469) q[0];
sx q[0];
rz(0.47250611) q[0];
rz(-pi) q[1];
rz(0.36608669) q[2];
sx q[2];
rz(-1.5913858) q[2];
sx q[2];
rz(2.1799157) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96129744) q[1];
sx q[1];
rz(-0.63034454) q[1];
sx q[1];
rz(2.784009) q[1];
x q[2];
rz(-2.3634745) q[3];
sx q[3];
rz(-1.3740842) q[3];
sx q[3];
rz(2.9784378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9951524) q[2];
sx q[2];
rz(-1.012994) q[2];
sx q[2];
rz(0.97274441) q[2];
rz(-1.3257596) q[3];
sx q[3];
rz(-1.1498068) q[3];
sx q[3];
rz(-2.5241234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72508088) q[0];
sx q[0];
rz(-1.4840115) q[0];
sx q[0];
rz(3.1120279) q[0];
rz(-1.0211396) q[1];
sx q[1];
rz(-2.857326) q[1];
sx q[1];
rz(2.8254642) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97702867) q[0];
sx q[0];
rz(-0.029677756) q[0];
sx q[0];
rz(2.3182431) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7720376) q[2];
sx q[2];
rz(-1.4788718) q[2];
sx q[2];
rz(-0.52009174) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7252032) q[1];
sx q[1];
rz(-1.3969743) q[1];
sx q[1];
rz(-1.9085788) q[1];
rz(-1.49316) q[3];
sx q[3];
rz(-0.85143733) q[3];
sx q[3];
rz(0.84859728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0315087) q[2];
sx q[2];
rz(-1.2981334) q[2];
sx q[2];
rz(-2.617344) q[2];
rz(0.60018572) q[3];
sx q[3];
rz(-0.46349183) q[3];
sx q[3];
rz(-2.0704827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6780739) q[0];
sx q[0];
rz(-1.2200032) q[0];
sx q[0];
rz(2.779261) q[0];
rz(-0.70829567) q[1];
sx q[1];
rz(-0.34168044) q[1];
sx q[1];
rz(1.1452311) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6041421) q[0];
sx q[0];
rz(-1.6509127) q[0];
sx q[0];
rz(-1.6782922) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6537204) q[2];
sx q[2];
rz(-2.2053501) q[2];
sx q[2];
rz(1.799859) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2140183) q[1];
sx q[1];
rz(-1.8248789) q[1];
sx q[1];
rz(2.2831282) q[1];
x q[2];
rz(-1.5208552) q[3];
sx q[3];
rz(-1.7347214) q[3];
sx q[3];
rz(2.0633179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5273253) q[2];
sx q[2];
rz(-2.3136487) q[2];
sx q[2];
rz(0.0086199363) q[2];
rz(0.65420592) q[3];
sx q[3];
rz(-2.4253186) q[3];
sx q[3];
rz(0.27230984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24097405) q[0];
sx q[0];
rz(-0.32855496) q[0];
sx q[0];
rz(0.7884489) q[0];
rz(2.12766) q[1];
sx q[1];
rz(-1.8221816) q[1];
sx q[1];
rz(0.088277146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7006314) q[0];
sx q[0];
rz(-0.92871237) q[0];
sx q[0];
rz(-0.93812801) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6730673) q[2];
sx q[2];
rz(-2.132093) q[2];
sx q[2];
rz(-3.1401538) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4406529) q[1];
sx q[1];
rz(-1.0723317) q[1];
sx q[1];
rz(2.8200214) q[1];
x q[2];
rz(-1.621219) q[3];
sx q[3];
rz(-2.6459624) q[3];
sx q[3];
rz(2.4448066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91904116) q[2];
sx q[2];
rz(-1.4107979) q[2];
sx q[2];
rz(0.50773531) q[2];
rz(3.0159085) q[3];
sx q[3];
rz(-1.9583743) q[3];
sx q[3];
rz(2.3465033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.3313726) q[0];
sx q[0];
rz(-0.75646821) q[0];
sx q[0];
rz(-0.08547011) q[0];
rz(-2.5147009) q[1];
sx q[1];
rz(-1.1566999) q[1];
sx q[1];
rz(-1.8524106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3835433) q[0];
sx q[0];
rz(-1.6746759) q[0];
sx q[0];
rz(-1.3707177) q[0];
rz(-2.9741227) q[2];
sx q[2];
rz(-2.390927) q[2];
sx q[2];
rz(1.0274061) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7886145) q[1];
sx q[1];
rz(-1.7403688) q[1];
sx q[1];
rz(1.3296849) q[1];
x q[2];
rz(-1.3335532) q[3];
sx q[3];
rz(-1.5565775) q[3];
sx q[3];
rz(-2.1923101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.36279303) q[2];
sx q[2];
rz(-2.1976566) q[2];
sx q[2];
rz(-1.8865406) q[2];
rz(3.0826027) q[3];
sx q[3];
rz(-1.9374282) q[3];
sx q[3];
rz(0.98439938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.4208218) q[0];
sx q[0];
rz(-1.868792) q[0];
sx q[0];
rz(-0.76501784) q[0];
rz(1.4014686) q[1];
sx q[1];
rz(-0.88686371) q[1];
sx q[1];
rz(2.9041362) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260587) q[0];
sx q[0];
rz(-2.1977673) q[0];
sx q[0];
rz(-3.0978908) q[0];
rz(-pi) q[1];
rz(3.0661707) q[2];
sx q[2];
rz(-0.89176501) q[2];
sx q[2];
rz(-1.1231658) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85126696) q[1];
sx q[1];
rz(-1.2801761) q[1];
sx q[1];
rz(-1.5253011) q[1];
rz(-pi) q[2];
rz(-1.6066437) q[3];
sx q[3];
rz(-2.4447943) q[3];
sx q[3];
rz(-3.0655332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.456363) q[2];
sx q[2];
rz(-2.3380029) q[2];
sx q[2];
rz(-0.85344273) q[2];
rz(1.338909) q[3];
sx q[3];
rz(-1.5693376) q[3];
sx q[3];
rz(-2.9760823) q[3];
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
rz(2.6859739) q[0];
sx q[0];
rz(-1.2116665) q[0];
sx q[0];
rz(-0.18534216) q[0];
rz(-2.7475884) q[1];
sx q[1];
rz(-1.3961671) q[1];
sx q[1];
rz(-1.0967163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1857796) q[0];
sx q[0];
rz(-2.5086864) q[0];
sx q[0];
rz(-0.88297412) q[0];
rz(-pi) q[1];
rz(1.4273066) q[2];
sx q[2];
rz(-2.4713785) q[2];
sx q[2];
rz(-2.0795979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85307676) q[1];
sx q[1];
rz(-2.4429053) q[1];
sx q[1];
rz(-2.6050911) q[1];
x q[2];
rz(0.53727178) q[3];
sx q[3];
rz(-2.4747765) q[3];
sx q[3];
rz(1.3992573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77157053) q[2];
sx q[2];
rz(-3.0574419) q[2];
sx q[2];
rz(-0.30878511) q[2];
rz(1.8874946) q[3];
sx q[3];
rz(-1.8697238) q[3];
sx q[3];
rz(2.0103644) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18987385) q[0];
sx q[0];
rz(-1.254344) q[0];
sx q[0];
rz(0.21981123) q[0];
rz(1.9728569) q[1];
sx q[1];
rz(-1.7857779) q[1];
sx q[1];
rz(-1.2723602) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4502122) q[0];
sx q[0];
rz(-2.5369329) q[0];
sx q[0];
rz(1.613841) q[0];
rz(-pi) q[1];
rz(2.9909439) q[2];
sx q[2];
rz(-0.51567557) q[2];
sx q[2];
rz(-0.0055238481) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2337964) q[1];
sx q[1];
rz(-1.7410454) q[1];
sx q[1];
rz(-2.7336804) q[1];
rz(0.50140372) q[3];
sx q[3];
rz(-2.025017) q[3];
sx q[3];
rz(2.8387031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9610338) q[2];
sx q[2];
rz(-2.2982633) q[2];
sx q[2];
rz(-2.6832704) q[2];
rz(-2.787163) q[3];
sx q[3];
rz(-1.7731881) q[3];
sx q[3];
rz(1.1752769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71427041) q[0];
sx q[0];
rz(-1.9388119) q[0];
sx q[0];
rz(0.74139968) q[0];
rz(-1.4330014) q[1];
sx q[1];
rz(-2.7129136) q[1];
sx q[1];
rz(-3.0523849) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36486711) q[0];
sx q[0];
rz(-1.3898384) q[0];
sx q[0];
rz(-1.0904161) q[0];
rz(-pi) q[1];
rz(-1.9544425) q[2];
sx q[2];
rz(-2.128278) q[2];
sx q[2];
rz(-1.5889744) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9392689) q[1];
sx q[1];
rz(-1.2319733) q[1];
sx q[1];
rz(0.044087709) q[1];
rz(-pi) q[2];
rz(2.6939635) q[3];
sx q[3];
rz(-2.3530686) q[3];
sx q[3];
rz(-2.727689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0006813) q[2];
sx q[2];
rz(-3.0912283) q[2];
sx q[2];
rz(0.64765206) q[2];
rz(2.7378313) q[3];
sx q[3];
rz(-1.8881366) q[3];
sx q[3];
rz(-3.0103179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.1972926) q[0];
sx q[0];
rz(-2.1680752) q[0];
sx q[0];
rz(-2.0808921) q[0];
rz(0.79509673) q[1];
sx q[1];
rz(-0.92082321) q[1];
sx q[1];
rz(-0.77650741) q[1];
rz(2.2262103) q[2];
sx q[2];
rz(-1.4861098) q[2];
sx q[2];
rz(-2.8783023) q[2];
rz(3.0632302) q[3];
sx q[3];
rz(-1.4808803) q[3];
sx q[3];
rz(-0.63054569) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
