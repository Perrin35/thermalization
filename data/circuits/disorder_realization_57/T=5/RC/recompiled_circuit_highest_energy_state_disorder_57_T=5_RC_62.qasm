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
rz(-1.3402101) q[0];
sx q[0];
rz(3.4184472) q[0];
sx q[0];
rz(10.463538) q[0];
rz(2.7998595) q[1];
sx q[1];
rz(-0.83655292) q[1];
sx q[1];
rz(-0.41681448) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3562389) q[0];
sx q[0];
rz(-0.84945852) q[0];
sx q[0];
rz(1.712404) q[0];
x q[1];
rz(-2.5854163) q[2];
sx q[2];
rz(-1.3601298) q[2];
sx q[2];
rz(-2.3245101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.15910251) q[1];
sx q[1];
rz(-2.1381408) q[1];
sx q[1];
rz(0.87615396) q[1];
rz(-pi) q[2];
rz(2.9800426) q[3];
sx q[3];
rz(-1.6027502) q[3];
sx q[3];
rz(2.3475501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1699528) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(-1.8876342) q[2];
rz(1.5597255) q[3];
sx q[3];
rz(-0.73837787) q[3];
sx q[3];
rz(1.1221277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9480243) q[0];
sx q[0];
rz(-1.3688315) q[0];
sx q[0];
rz(0.76876202) q[0];
rz(1.7747152) q[1];
sx q[1];
rz(-1.3221075) q[1];
sx q[1];
rz(-2.1034525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4878142) q[0];
sx q[0];
rz(-0.013049203) q[0];
sx q[0];
rz(-2.292146) q[0];
rz(-0.72121303) q[2];
sx q[2];
rz(-1.8497457) q[2];
sx q[2];
rz(-2.9245289) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1662054) q[1];
sx q[1];
rz(-2.176099) q[1];
sx q[1];
rz(-2.0209842) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0821876) q[3];
sx q[3];
rz(-1.5540136) q[3];
sx q[3];
rz(-2.8458561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.495503) q[2];
sx q[2];
rz(-1.8547408) q[2];
sx q[2];
rz(-0.82026473) q[2];
rz(2.8881554) q[3];
sx q[3];
rz(-2.7866252) q[3];
sx q[3];
rz(-0.36111116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28061098) q[0];
sx q[0];
rz(-2.6237223) q[0];
sx q[0];
rz(-2.2542727) q[0];
rz(0.53030983) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(-1.1393772) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054087669) q[0];
sx q[0];
rz(-1.2257135) q[0];
sx q[0];
rz(-1.7214543) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33311756) q[2];
sx q[2];
rz(-2.0764167) q[2];
sx q[2];
rz(2.5207455) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4344221) q[1];
sx q[1];
rz(-1.5853197) q[1];
sx q[1];
rz(2.8959031) q[1];
rz(-pi) q[2];
x q[2];
rz(0.091413012) q[3];
sx q[3];
rz(-2.6210636) q[3];
sx q[3];
rz(-2.6176069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.016971074) q[2];
sx q[2];
rz(-1.4554224) q[2];
sx q[2];
rz(-2.7023081) q[2];
rz(-0.47075054) q[3];
sx q[3];
rz(-3.0622523) q[3];
sx q[3];
rz(-1.9973756) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73862326) q[0];
sx q[0];
rz(-0.93337494) q[0];
sx q[0];
rz(-0.11548197) q[0];
rz(-3.0237517) q[1];
sx q[1];
rz(-1.9385447) q[1];
sx q[1];
rz(-1.3062564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7457434) q[0];
sx q[0];
rz(-1.9877399) q[0];
sx q[0];
rz(2.8265619) q[0];
rz(2.6682786) q[2];
sx q[2];
rz(-0.71719786) q[2];
sx q[2];
rz(1.0953449) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0769121) q[1];
sx q[1];
rz(-2.2222328) q[1];
sx q[1];
rz(2.3802735) q[1];
rz(0.0751817) q[3];
sx q[3];
rz(-2.1422221) q[3];
sx q[3];
rz(0.21440766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8370886) q[2];
sx q[2];
rz(-1.8279165) q[2];
sx q[2];
rz(-2.9869249) q[2];
rz(1.539544) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(-0.60607564) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57347572) q[0];
sx q[0];
rz(-2.0982168) q[0];
sx q[0];
rz(2.306275) q[0];
rz(2.4256445) q[1];
sx q[1];
rz(-2.7138111) q[1];
sx q[1];
rz(0.11944019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41573856) q[0];
sx q[0];
rz(-1.6266168) q[0];
sx q[0];
rz(0.016896642) q[0];
rz(-pi) q[1];
rz(-0.30952366) q[2];
sx q[2];
rz(-0.46876568) q[2];
sx q[2];
rz(-0.13407198) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1491974) q[1];
sx q[1];
rz(-2.1572504) q[1];
sx q[1];
rz(-1.0066562) q[1];
x q[2];
rz(1.1501081) q[3];
sx q[3];
rz(-1.0771695) q[3];
sx q[3];
rz(-2.8378262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9014088) q[2];
sx q[2];
rz(-2.8808424) q[2];
sx q[2];
rz(0.060997941) q[2];
rz(0.048246233) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(-2.4129996) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7399087) q[0];
sx q[0];
rz(-2.4329199) q[0];
sx q[0];
rz(-2.0106864) q[0];
rz(-0.4785969) q[1];
sx q[1];
rz(-2.5391948) q[1];
sx q[1];
rz(-1.4071677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81086195) q[0];
sx q[0];
rz(-1.5172345) q[0];
sx q[0];
rz(-1.5723482) q[0];
rz(-pi) q[1];
rz(2.2931104) q[2];
sx q[2];
rz(-1.7847848) q[2];
sx q[2];
rz(0.011836476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.045627198) q[1];
sx q[1];
rz(-1.4573604) q[1];
sx q[1];
rz(1.3413702) q[1];
x q[2];
rz(1.3899441) q[3];
sx q[3];
rz(-1.3612124) q[3];
sx q[3];
rz(1.7169894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56167928) q[2];
sx q[2];
rz(-0.40739569) q[2];
sx q[2];
rz(1.178406) q[2];
rz(-1.0493578) q[3];
sx q[3];
rz(-0.5368084) q[3];
sx q[3];
rz(-1.2843081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.9310164) q[0];
sx q[0];
rz(-1.9518305) q[0];
sx q[0];
rz(1.806102) q[0];
rz(-1.3212851) q[1];
sx q[1];
rz(-1.0711461) q[1];
sx q[1];
rz(-0.30805045) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0548693) q[0];
sx q[0];
rz(-0.44741524) q[0];
sx q[0];
rz(0.50680508) q[0];
rz(1.9676061) q[2];
sx q[2];
rz(-1.0894629) q[2];
sx q[2];
rz(2.2725042) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2854453) q[1];
sx q[1];
rz(-0.79570192) q[1];
sx q[1];
rz(-0.85983069) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9568029) q[3];
sx q[3];
rz(-1.6148477) q[3];
sx q[3];
rz(0.06202997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7293952) q[2];
sx q[2];
rz(-0.9407548) q[2];
sx q[2];
rz(-0.13128734) q[2];
rz(2.4042551) q[3];
sx q[3];
rz(-1.3056583) q[3];
sx q[3];
rz(2.9837515) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8312296) q[0];
sx q[0];
rz(-1.0897626) q[0];
sx q[0];
rz(0.53037733) q[0];
rz(-1.7165548) q[1];
sx q[1];
rz(-1.1385463) q[1];
sx q[1];
rz(-2.4200965) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5576241) q[0];
sx q[0];
rz(-2.1707105) q[0];
sx q[0];
rz(1.0217838) q[0];
x q[1];
rz(-1.7185831) q[2];
sx q[2];
rz(-1.1061386) q[2];
sx q[2];
rz(2.8140659) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87807314) q[1];
sx q[1];
rz(-1.7612447) q[1];
sx q[1];
rz(-1.7909122) q[1];
x q[2];
rz(-0.78254487) q[3];
sx q[3];
rz(-1.6492606) q[3];
sx q[3];
rz(0.162938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7044652) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(1.6638157) q[2];
rz(-1.1635228) q[3];
sx q[3];
rz(-2.0587557) q[3];
sx q[3];
rz(1.6400953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32996938) q[0];
sx q[0];
rz(-1.4412619) q[0];
sx q[0];
rz(0.60923088) q[0];
rz(1.5902663) q[1];
sx q[1];
rz(-2.8158999) q[1];
sx q[1];
rz(-1.685166) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5824175) q[0];
sx q[0];
rz(-1.32138) q[0];
sx q[0];
rz(1.2024058) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6316192) q[2];
sx q[2];
rz(-1.0863388) q[2];
sx q[2];
rz(-1.4549507) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.85553414) q[1];
sx q[1];
rz(-1.1127591) q[1];
sx q[1];
rz(-0.18064255) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84734579) q[3];
sx q[3];
rz(-2.06978) q[3];
sx q[3];
rz(2.8049198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7234601) q[2];
sx q[2];
rz(-0.89833608) q[2];
sx q[2];
rz(2.2612803) q[2];
rz(2.9355925) q[3];
sx q[3];
rz(-0.42492953) q[3];
sx q[3];
rz(0.4900842) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77077615) q[0];
sx q[0];
rz(-0.66148615) q[0];
sx q[0];
rz(0.01195512) q[0];
rz(-1.6318343) q[1];
sx q[1];
rz(-0.44523528) q[1];
sx q[1];
rz(0.59648046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79967989) q[0];
sx q[0];
rz(-1.1014859) q[0];
sx q[0];
rz(1.5921568) q[0];
rz(-1.8305186) q[2];
sx q[2];
rz(-1.9225645) q[2];
sx q[2];
rz(3.0992257) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0833019) q[1];
sx q[1];
rz(-1.4920248) q[1];
sx q[1];
rz(2.0473256) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29774547) q[3];
sx q[3];
rz(-1.4199054) q[3];
sx q[3];
rz(-2.6835364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.163588) q[2];
sx q[2];
rz(-2.1808193) q[2];
sx q[2];
rz(0.19700024) q[2];
rz(1.9893076) q[3];
sx q[3];
rz(-0.14324337) q[3];
sx q[3];
rz(0.44767374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.2780509) q[0];
sx q[0];
rz(-1.0095689) q[0];
sx q[0];
rz(2.9472245) q[0];
rz(-2.9327783) q[1];
sx q[1];
rz(-1.6046235) q[1];
sx q[1];
rz(-1.0135289) q[1];
rz(-2.2603358) q[2];
sx q[2];
rz(-1.643558) q[2];
sx q[2];
rz(0.087842077) q[2];
rz(0.50894751) q[3];
sx q[3];
rz(-1.043368) q[3];
sx q[3];
rz(2.1989078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
