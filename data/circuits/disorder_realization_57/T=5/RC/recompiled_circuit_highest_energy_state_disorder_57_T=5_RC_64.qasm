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
rz(-0.34173319) q[1];
sx q[1];
rz(-2.3050397) q[1];
sx q[1];
rz(0.41681448) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57274103) q[0];
sx q[0];
rz(-2.4089455) q[0];
sx q[0];
rz(-0.1591263) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55617638) q[2];
sx q[2];
rz(-1.3601298) q[2];
sx q[2];
rz(0.81708252) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3089247) q[1];
sx q[1];
rz(-1.0005669) q[1];
sx q[1];
rz(0.69242386) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5384212) q[3];
sx q[3];
rz(-1.7322632) q[3];
sx q[3];
rz(0.77154713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1699528) q[2];
sx q[2];
rz(-1.4687186) q[2];
sx q[2];
rz(-1.8876342) q[2];
rz(1.5597255) q[3];
sx q[3];
rz(-0.73837787) q[3];
sx q[3];
rz(-2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.1935683) q[0];
sx q[0];
rz(-1.7727611) q[0];
sx q[0];
rz(2.3728306) q[0];
rz(-1.7747152) q[1];
sx q[1];
rz(-1.8194852) q[1];
sx q[1];
rz(1.0381402) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067613451) q[0];
sx q[0];
rz(-1.5609976) q[0];
sx q[0];
rz(3.1329747) q[0];
rz(-2.7323167) q[2];
sx q[2];
rz(-0.76414327) q[2];
sx q[2];
rz(2.0913569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8703192) q[1];
sx q[1];
rz(-2.4044577) q[1];
sx q[1];
rz(-2.5802274) q[1];
rz(1.6065381) q[3];
sx q[3];
rz(-0.48887353) q[3];
sx q[3];
rz(-1.8349748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.495503) q[2];
sx q[2];
rz(-1.8547408) q[2];
sx q[2];
rz(-0.82026473) q[2];
rz(-2.8881554) q[3];
sx q[3];
rz(-2.7866252) q[3];
sx q[3];
rz(-2.7804815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28061098) q[0];
sx q[0];
rz(-0.51787037) q[0];
sx q[0];
rz(0.88731998) q[0];
rz(0.53030983) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(-1.1393772) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47592012) q[0];
sx q[0];
rz(-2.7662656) q[0];
sx q[0];
rz(-0.39552839) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33311756) q[2];
sx q[2];
rz(-1.0651759) q[2];
sx q[2];
rz(-2.5207455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2816086) q[1];
sx q[1];
rz(-1.8164595) q[1];
sx q[1];
rz(-1.5857693) q[1];
rz(3.0501796) q[3];
sx q[3];
rz(-0.52052906) q[3];
sx q[3];
rz(0.5239858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.016971074) q[2];
sx q[2];
rz(-1.6861702) q[2];
sx q[2];
rz(2.7023081) q[2];
rz(2.6708421) q[3];
sx q[3];
rz(-3.0622523) q[3];
sx q[3];
rz(1.144217) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73862326) q[0];
sx q[0];
rz(-2.2082177) q[0];
sx q[0];
rz(-0.11548197) q[0];
rz(-0.11784095) q[1];
sx q[1];
rz(-1.203048) q[1];
sx q[1];
rz(-1.3062564) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3958493) q[0];
sx q[0];
rz(-1.9877399) q[0];
sx q[0];
rz(2.8265619) q[0];
x q[1];
rz(-1.9491862) q[2];
sx q[2];
rz(-0.94589627) q[2];
sx q[2];
rz(-1.6922127) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0769121) q[1];
sx q[1];
rz(-2.2222328) q[1];
sx q[1];
rz(-0.76131911) q[1];
rz(3.066411) q[3];
sx q[3];
rz(-0.99937056) q[3];
sx q[3];
rz(-2.927185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8370886) q[2];
sx q[2];
rz(-1.3136761) q[2];
sx q[2];
rz(-2.9869249) q[2];
rz(-1.6020487) q[3];
sx q[3];
rz(-1.3402901) q[3];
sx q[3];
rz(-2.535517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.57347572) q[0];
sx q[0];
rz(-1.0433759) q[0];
sx q[0];
rz(2.306275) q[0];
rz(-2.4256445) q[1];
sx q[1];
rz(-2.7138111) q[1];
sx q[1];
rz(-0.11944019) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7258541) q[0];
sx q[0];
rz(-1.5149759) q[0];
sx q[0];
rz(-3.124696) q[0];
rz(-2.832069) q[2];
sx q[2];
rz(-2.672827) q[2];
sx q[2];
rz(-0.13407198) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.24158289) q[1];
sx q[1];
rz(-2.0323583) q[1];
sx q[1];
rz(-2.475283) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9914845) q[3];
sx q[3];
rz(-1.0771695) q[3];
sx q[3];
rz(-2.8378262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2401838) q[2];
sx q[2];
rz(-2.8808424) q[2];
sx q[2];
rz(0.060997941) q[2];
rz(-3.0933464) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(0.72859305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7399087) q[0];
sx q[0];
rz(-0.70867276) q[0];
sx q[0];
rz(-2.0106864) q[0];
rz(2.6629958) q[1];
sx q[1];
rz(-0.60239783) q[1];
sx q[1];
rz(-1.7344249) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3815752) q[0];
sx q[0];
rz(-1.572346) q[0];
sx q[0];
rz(-3.0880307) q[0];
rz(-pi) q[1];
rz(-0.28193177) q[2];
sx q[2];
rz(-2.273186) q[2];
sx q[2];
rz(1.7676181) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6428522) q[1];
sx q[1];
rz(-1.7987218) q[1];
sx q[1];
rz(-3.0251316) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70194093) q[3];
sx q[3];
rz(-2.8656368) q[3];
sx q[3];
rz(-0.70358932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56167928) q[2];
sx q[2];
rz(-2.734197) q[2];
sx q[2];
rz(-1.178406) q[2];
rz(2.0922349) q[3];
sx q[3];
rz(-2.6047843) q[3];
sx q[3];
rz(-1.8572846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2105763) q[0];
sx q[0];
rz(-1.1897621) q[0];
sx q[0];
rz(1.3354906) q[0];
rz(1.8203075) q[1];
sx q[1];
rz(-1.0711461) q[1];
sx q[1];
rz(2.8335422) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0548693) q[0];
sx q[0];
rz(-0.44741524) q[0];
sx q[0];
rz(2.6347876) q[0];
rz(2.5045583) q[2];
sx q[2];
rz(-2.5278628) q[2];
sx q[2];
rz(-1.5368652) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1738893) q[1];
sx q[1];
rz(-0.99882616) q[1];
sx q[1];
rz(-2.553945) q[1];
rz(3.0940476) q[3];
sx q[3];
rz(-1.9564087) q[3];
sx q[3];
rz(-1.4908718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41219741) q[2];
sx q[2];
rz(-2.2008379) q[2];
sx q[2];
rz(-3.0103053) q[2];
rz(-0.73733759) q[3];
sx q[3];
rz(-1.3056583) q[3];
sx q[3];
rz(2.9837515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103631) q[0];
sx q[0];
rz(-1.0897626) q[0];
sx q[0];
rz(-2.6112153) q[0];
rz(-1.7165548) q[1];
sx q[1];
rz(-2.0030463) q[1];
sx q[1];
rz(2.4200965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58396858) q[0];
sx q[0];
rz(-0.97088214) q[0];
sx q[0];
rz(-2.1198089) q[0];
rz(-0.28571911) q[2];
sx q[2];
rz(-0.48595324) q[2];
sx q[2];
rz(0.0067809502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87807314) q[1];
sx q[1];
rz(-1.7612447) q[1];
sx q[1];
rz(-1.3506804) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0305392) q[3];
sx q[3];
rz(-2.3559582) q[3];
sx q[3];
rz(-1.812404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7044652) q[2];
sx q[2];
rz(-0.4054873) q[2];
sx q[2];
rz(-1.6638157) q[2];
rz(-1.9780698) q[3];
sx q[3];
rz(-2.0587557) q[3];
sx q[3];
rz(1.5014974) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8116233) q[0];
sx q[0];
rz(-1.7003308) q[0];
sx q[0];
rz(-0.60923088) q[0];
rz(1.5513264) q[1];
sx q[1];
rz(-2.8158999) q[1];
sx q[1];
rz(-1.4564266) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0349755) q[0];
sx q[0];
rz(-1.2143232) q[0];
sx q[0];
rz(0.26654213) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6316192) q[2];
sx q[2];
rz(-1.0863388) q[2];
sx q[2];
rz(1.6866419) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2472154) q[1];
sx q[1];
rz(-0.4900107) q[1];
sx q[1];
rz(1.2213329) q[1];
x q[2];
rz(-0.88202694) q[3];
sx q[3];
rz(-2.2891364) q[3];
sx q[3];
rz(0.73757987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7234601) q[2];
sx q[2];
rz(-0.89833608) q[2];
sx q[2];
rz(2.2612803) q[2];
rz(2.9355925) q[3];
sx q[3];
rz(-0.42492953) q[3];
sx q[3];
rz(-2.6515085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.3708165) q[0];
sx q[0];
rz(-0.66148615) q[0];
sx q[0];
rz(-3.1296375) q[0];
rz(1.5097584) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(2.5451122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3801381) q[0];
sx q[0];
rz(-1.5517457) q[0];
sx q[0];
rz(-2.6721902) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.311074) q[2];
sx q[2];
rz(-1.9225645) q[2];
sx q[2];
rz(-3.0992257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0582907) q[1];
sx q[1];
rz(-1.6495678) q[1];
sx q[1];
rz(-1.0942671) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4130728) q[3];
sx q[3];
rz(-1.2765358) q[3];
sx q[3];
rz(2.074948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97800469) q[2];
sx q[2];
rz(-2.1808193) q[2];
sx q[2];
rz(-0.19700024) q[2];
rz(1.9893076) q[3];
sx q[3];
rz(-2.9983493) q[3];
sx q[3];
rz(-0.44767374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86354179) q[0];
sx q[0];
rz(-2.1320237) q[0];
sx q[0];
rz(-0.19436819) q[0];
rz(2.9327783) q[1];
sx q[1];
rz(-1.5369692) q[1];
sx q[1];
rz(2.1280638) q[1];
rz(-3.0473982) q[2];
sx q[2];
rz(-2.258156) q[2];
sx q[2];
rz(-1.5428262) q[2];
rz(2.6326451) q[3];
sx q[3];
rz(-2.0982246) q[3];
sx q[3];
rz(-0.94268483) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
