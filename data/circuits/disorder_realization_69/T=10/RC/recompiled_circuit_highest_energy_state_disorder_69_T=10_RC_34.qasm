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
rz(-2.3186853) q[0];
sx q[0];
rz(3.5003852) q[0];
sx q[0];
rz(8.556463) q[0];
rz(-1.4153642) q[1];
sx q[1];
rz(-1.1338898) q[1];
sx q[1];
rz(2.147832) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.225634) q[0];
sx q[0];
rz(-1.4580112) q[0];
sx q[0];
rz(-2.7382572) q[0];
x q[1];
rz(-1.1311985) q[2];
sx q[2];
rz(-1.1102668) q[2];
sx q[2];
rz(1.4714413) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43683166) q[1];
sx q[1];
rz(-1.1306354) q[1];
sx q[1];
rz(2.6635567) q[1];
rz(-1.1339784) q[3];
sx q[3];
rz(-1.7941495) q[3];
sx q[3];
rz(-2.2586266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22780861) q[2];
sx q[2];
rz(-1.8081534) q[2];
sx q[2];
rz(-2.559973) q[2];
rz(2.4070814) q[3];
sx q[3];
rz(-1.651265) q[3];
sx q[3];
rz(0.41828004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89473474) q[0];
sx q[0];
rz(-1.7623836) q[0];
sx q[0];
rz(-0.49767622) q[0];
rz(1.0408164) q[1];
sx q[1];
rz(-2.7383995) q[1];
sx q[1];
rz(-1.9042447) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055641551) q[0];
sx q[0];
rz(-2.4179672) q[0];
sx q[0];
rz(3.0776204) q[0];
x q[1];
rz(-1.2096268) q[2];
sx q[2];
rz(-1.3945082) q[2];
sx q[2];
rz(0.90106264) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1488545) q[1];
sx q[1];
rz(-1.9991367) q[1];
sx q[1];
rz(-0.15495877) q[1];
rz(-pi) q[2];
rz(-2.4604843) q[3];
sx q[3];
rz(-2.8048385) q[3];
sx q[3];
rz(0.82267534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.70025468) q[2];
sx q[2];
rz(-2.4105218) q[2];
sx q[2];
rz(-0.31044427) q[2];
rz(-1.1566409) q[3];
sx q[3];
rz(-1.1247331) q[3];
sx q[3];
rz(-1.9748851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0079086) q[0];
sx q[0];
rz(-2.5569361) q[0];
sx q[0];
rz(0.34580082) q[0];
rz(-2.8969104) q[1];
sx q[1];
rz(-0.9318277) q[1];
sx q[1];
rz(1.7049047) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0999683) q[0];
sx q[0];
rz(-1.9910553) q[0];
sx q[0];
rz(-0.34456518) q[0];
rz(-pi) q[1];
rz(2.56836) q[2];
sx q[2];
rz(-0.65224183) q[2];
sx q[2];
rz(-2.531372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50592283) q[1];
sx q[1];
rz(-1.782592) q[1];
sx q[1];
rz(2.3426272) q[1];
rz(-1.9448024) q[3];
sx q[3];
rz(-2.0197649) q[3];
sx q[3];
rz(-1.6351007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1338542) q[2];
sx q[2];
rz(-1.1275007) q[2];
sx q[2];
rz(-2.5109042) q[2];
rz(0.33356365) q[3];
sx q[3];
rz(-2.1400698) q[3];
sx q[3];
rz(1.001531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.075060464) q[0];
sx q[0];
rz(-0.94828951) q[0];
sx q[0];
rz(-1.4045658) q[0];
rz(0.36918494) q[1];
sx q[1];
rz(-1.7351979) q[1];
sx q[1];
rz(-3.0558443) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8456025) q[0];
sx q[0];
rz(-1.9787425) q[0];
sx q[0];
rz(1.8970117) q[0];
x q[1];
rz(-1.9783114) q[2];
sx q[2];
rz(-1.9980944) q[2];
sx q[2];
rz(-2.3321762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4032674) q[1];
sx q[1];
rz(-1.333263) q[1];
sx q[1];
rz(-2.0473785) q[1];
x q[2];
rz(-1.2530009) q[3];
sx q[3];
rz(-1.9488584) q[3];
sx q[3];
rz(1.2757358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3021476) q[2];
sx q[2];
rz(-1.9044694) q[2];
sx q[2];
rz(0.5298003) q[2];
rz(3.1032622) q[3];
sx q[3];
rz(-0.72961346) q[3];
sx q[3];
rz(-1.5420325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2235276) q[0];
sx q[0];
rz(-0.46850884) q[0];
sx q[0];
rz(1.9388306) q[0];
rz(0.27944061) q[1];
sx q[1];
rz(-1.025082) q[1];
sx q[1];
rz(2.3675945) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1647332) q[0];
sx q[0];
rz(-0.79553878) q[0];
sx q[0];
rz(1.1245994) q[0];
rz(-pi) q[1];
rz(-0.33542893) q[2];
sx q[2];
rz(-2.895439) q[2];
sx q[2];
rz(-0.85573643) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8786278) q[1];
sx q[1];
rz(-1.557918) q[1];
sx q[1];
rz(-1.4317382) q[1];
x q[2];
rz(-1.0353885) q[3];
sx q[3];
rz(-1.5683953) q[3];
sx q[3];
rz(-2.6794499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1191795) q[2];
sx q[2];
rz(-3.000562) q[2];
sx q[2];
rz(0.5640344) q[2];
rz(0.65308475) q[3];
sx q[3];
rz(-1.1354732) q[3];
sx q[3];
rz(-2.691332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3611203) q[0];
sx q[0];
rz(-0.55279624) q[0];
sx q[0];
rz(-2.4984388) q[0];
rz(1.9505352) q[1];
sx q[1];
rz(-1.4570313) q[1];
sx q[1];
rz(-2.4868884) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768652) q[0];
sx q[0];
rz(-0.18023364) q[0];
sx q[0];
rz(-0.36156543) q[0];
x q[1];
rz(0.68393647) q[2];
sx q[2];
rz(-1.5412871) q[2];
sx q[2];
rz(-2.3060407) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5712329) q[1];
sx q[1];
rz(-2.1132937) q[1];
sx q[1];
rz(-0.076537655) q[1];
rz(-1.7406171) q[3];
sx q[3];
rz(-0.39547503) q[3];
sx q[3];
rz(2.2251289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6341298) q[2];
sx q[2];
rz(-0.37709388) q[2];
sx q[2];
rz(1.4748352) q[2];
rz(0.48464388) q[3];
sx q[3];
rz(-2.1438997) q[3];
sx q[3];
rz(1.8528329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5019048) q[0];
sx q[0];
rz(-0.26308331) q[0];
sx q[0];
rz(-0.68341533) q[0];
rz(3.0112093) q[1];
sx q[1];
rz(-1.5510635) q[1];
sx q[1];
rz(-0.032616671) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0158247) q[0];
sx q[0];
rz(-0.59249632) q[0];
sx q[0];
rz(-1.8886719) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7702297) q[2];
sx q[2];
rz(-0.89493361) q[2];
sx q[2];
rz(0.95914721) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72717818) q[1];
sx q[1];
rz(-1.4041931) q[1];
sx q[1];
rz(1.8397306) q[1];
rz(-pi) q[2];
rz(1.2682876) q[3];
sx q[3];
rz(-0.72495809) q[3];
sx q[3];
rz(2.8444949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42328295) q[2];
sx q[2];
rz(-1.1085359) q[2];
sx q[2];
rz(-2.5763467) q[2];
rz(1.7806753) q[3];
sx q[3];
rz(-2.8748685) q[3];
sx q[3];
rz(-0.81418532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3703506) q[0];
sx q[0];
rz(-3.0841565) q[0];
sx q[0];
rz(-2.8420319) q[0];
rz(-1.4467422) q[1];
sx q[1];
rz(-0.87211496) q[1];
sx q[1];
rz(2.1902693) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.882269) q[0];
sx q[0];
rz(-2.933942) q[0];
sx q[0];
rz(-1.4126247) q[0];
x q[1];
rz(-1.5278242) q[2];
sx q[2];
rz(-2.9298721) q[2];
sx q[2];
rz(2.4076902) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.42859498) q[1];
sx q[1];
rz(-0.42474898) q[1];
sx q[1];
rz(0.049475706) q[1];
rz(-1.3090012) q[3];
sx q[3];
rz(-1.6556532) q[3];
sx q[3];
rz(-3.121162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99913725) q[2];
sx q[2];
rz(-0.74644011) q[2];
sx q[2];
rz(2.8738521) q[2];
rz(-1.0033876) q[3];
sx q[3];
rz(-2.181874) q[3];
sx q[3];
rz(2.3333534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
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
rz(0.35850152) q[0];
sx q[0];
rz(-2.6434904) q[0];
sx q[0];
rz(0.95712334) q[0];
rz(2.762291) q[1];
sx q[1];
rz(-0.4117659) q[1];
sx q[1];
rz(-2.9764825) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9921761) q[0];
sx q[0];
rz(-2.3603129) q[0];
sx q[0];
rz(0.4968471) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4719285) q[2];
sx q[2];
rz(-2.6672088) q[2];
sx q[2];
rz(2.6883467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8891661) q[1];
sx q[1];
rz(-2.5339087) q[1];
sx q[1];
rz(-1.8914188) q[1];
rz(-1.9103138) q[3];
sx q[3];
rz(-0.63106189) q[3];
sx q[3];
rz(-1.5724044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2740606) q[2];
sx q[2];
rz(-1.910285) q[2];
sx q[2];
rz(-0.77862281) q[2];
rz(2.8042931) q[3];
sx q[3];
rz(-1.0236579) q[3];
sx q[3];
rz(1.5571099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.857665) q[0];
sx q[0];
rz(-1.090467) q[0];
sx q[0];
rz(2.0528059) q[0];
rz(-0.42824832) q[1];
sx q[1];
rz(-2.1183522) q[1];
sx q[1];
rz(-0.64839378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26864949) q[0];
sx q[0];
rz(-2.6336484) q[0];
sx q[0];
rz(2.5769039) q[0];
x q[1];
rz(-0.044610046) q[2];
sx q[2];
rz(-1.9084435) q[2];
sx q[2];
rz(-2.5672947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9251483) q[1];
sx q[1];
rz(-2.2183462) q[1];
sx q[1];
rz(-1.5179894) q[1];
rz(-1.6850059) q[3];
sx q[3];
rz(-1.0100967) q[3];
sx q[3];
rz(1.4637092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8757561) q[2];
sx q[2];
rz(-1.0886322) q[2];
sx q[2];
rz(-0.17975532) q[2];
rz(-1.9836551) q[3];
sx q[3];
rz(-1.4561184) q[3];
sx q[3];
rz(1.2402844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801878) q[0];
sx q[0];
rz(-2.9869933) q[0];
sx q[0];
rz(-2.3416478) q[0];
rz(-1.1663306) q[1];
sx q[1];
rz(-0.98465289) q[1];
sx q[1];
rz(-2.226895) q[1];
rz(-0.24576743) q[2];
sx q[2];
rz(-1.859833) q[2];
sx q[2];
rz(-1.8219994) q[2];
rz(-2.5869681) q[3];
sx q[3];
rz(-2.7477874) q[3];
sx q[3];
rz(1.2870233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
