OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6405606) q[0];
sx q[0];
rz(-0.80067331) q[0];
sx q[0];
rz(-0.14525695) q[0];
rz(-0.89291209) q[1];
sx q[1];
rz(3.555759) q[1];
sx q[1];
rz(10.531737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1205638) q[0];
sx q[0];
rz(-1.0146838) q[0];
sx q[0];
rz(-0.70379852) q[0];
rz(0.85914454) q[2];
sx q[2];
rz(-2.7224053) q[2];
sx q[2];
rz(2.223658) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8095137) q[1];
sx q[1];
rz(-2.4831536) q[1];
sx q[1];
rz(2.0488942) q[1];
rz(-pi) q[2];
rz(0.21671076) q[3];
sx q[3];
rz(-1.9418849) q[3];
sx q[3];
rz(-2.3909274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60078159) q[2];
sx q[2];
rz(-1.3593707) q[2];
sx q[2];
rz(-0.092279807) q[2];
rz(-1.4394834) q[3];
sx q[3];
rz(-1.9974134) q[3];
sx q[3];
rz(-2.2030742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0094902078) q[0];
sx q[0];
rz(-1.1630031) q[0];
sx q[0];
rz(-0.60390419) q[0];
rz(1.6101135) q[1];
sx q[1];
rz(-1.8096626) q[1];
sx q[1];
rz(1.5276705) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92787111) q[0];
sx q[0];
rz(-1.8851213) q[0];
sx q[0];
rz(2.1840591) q[0];
rz(0.98601922) q[2];
sx q[2];
rz(-0.9169609) q[2];
sx q[2];
rz(-0.80160917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.70992596) q[1];
sx q[1];
rz(-1.5824236) q[1];
sx q[1];
rz(-2.8678368) q[1];
rz(0.0068596938) q[3];
sx q[3];
rz(-0.09819542) q[3];
sx q[3];
rz(-2.6650037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4332726) q[2];
sx q[2];
rz(-1.4488139) q[2];
sx q[2];
rz(-1.0478919) q[2];
rz(-0.6692872) q[3];
sx q[3];
rz(-0.71988121) q[3];
sx q[3];
rz(1.1650813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79353756) q[0];
sx q[0];
rz(-2.3568643) q[0];
sx q[0];
rz(-2.134557) q[0];
rz(-0.29193613) q[1];
sx q[1];
rz(-2.597229) q[1];
sx q[1];
rz(0.67063355) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48146433) q[0];
sx q[0];
rz(-0.95230904) q[0];
sx q[0];
rz(3.0342558) q[0];
x q[1];
rz(0.67523662) q[2];
sx q[2];
rz(-0.91885447) q[2];
sx q[2];
rz(2.3777131) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.933316) q[1];
sx q[1];
rz(-1.2348935) q[1];
sx q[1];
rz(-2.4927054) q[1];
rz(-pi) q[2];
rz(0.11588736) q[3];
sx q[3];
rz(-1.8975583) q[3];
sx q[3];
rz(2.9123757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2900419) q[2];
sx q[2];
rz(-2.2883577) q[2];
sx q[2];
rz(-2.2500989) q[2];
rz(-1.1024891) q[3];
sx q[3];
rz(-2.1474371) q[3];
sx q[3];
rz(1.4055143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3201228) q[0];
sx q[0];
rz(-1.7550884) q[0];
sx q[0];
rz(2.084305) q[0];
rz(0.0013466324) q[1];
sx q[1];
rz(-0.82095447) q[1];
sx q[1];
rz(0.79426208) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43633902) q[0];
sx q[0];
rz(-0.90898517) q[0];
sx q[0];
rz(1.5940985) q[0];
rz(-pi) q[1];
rz(0.67790548) q[2];
sx q[2];
rz(-2.6982582) q[2];
sx q[2];
rz(-1.6010798) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9201607) q[1];
sx q[1];
rz(-0.86593141) q[1];
sx q[1];
rz(-0.95734289) q[1];
rz(-pi) q[2];
rz(-0.36093386) q[3];
sx q[3];
rz(-0.31692255) q[3];
sx q[3];
rz(1.288687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1069676) q[2];
sx q[2];
rz(-0.62109533) q[2];
sx q[2];
rz(-2.1194439) q[2];
rz(-1.243783) q[3];
sx q[3];
rz(-0.94977489) q[3];
sx q[3];
rz(1.7826084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6720402) q[0];
sx q[0];
rz(-1.7615027) q[0];
sx q[0];
rz(-0.45355466) q[0];
rz(2.1039311) q[1];
sx q[1];
rz(-1.1353759) q[1];
sx q[1];
rz(0.79197788) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29402367) q[0];
sx q[0];
rz(-1.5269465) q[0];
sx q[0];
rz(1.8490318) q[0];
rz(2.1395965) q[2];
sx q[2];
rz(-0.65868568) q[2];
sx q[2];
rz(1.9994761) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9336613) q[1];
sx q[1];
rz(-1.232389) q[1];
sx q[1];
rz(1.4470571) q[1];
x q[2];
rz(1.4384934) q[3];
sx q[3];
rz(-0.13935329) q[3];
sx q[3];
rz(2.9840368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3114634) q[2];
sx q[2];
rz(-0.66581231) q[2];
sx q[2];
rz(-0.30544454) q[2];
rz(0.18181248) q[3];
sx q[3];
rz(-1.5287377) q[3];
sx q[3];
rz(2.4055433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.55564725) q[0];
sx q[0];
rz(-0.82252994) q[0];
sx q[0];
rz(0.12538759) q[0];
rz(1.5749982) q[1];
sx q[1];
rz(-1.4488723) q[1];
sx q[1];
rz(-0.032141846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9956995) q[0];
sx q[0];
rz(-2.3902378) q[0];
sx q[0];
rz(0.32352792) q[0];
rz(1.6400385) q[2];
sx q[2];
rz(-2.4089795) q[2];
sx q[2];
rz(-2.46012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.39592) q[1];
sx q[1];
rz(-2.2587639) q[1];
sx q[1];
rz(-2.2017971) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6866845) q[3];
sx q[3];
rz(-0.82852302) q[3];
sx q[3];
rz(-2.7477086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0142168) q[2];
sx q[2];
rz(-1.2140423) q[2];
sx q[2];
rz(1.1227013) q[2];
rz(0.073401062) q[3];
sx q[3];
rz(-0.95728907) q[3];
sx q[3];
rz(-2.5286123) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4323394) q[0];
sx q[0];
rz(-1.4839577) q[0];
sx q[0];
rz(-2.6313229) q[0];
rz(0.36422745) q[1];
sx q[1];
rz(-2.7204456) q[1];
sx q[1];
rz(1.4998923) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9534665) q[0];
sx q[0];
rz(-2.0333485) q[0];
sx q[0];
rz(-3.0234087) q[0];
rz(-pi) q[1];
rz(1.4128311) q[2];
sx q[2];
rz(-1.6876432) q[2];
sx q[2];
rz(-0.37010461) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0081029) q[1];
sx q[1];
rz(-1.7945707) q[1];
sx q[1];
rz(2.6340805) q[1];
x q[2];
rz(-0.73697258) q[3];
sx q[3];
rz(-1.7993357) q[3];
sx q[3];
rz(-1.6933954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.69786543) q[2];
sx q[2];
rz(-1.5865734) q[2];
sx q[2];
rz(-2.2743684) q[2];
rz(-1.5602268) q[3];
sx q[3];
rz(-2.9428704) q[3];
sx q[3];
rz(-1.8038512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7241868) q[0];
sx q[0];
rz(-0.066411821) q[0];
sx q[0];
rz(-2.7253286) q[0];
rz(1.9789713) q[1];
sx q[1];
rz(-1.1888209) q[1];
sx q[1];
rz(0.75688854) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1838186) q[0];
sx q[0];
rz(-1.2831492) q[0];
sx q[0];
rz(-1.1270866) q[0];
rz(0.9932809) q[2];
sx q[2];
rz(-1.2009635) q[2];
sx q[2];
rz(0.12060697) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6357506) q[1];
sx q[1];
rz(-1.381784) q[1];
sx q[1];
rz(-1.1357978) q[1];
x q[2];
rz(-2.5439436) q[3];
sx q[3];
rz(-0.72859287) q[3];
sx q[3];
rz(0.33580175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4550712) q[2];
sx q[2];
rz(-1.4342118) q[2];
sx q[2];
rz(-1.3580458) q[2];
rz(-2.3250735) q[3];
sx q[3];
rz(-1.2314545) q[3];
sx q[3];
rz(-0.16286287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002976) q[0];
sx q[0];
rz(-1.6930027) q[0];
sx q[0];
rz(1.7342389) q[0];
rz(0.74527144) q[1];
sx q[1];
rz(-0.73423568) q[1];
sx q[1];
rz(-1.5819246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9055515) q[0];
sx q[0];
rz(-2.2187382) q[0];
sx q[0];
rz(-2.7929162) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4938131) q[2];
sx q[2];
rz(-1.544892) q[2];
sx q[2];
rz(1.6960953) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.24589495) q[1];
sx q[1];
rz(-2.0132397) q[1];
sx q[1];
rz(2.9986283) q[1];
rz(-pi) q[2];
rz(-1.6635062) q[3];
sx q[3];
rz(-2.5639236) q[3];
sx q[3];
rz(-2.3232164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4173296) q[2];
sx q[2];
rz(-1.686692) q[2];
sx q[2];
rz(1.1154741) q[2];
rz(-2.8880902) q[3];
sx q[3];
rz(-1.8799672) q[3];
sx q[3];
rz(-0.0089664627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.020092) q[0];
sx q[0];
rz(-1.2637063) q[0];
sx q[0];
rz(2.3790835) q[0];
rz(2.1992042) q[1];
sx q[1];
rz(-2.0768879) q[1];
sx q[1];
rz(1.2368088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13888966) q[0];
sx q[0];
rz(-1.6322487) q[0];
sx q[0];
rz(2.6451254) q[0];
rz(2.5618844) q[2];
sx q[2];
rz(-2.1695608) q[2];
sx q[2];
rz(2.6331462) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2357465) q[1];
sx q[1];
rz(-1.5457834) q[1];
sx q[1];
rz(-0.025084875) q[1];
x q[2];
rz(-3.0701261) q[3];
sx q[3];
rz(-1.9744919) q[3];
sx q[3];
rz(1.57406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8952055) q[2];
sx q[2];
rz(-0.095194101) q[2];
sx q[2];
rz(0.0072366317) q[2];
rz(2.9712408) q[3];
sx q[3];
rz(-1.6536313) q[3];
sx q[3];
rz(0.65792221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8728747) q[0];
sx q[0];
rz(-1.3997411) q[0];
sx q[0];
rz(1.1135143) q[0];
rz(0.089182236) q[1];
sx q[1];
rz(-1.5166278) q[1];
sx q[1];
rz(1.2022432) q[1];
rz(-2.4287672) q[2];
sx q[2];
rz(-0.20607866) q[2];
sx q[2];
rz(-1.1417749) q[2];
rz(-0.30059697) q[3];
sx q[3];
rz(-1.3924122) q[3];
sx q[3];
rz(-2.198749) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
