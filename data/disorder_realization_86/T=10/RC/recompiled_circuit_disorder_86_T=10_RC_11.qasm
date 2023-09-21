OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(1.4120742) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(1.6593978) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94443653) q[0];
sx q[0];
rz(-1.4154139) q[0];
sx q[0];
rz(-2.7128501) q[0];
rz(-2.8843845) q[2];
sx q[2];
rz(-1.6991985) q[2];
sx q[2];
rz(0.53127015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9285779) q[1];
sx q[1];
rz(-1.9933812) q[1];
sx q[1];
rz(-1.0282474) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0406038) q[3];
sx q[3];
rz(-1.0102934) q[3];
sx q[3];
rz(1.8141754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1564864) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(-2.1872897) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99825478) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(1.5401309) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-2.1781133) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.119495) q[0];
sx q[0];
rz(-0.95675981) q[0];
sx q[0];
rz(-3.1387781) q[0];
rz(-2.1396779) q[2];
sx q[2];
rz(-1.8520253) q[2];
sx q[2];
rz(0.98904726) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3669489) q[1];
sx q[1];
rz(-2.0683214) q[1];
sx q[1];
rz(2.7760387) q[1];
x q[2];
rz(0.24641896) q[3];
sx q[3];
rz(-1.2296457) q[3];
sx q[3];
rz(0.44647549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5144689) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(-2.3965805) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(-2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2117675) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(0.79743687) q[0];
rz(2.0939317) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(-0.55999666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13941923) q[0];
sx q[0];
rz(-1.1886485) q[0];
sx q[0];
rz(-2.9234773) q[0];
rz(-2.838344) q[2];
sx q[2];
rz(-1.5911284) q[2];
sx q[2];
rz(-0.081239935) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3126038) q[1];
sx q[1];
rz(-1.3723515) q[1];
sx q[1];
rz(2.7752084) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47645724) q[3];
sx q[3];
rz(-2.0072848) q[3];
sx q[3];
rz(-0.95526327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.75227633) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(0.28451434) q[0];
rz(-2.8248887) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(-1.2987312) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1176778) q[0];
sx q[0];
rz(-1.0751343) q[0];
sx q[0];
rz(-0.25650521) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3372075) q[2];
sx q[2];
rz(-1.5717236) q[2];
sx q[2];
rz(-1.5915807) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97949308) q[1];
sx q[1];
rz(-0.59826189) q[1];
sx q[1];
rz(1.23566) q[1];
x q[2];
rz(1.2508568) q[3];
sx q[3];
rz(-0.3443998) q[3];
sx q[3];
rz(-2.875945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(-0.38468012) q[2];
rz(-2.3875333) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(-1.4543021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5383179) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(-1.3866562) q[0];
rz(2.9105913) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(-0.2968266) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41243991) q[0];
sx q[0];
rz(-1.6099596) q[0];
sx q[0];
rz(-0.57106437) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7237687) q[2];
sx q[2];
rz(-2.3077871) q[2];
sx q[2];
rz(-1.0066777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8639212) q[1];
sx q[1];
rz(-1.9660945) q[1];
sx q[1];
rz(-1.0635832) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6597219) q[3];
sx q[3];
rz(-0.85907912) q[3];
sx q[3];
rz(-2.5951648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(-2.7491167) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0157938) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(1.8136576) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(-0.60633916) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7327001) q[0];
sx q[0];
rz(-0.048763976) q[0];
sx q[0];
rz(-2.7932037) q[0];
rz(-2.7251284) q[2];
sx q[2];
rz(-2.205924) q[2];
sx q[2];
rz(-0.093402775) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50748435) q[1];
sx q[1];
rz(-2.1180696) q[1];
sx q[1];
rz(0.21168153) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.074999853) q[3];
sx q[3];
rz(-1.787775) q[3];
sx q[3];
rz(-2.0892339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63885826) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(1.139337) q[2];
rz(-1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(3.0373354) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6095603) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(-1.5628901) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(0.79024822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9433141) q[0];
sx q[0];
rz(-0.6753079) q[0];
sx q[0];
rz(-0.4839464) q[0];
rz(-pi) q[1];
rz(3.087567) q[2];
sx q[2];
rz(-2.9276491) q[2];
sx q[2];
rz(-2.8097048) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0223169) q[1];
sx q[1];
rz(-1.9609309) q[1];
sx q[1];
rz(-1.2783865) q[1];
rz(-pi) q[2];
rz(0.4831794) q[3];
sx q[3];
rz(-2.4348767) q[3];
sx q[3];
rz(0.28373517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(1.7283758) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(-0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(2.5324902) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1202639) q[0];
sx q[0];
rz(-1.6006032) q[0];
sx q[0];
rz(-1.0111615) q[0];
rz(2.6128204) q[2];
sx q[2];
rz(-1.9430338) q[2];
sx q[2];
rz(-2.1793274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3125004) q[1];
sx q[1];
rz(-1.6626076) q[1];
sx q[1];
rz(-1.5593668) q[1];
x q[2];
rz(-0.96447585) q[3];
sx q[3];
rz(-1.8165605) q[3];
sx q[3];
rz(2.5442459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(0.88225538) q[2];
rz(-1.7404209) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700478) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(-1.7154988) q[0];
rz(-3.0601314) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(0.55823278) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0575858) q[0];
sx q[0];
rz(-1.9965729) q[0];
sx q[0];
rz(0.42752479) q[0];
x q[1];
rz(-1.4633281) q[2];
sx q[2];
rz(-1.5170013) q[2];
sx q[2];
rz(-1.4449643) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3840752) q[1];
sx q[1];
rz(-1.5069403) q[1];
sx q[1];
rz(-2.3149895) q[1];
rz(-pi) q[2];
rz(-1.3853119) q[3];
sx q[3];
rz(-1.0914601) q[3];
sx q[3];
rz(2.6244147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(-1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.50487173) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(1.6037534) q[0];
rz(0.82540712) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(2.6182981) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1062746) q[0];
sx q[0];
rz(-2.1786852) q[0];
sx q[0];
rz(-1.4950698) q[0];
x q[1];
rz(-0.043838219) q[2];
sx q[2];
rz(-2.1423116) q[2];
sx q[2];
rz(0.29585719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0303505) q[1];
sx q[1];
rz(-2.6551464) q[1];
sx q[1];
rz(-2.4187947) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2449805) q[3];
sx q[3];
rz(-1.3462726) q[3];
sx q[3];
rz(-2.6905439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.837073) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(1.0860898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(-1.6745463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(0.011209839) q[2];
sx q[2];
rz(-1.8335473) q[2];
sx q[2];
rz(2.0469472) q[2];
rz(0.23444093) q[3];
sx q[3];
rz(-0.81080484) q[3];
sx q[3];
rz(-0.57639359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
