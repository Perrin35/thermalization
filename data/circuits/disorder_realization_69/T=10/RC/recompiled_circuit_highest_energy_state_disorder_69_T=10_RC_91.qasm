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
rz(0.82290736) q[0];
sx q[0];
rz(-0.35879254) q[0];
sx q[0];
rz(0.86831492) q[0];
rz(-1.4153642) q[1];
sx q[1];
rz(-1.1338898) q[1];
sx q[1];
rz(-0.99376065) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.438765) q[0];
sx q[0];
rz(-1.1701705) q[0];
sx q[0];
rz(1.6933269) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6401071) q[2];
sx q[2];
rz(-1.9619313) q[2];
sx q[2];
rz(-2.8361965) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3512416) q[1];
sx q[1];
rz(-1.1415485) q[1];
sx q[1];
rz(1.0830948) q[1];
rz(-2.0076143) q[3];
sx q[3];
rz(-1.7941495) q[3];
sx q[3];
rz(-0.88296605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.913784) q[2];
sx q[2];
rz(-1.3334393) q[2];
sx q[2];
rz(2.559973) q[2];
rz(-2.4070814) q[3];
sx q[3];
rz(-1.4903277) q[3];
sx q[3];
rz(0.41828004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89473474) q[0];
sx q[0];
rz(-1.3792091) q[0];
sx q[0];
rz(-2.6439164) q[0];
rz(-2.1007762) q[1];
sx q[1];
rz(-2.7383995) q[1];
sx q[1];
rz(1.2373479) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0859511) q[0];
sx q[0];
rz(-0.72362542) q[0];
sx q[0];
rz(3.0776204) q[0];
rz(-2.9534229) q[2];
sx q[2];
rz(-1.9261179) q[2];
sx q[2];
rz(-0.73588348) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50857568) q[1];
sx q[1];
rz(-2.687722) q[1];
sx q[1];
rz(1.244844) q[1];
rz(-pi) q[2];
rz(-2.4604843) q[3];
sx q[3];
rz(-0.33675413) q[3];
sx q[3];
rz(2.3189173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.70025468) q[2];
sx q[2];
rz(-2.4105218) q[2];
sx q[2];
rz(0.31044427) q[2];
rz(-1.9849518) q[3];
sx q[3];
rz(-1.1247331) q[3];
sx q[3];
rz(1.9748851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0079086) q[0];
sx q[0];
rz(-0.5846566) q[0];
sx q[0];
rz(0.34580082) q[0];
rz(2.8969104) q[1];
sx q[1];
rz(-0.9318277) q[1];
sx q[1];
rz(-1.7049047) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0416243) q[0];
sx q[0];
rz(-1.1505373) q[0];
sx q[0];
rz(2.7970275) q[0];
rz(-2.5710996) q[2];
sx q[2];
rz(-1.2353503) q[2];
sx q[2];
rz(-2.6550309) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2785221) q[1];
sx q[1];
rz(-0.82048847) q[1];
sx q[1];
rz(0.29149518) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9448024) q[3];
sx q[3];
rz(-1.1218277) q[3];
sx q[3];
rz(1.6351007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1338542) q[2];
sx q[2];
rz(-2.014092) q[2];
sx q[2];
rz(0.63068843) q[2];
rz(2.808029) q[3];
sx q[3];
rz(-2.1400698) q[3];
sx q[3];
rz(2.1400616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0665322) q[0];
sx q[0];
rz(-2.1933031) q[0];
sx q[0];
rz(1.4045658) q[0];
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
rz(2.7333711) q[0];
sx q[0];
rz(-1.2722172) q[0];
sx q[0];
rz(2.7135506) q[0];
x q[1];
rz(-2.4253885) q[2];
sx q[2];
rz(-2.5599481) q[2];
sx q[2];
rz(0.003613506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4300116) q[1];
sx q[1];
rz(-2.0329355) q[1];
sx q[1];
rz(-0.26600809) q[1];
x q[2];
rz(-2.4749807) q[3];
sx q[3];
rz(-0.48891196) q[3];
sx q[3];
rz(2.0036445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.839445) q[2];
sx q[2];
rz(-1.2371233) q[2];
sx q[2];
rz(-0.5298003) q[2];
rz(3.1032622) q[3];
sx q[3];
rz(-0.72961346) q[3];
sx q[3];
rz(-1.5420325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2235276) q[0];
sx q[0];
rz(-0.46850884) q[0];
sx q[0];
rz(1.9388306) q[0];
rz(-0.27944061) q[1];
sx q[1];
rz(-2.1165106) q[1];
sx q[1];
rz(-0.7739982) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27085486) q[0];
sx q[0];
rz(-1.2574728) q[0];
sx q[0];
rz(-2.3148651) q[0];
rz(2.9086529) q[2];
sx q[2];
rz(-1.4904984) q[2];
sx q[2];
rz(-0.38903076) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2629649) q[1];
sx q[1];
rz(-1.5836746) q[1];
sx q[1];
rz(1.7098544) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5755025) q[3];
sx q[3];
rz(-0.53541267) q[3];
sx q[3];
rz(1.104606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1191795) q[2];
sx q[2];
rz(-0.1410307) q[2];
sx q[2];
rz(2.5775583) q[2];
rz(2.4885079) q[3];
sx q[3];
rz(-1.1354732) q[3];
sx q[3];
rz(2.691332) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7804724) q[0];
sx q[0];
rz(-0.55279624) q[0];
sx q[0];
rz(-2.4984388) q[0];
rz(1.9505352) q[1];
sx q[1];
rz(-1.6845614) q[1];
sx q[1];
rz(2.4868884) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7793286) q[0];
sx q[0];
rz(-1.6342499) q[0];
sx q[0];
rz(-0.16880798) q[0];
rz(-1.5327318) q[2];
sx q[2];
rz(-0.88721472) q[2];
sx q[2];
rz(-2.3823007) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7186942) q[1];
sx q[1];
rz(-2.5942583) q[1];
sx q[1];
rz(1.4446299) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.07043802) q[3];
sx q[3];
rz(-1.1813191) q[3];
sx q[3];
rz(0.73274437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6341298) q[2];
sx q[2];
rz(-2.7644988) q[2];
sx q[2];
rz(-1.6667574) q[2];
rz(-0.48464388) q[3];
sx q[3];
rz(-2.1438997) q[3];
sx q[3];
rz(-1.8528329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5019048) q[0];
sx q[0];
rz(-2.8785093) q[0];
sx q[0];
rz(-0.68341533) q[0];
rz(3.0112093) q[1];
sx q[1];
rz(-1.5905292) q[1];
sx q[1];
rz(0.032616671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0158247) q[0];
sx q[0];
rz(-0.59249632) q[0];
sx q[0];
rz(-1.2529208) q[0];
rz(-pi) q[1];
rz(-0.68565418) q[2];
sx q[2];
rz(-1.4156121) q[2];
sx q[2];
rz(-2.6557166) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8396388) q[1];
sx q[1];
rz(-0.31529135) q[1];
sx q[1];
rz(-2.1350506) q[1];
rz(-pi) q[2];
rz(2.8835589) q[3];
sx q[3];
rz(-0.88538187) q[3];
sx q[3];
rz(3.0436707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7183097) q[2];
sx q[2];
rz(-2.0330567) q[2];
sx q[2];
rz(-0.56524593) q[2];
rz(1.3609173) q[3];
sx q[3];
rz(-2.8748685) q[3];
sx q[3];
rz(-2.3274073) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3703506) q[0];
sx q[0];
rz(-0.057436198) q[0];
sx q[0];
rz(-0.29956079) q[0];
rz(-1.6948505) q[1];
sx q[1];
rz(-2.2694777) q[1];
sx q[1];
rz(2.1902693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46630105) q[0];
sx q[0];
rz(-1.5383175) q[0];
sx q[0];
rz(-1.3656653) q[0];
x q[1];
rz(-0.0092333992) q[2];
sx q[2];
rz(-1.3592741) q[2];
sx q[2];
rz(0.68995014) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.42859498) q[1];
sx q[1];
rz(-2.7168437) q[1];
sx q[1];
rz(-3.0921169) q[1];
rz(-pi) q[2];
rz(-1.8883315) q[3];
sx q[3];
rz(-2.866689) q[3];
sx q[3];
rz(1.8567228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99913725) q[2];
sx q[2];
rz(-2.3951525) q[2];
sx q[2];
rz(-2.8738521) q[2];
rz(1.0033876) q[3];
sx q[3];
rz(-2.181874) q[3];
sx q[3];
rz(0.80823922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35850152) q[0];
sx q[0];
rz(-2.6434904) q[0];
sx q[0];
rz(2.1844693) q[0];
rz(0.3793017) q[1];
sx q[1];
rz(-0.4117659) q[1];
sx q[1];
rz(2.9764825) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6388164) q[0];
sx q[0];
rz(-0.90327016) q[0];
sx q[0];
rz(2.0124042) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0909507) q[2];
sx q[2];
rz(-1.0989185) q[2];
sx q[2];
rz(2.7993921) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2735426) q[1];
sx q[1];
rz(-0.99813491) q[1];
sx q[1];
rz(-2.9258201) q[1];
rz(-2.1741207) q[3];
sx q[3];
rz(-1.3730197) q[3];
sx q[3];
rz(-2.8654049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86753201) q[2];
sx q[2];
rz(-1.910285) q[2];
sx q[2];
rz(2.3629698) q[2];
rz(2.8042931) q[3];
sx q[3];
rz(-1.0236579) q[3];
sx q[3];
rz(1.5571099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.857665) q[0];
sx q[0];
rz(-2.0511257) q[0];
sx q[0];
rz(-1.0887867) q[0];
rz(-0.42824832) q[1];
sx q[1];
rz(-2.1183522) q[1];
sx q[1];
rz(-0.64839378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3449865) q[0];
sx q[0];
rz(-1.8341176) q[0];
sx q[0];
rz(-0.4395606) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6971385) q[2];
sx q[2];
rz(-2.8011232) q[2];
sx q[2];
rz(2.7012417) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9251483) q[1];
sx q[1];
rz(-2.2183462) q[1];
sx q[1];
rz(1.5179894) q[1];
rz(-pi) q[2];
rz(-0.56364735) q[3];
sx q[3];
rz(-1.6674588) q[3];
sx q[3];
rz(-0.046162995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8757561) q[2];
sx q[2];
rz(-2.0529604) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.6801878) q[0];
sx q[0];
rz(-2.9869933) q[0];
sx q[0];
rz(-2.3416478) q[0];
rz(-1.9752621) q[1];
sx q[1];
rz(-2.1569398) q[1];
sx q[1];
rz(0.91469761) q[1];
rz(2.2565319) q[2];
sx q[2];
rz(-2.764438) q[2];
sx q[2];
rz(0.59791313) q[2];
rz(-0.33954444) q[3];
sx q[3];
rz(-1.3673269) q[3];
sx q[3];
rz(0.23585503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
