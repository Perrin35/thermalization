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
rz(-1.4827363) q[0];
sx q[0];
rz(-2.1602477) q[0];
sx q[0];
rz(2.0438097) q[0];
rz(-0.78805796) q[1];
sx q[1];
rz(-2.0907953) q[1];
sx q[1];
rz(-0.0069590574) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58034694) q[0];
sx q[0];
rz(-2.2202507) q[0];
sx q[0];
rz(-3.0710996) q[0];
rz(-pi) q[1];
rz(0.9240146) q[2];
sx q[2];
rz(-1.7597464) q[2];
sx q[2];
rz(-1.3786157) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4930058) q[1];
sx q[1];
rz(-1.8064335) q[1];
sx q[1];
rz(-0.034502397) q[1];
rz(-pi) q[2];
rz(2.801729) q[3];
sx q[3];
rz(-1.6714665) q[3];
sx q[3];
rz(-1.0235746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9030582) q[2];
sx q[2];
rz(-1.3471341) q[2];
sx q[2];
rz(-1.4702338) q[2];
rz(3.0155731) q[3];
sx q[3];
rz(-2.6991548) q[3];
sx q[3];
rz(2.457705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0266492) q[0];
sx q[0];
rz(-0.4489972) q[0];
sx q[0];
rz(0.63585109) q[0];
rz(0.77330971) q[1];
sx q[1];
rz(-2.4644303) q[1];
sx q[1];
rz(-1.4453452) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34958865) q[0];
sx q[0];
rz(-0.82539648) q[0];
sx q[0];
rz(2.1950153) q[0];
x q[1];
rz(2.0436297) q[2];
sx q[2];
rz(-0.84542984) q[2];
sx q[2];
rz(0.54976094) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.066582908) q[1];
sx q[1];
rz(-0.72117469) q[1];
sx q[1];
rz(0.23951342) q[1];
rz(-pi) q[2];
rz(1.0926756) q[3];
sx q[3];
rz(-2.993686) q[3];
sx q[3];
rz(0.92648244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.14629743) q[2];
sx q[2];
rz(-1.7701021) q[2];
sx q[2];
rz(-0.13776097) q[2];
rz(2.6855101) q[3];
sx q[3];
rz(-2.5465953) q[3];
sx q[3];
rz(-2.6661787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.544203) q[0];
sx q[0];
rz(-1.5837639) q[0];
sx q[0];
rz(0.57918817) q[0];
rz(3.0384565) q[1];
sx q[1];
rz(-1.8201647) q[1];
sx q[1];
rz(0.78027049) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0294757) q[0];
sx q[0];
rz(-1.5896086) q[0];
sx q[0];
rz(0.095186724) q[0];
rz(-0.19635503) q[2];
sx q[2];
rz(-2.4508424) q[2];
sx q[2];
rz(-0.029411246) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0611813) q[1];
sx q[1];
rz(-1.965251) q[1];
sx q[1];
rz(2.9030672) q[1];
x q[2];
rz(1.8336201) q[3];
sx q[3];
rz(-0.71280957) q[3];
sx q[3];
rz(-1.3808911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66037336) q[2];
sx q[2];
rz(-1.4619091) q[2];
sx q[2];
rz(0.090864651) q[2];
rz(1.6615435) q[3];
sx q[3];
rz(-1.816498) q[3];
sx q[3];
rz(0.81502325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12260967) q[0];
sx q[0];
rz(-2.9538587) q[0];
sx q[0];
rz(-2.6222099) q[0];
rz(1.9620365) q[1];
sx q[1];
rz(-2.4603381) q[1];
sx q[1];
rz(-3.093241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7413986) q[0];
sx q[0];
rz(-2.1693008) q[0];
sx q[0];
rz(0.6299751) q[0];
x q[1];
rz(-0.3860571) q[2];
sx q[2];
rz(-2.0252844) q[2];
sx q[2];
rz(0.41815652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0569339) q[1];
sx q[1];
rz(-0.48843436) q[1];
sx q[1];
rz(3.0822251) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3062258) q[3];
sx q[3];
rz(-0.79341429) q[3];
sx q[3];
rz(-2.1891264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4590596) q[2];
sx q[2];
rz(-0.30632633) q[2];
sx q[2];
rz(-1.4502067) q[2];
rz(-1.1059443) q[3];
sx q[3];
rz(-1.148843) q[3];
sx q[3];
rz(-0.86627427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76097101) q[0];
sx q[0];
rz(-0.80046099) q[0];
sx q[0];
rz(2.3642484) q[0];
rz(-2.6457973) q[1];
sx q[1];
rz(-2.7340041) q[1];
sx q[1];
rz(2.9487603) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7161761) q[0];
sx q[0];
rz(-1.2649853) q[0];
sx q[0];
rz(-2.5772463) q[0];
rz(-pi) q[1];
rz(2.3740511) q[2];
sx q[2];
rz(-1.3965675) q[2];
sx q[2];
rz(1.7310639) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7954171) q[1];
sx q[1];
rz(-2.8320791) q[1];
sx q[1];
rz(1.2796106) q[1];
rz(-1.9726095) q[3];
sx q[3];
rz(-1.4207559) q[3];
sx q[3];
rz(0.91317859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5567646) q[2];
sx q[2];
rz(-0.77797055) q[2];
sx q[2];
rz(0.83621109) q[2];
rz(-0.95544514) q[3];
sx q[3];
rz(-1.7183869) q[3];
sx q[3];
rz(0.53171617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29109508) q[0];
sx q[0];
rz(-1.2245155) q[0];
sx q[0];
rz(-3.1078872) q[0];
rz(-2.2233502) q[1];
sx q[1];
rz(-0.70437175) q[1];
sx q[1];
rz(-2.4423626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3501773) q[0];
sx q[0];
rz(-2.3473513) q[0];
sx q[0];
rz(0.82483952) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9543058) q[2];
sx q[2];
rz(-1.6237055) q[2];
sx q[2];
rz(-2.1250181) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4650326) q[1];
sx q[1];
rz(-1.1313712) q[1];
sx q[1];
rz(-1.7656754) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27856234) q[3];
sx q[3];
rz(-1.8980366) q[3];
sx q[3];
rz(-0.95765169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0596727) q[2];
sx q[2];
rz(-0.6531738) q[2];
sx q[2];
rz(-3.0461404) q[2];
rz(-2.0898315) q[3];
sx q[3];
rz(-1.2442518) q[3];
sx q[3];
rz(2.2473647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28602257) q[0];
sx q[0];
rz(-0.48208553) q[0];
sx q[0];
rz(1.1676189) q[0];
rz(-2.7883912) q[1];
sx q[1];
rz(-0.51274061) q[1];
sx q[1];
rz(-2.8599427) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5198811) q[0];
sx q[0];
rz(-2.6586652) q[0];
sx q[0];
rz(1.7596702) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4942964) q[2];
sx q[2];
rz(-1.4855097) q[2];
sx q[2];
rz(-2.6814987) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0822009) q[1];
sx q[1];
rz(-1.7052461) q[1];
sx q[1];
rz(-1.7198635) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1280888) q[3];
sx q[3];
rz(-1.1508599) q[3];
sx q[3];
rz(2.1518663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1436651) q[2];
sx q[2];
rz(-1.9196271) q[2];
sx q[2];
rz(1.8188875) q[2];
rz(-1.5023242) q[3];
sx q[3];
rz(-1.084525) q[3];
sx q[3];
rz(0.098202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99763501) q[0];
sx q[0];
rz(-1.8956381) q[0];
sx q[0];
rz(-0.27398807) q[0];
rz(2.5471089) q[1];
sx q[1];
rz(-1.7285873) q[1];
sx q[1];
rz(-2.6110173) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45121058) q[0];
sx q[0];
rz(-1.57461) q[0];
sx q[0];
rz(-0.03937748) q[0];
rz(-pi) q[1];
rz(0.12315788) q[2];
sx q[2];
rz(-2.5854857) q[2];
sx q[2];
rz(-1.7983939) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.952632) q[1];
sx q[1];
rz(-2.7285693) q[1];
sx q[1];
rz(2.7160591) q[1];
x q[2];
rz(-0.16333111) q[3];
sx q[3];
rz(-0.86465752) q[3];
sx q[3];
rz(2.284019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8990367) q[2];
sx q[2];
rz(-1.2071004) q[2];
sx q[2];
rz(-2.8846018) q[2];
rz(1.1993923) q[3];
sx q[3];
rz(-1.896984) q[3];
sx q[3];
rz(-1.6283584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5523819) q[0];
sx q[0];
rz(-0.95443812) q[0];
sx q[0];
rz(-2.6115665) q[0];
rz(-1.6389305) q[1];
sx q[1];
rz(-1.1549779) q[1];
sx q[1];
rz(0.93856215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28806557) q[0];
sx q[0];
rz(-2.1479283) q[0];
sx q[0];
rz(2.1398628) q[0];
rz(-2.4468719) q[2];
sx q[2];
rz(-2.955547) q[2];
sx q[2];
rz(-3.0228066) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1129866) q[1];
sx q[1];
rz(-2.6588106) q[1];
sx q[1];
rz(-2.042599) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7750562) q[3];
sx q[3];
rz(-2.5464749) q[3];
sx q[3];
rz(2.4581153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3878801) q[2];
sx q[2];
rz(-2.2215999) q[2];
sx q[2];
rz(-2.06125) q[2];
rz(-0.8148109) q[3];
sx q[3];
rz(-1.1956513) q[3];
sx q[3];
rz(1.4174392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1431047) q[0];
sx q[0];
rz(-1.657722) q[0];
sx q[0];
rz(-1.0888354) q[0];
rz(-0.067527436) q[1];
sx q[1];
rz(-1.2779526) q[1];
sx q[1];
rz(-2.3823104) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7831206) q[0];
sx q[0];
rz(-1.6896588) q[0];
sx q[0];
rz(1.3141743) q[0];
rz(-pi) q[1];
rz(1.6948872) q[2];
sx q[2];
rz(-1.3465836) q[2];
sx q[2];
rz(-1.2808764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85222178) q[1];
sx q[1];
rz(-1.3557649) q[1];
sx q[1];
rz(-0.76070666) q[1];
rz(2.5645026) q[3];
sx q[3];
rz(-2.7365757) q[3];
sx q[3];
rz(1.7187723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8074983) q[2];
sx q[2];
rz(-2.0457902) q[2];
sx q[2];
rz(-0.57541543) q[2];
rz(0.56257644) q[3];
sx q[3];
rz(-0.96499363) q[3];
sx q[3];
rz(1.3728728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(0.060487735) q[0];
sx q[0];
rz(-1.2860379) q[0];
sx q[0];
rz(2.2338569) q[0];
rz(2.0545215) q[1];
sx q[1];
rz(-2.0094951) q[1];
sx q[1];
rz(3.1241945) q[1];
rz(0.74093735) q[2];
sx q[2];
rz(-2.9334684) q[2];
sx q[2];
rz(2.8192782) q[2];
rz(-0.59845944) q[3];
sx q[3];
rz(-1.2597407) q[3];
sx q[3];
rz(1.5516439) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
