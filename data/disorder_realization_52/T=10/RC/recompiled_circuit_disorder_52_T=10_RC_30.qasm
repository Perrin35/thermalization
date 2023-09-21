OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5157226) q[0];
sx q[0];
rz(-0.54870257) q[0];
sx q[0];
rz(-0.8843511) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3724369) q[0];
sx q[0];
rz(-2.7363442) q[0];
sx q[0];
rz(2.0408819) q[0];
rz(-2.5764478) q[2];
sx q[2];
rz(-1.5663212) q[2];
sx q[2];
rz(2.7067513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86815392) q[1];
sx q[1];
rz(-1.3248797) q[1];
sx q[1];
rz(1.5032561) q[1];
x q[2];
rz(-0.35688551) q[3];
sx q[3];
rz(-0.89324739) q[3];
sx q[3];
rz(-0.68912904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3551336) q[2];
sx q[2];
rz(-2.3278475) q[2];
sx q[2];
rz(-2.4856429) q[2];
rz(-1.9338699) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(-2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.8939963) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(0.43352747) q[0];
rz(2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(-0.0072335009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4097071) q[0];
sx q[0];
rz(-1.1792372) q[0];
sx q[0];
rz(-3.0407228) q[0];
rz(0.8231926) q[2];
sx q[2];
rz(-0.41831145) q[2];
sx q[2];
rz(-0.38564607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0011598) q[1];
sx q[1];
rz(-1.3507441) q[1];
sx q[1];
rz(0.37186719) q[1];
x q[2];
rz(-2.238027) q[3];
sx q[3];
rz(-1.1840608) q[3];
sx q[3];
rz(-0.29153338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(-2.4439404) q[2];
rz(-3.0200322) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(0.30383032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4678629) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(1.7720222) q[0];
rz(1.9000152) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(-2.8799768) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81944377) q[0];
sx q[0];
rz(-1.561164) q[0];
sx q[0];
rz(1.5528414) q[0];
x q[1];
rz(0.82654731) q[2];
sx q[2];
rz(-1.6625704) q[2];
sx q[2];
rz(-1.2395791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97754543) q[1];
sx q[1];
rz(-2.3512212) q[1];
sx q[1];
rz(0.17069823) q[1];
rz(1.5393799) q[3];
sx q[3];
rz(-2.0835702) q[3];
sx q[3];
rz(-0.98785066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6802784) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(2.2198548) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(-2.8620499) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(-1.571636) q[0];
rz(-2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(-1.8932231) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6682537) q[0];
sx q[0];
rz(-0.11675294) q[0];
sx q[0];
rz(-2.1658685) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20935697) q[2];
sx q[2];
rz(-1.3931837) q[2];
sx q[2];
rz(-0.83724411) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8716988) q[1];
sx q[1];
rz(-1.6819685) q[1];
sx q[1];
rz(-2.486869) q[1];
rz(-pi) q[2];
rz(-0.29420935) q[3];
sx q[3];
rz(-2.3964786) q[3];
sx q[3];
rz(-1.7780768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0288329) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(0.83703414) q[2];
rz(-1.933243) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(-2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(-2.2391589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8474903) q[0];
sx q[0];
rz(-0.93203629) q[0];
sx q[0];
rz(1.627648) q[0];
x q[1];
rz(-1.8158185) q[2];
sx q[2];
rz(-1.1282215) q[2];
sx q[2];
rz(2.4545836) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.93859766) q[1];
sx q[1];
rz(-2.1122167) q[1];
sx q[1];
rz(2.2375537) q[1];
x q[2];
rz(1.2540713) q[3];
sx q[3];
rz(-1.5734908) q[3];
sx q[3];
rz(0.56841422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19501413) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(1.9449332) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(1.4590013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3301795) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(-0.50317558) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(-2.9398289) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87529463) q[0];
sx q[0];
rz(-1.6365956) q[0];
sx q[0];
rz(-1.468303) q[0];
rz(-2.0189507) q[2];
sx q[2];
rz(-1.1345703) q[2];
sx q[2];
rz(-3.089038) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4469874) q[1];
sx q[1];
rz(-1.1519377) q[1];
sx q[1];
rz(-2.6161731) q[1];
x q[2];
rz(-0.73268907) q[3];
sx q[3];
rz(-1.9273888) q[3];
sx q[3];
rz(2.4458812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21489828) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(0.26724896) q[2];
rz(-0.8231419) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(2.2657623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-3.0840432) q[0];
rz(-1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(0.94271359) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0038303) q[0];
sx q[0];
rz(-0.98235213) q[0];
sx q[0];
rz(1.0914735) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1812181) q[2];
sx q[2];
rz(-2.865961) q[2];
sx q[2];
rz(0.20197091) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49589866) q[1];
sx q[1];
rz(-1.2518479) q[1];
sx q[1];
rz(-2.0046528) q[1];
rz(0.25992486) q[3];
sx q[3];
rz(-2.0689031) q[3];
sx q[3];
rz(-0.92344027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1192347) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(-2.7589202) q[2];
rz(2.102397) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(0.97810811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.7169749) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(0.27012816) q[0];
rz(0.62942901) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(-2.8576635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2487508) q[0];
sx q[0];
rz(-2.1584956) q[0];
sx q[0];
rz(-2.6177004) q[0];
rz(-0.59992744) q[2];
sx q[2];
rz(-1.4495965) q[2];
sx q[2];
rz(1.8184513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.200951) q[1];
sx q[1];
rz(-2.0307699) q[1];
sx q[1];
rz(2.5775787) q[1];
rz(-pi) q[2];
rz(1.5408526) q[3];
sx q[3];
rz(-2.2303914) q[3];
sx q[3];
rz(-2.3823882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5876864) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(-1.2109057) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0652086) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(0.22924766) q[0];
rz(0.30300888) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(-1.4607666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0544573) q[0];
sx q[0];
rz(-1.5532877) q[0];
sx q[0];
rz(-2.8580335) q[0];
rz(-pi) q[1];
rz(0.42713366) q[2];
sx q[2];
rz(-2.7286227) q[2];
sx q[2];
rz(0.6461179) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8252392) q[1];
sx q[1];
rz(-1.5484527) q[1];
sx q[1];
rz(-0.52492001) q[1];
rz(-pi) q[2];
rz(-2.3567696) q[3];
sx q[3];
rz(-1.8514957) q[3];
sx q[3];
rz(2.7626038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71172697) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(-0.43186489) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(-0.7152043) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(-0.33690548) q[0];
rz(-2.9341872) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(-2.7609603) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779553) q[0];
sx q[0];
rz(-1.5008238) q[0];
sx q[0];
rz(-1.3689343) q[0];
rz(-pi) q[1];
rz(-2.400488) q[2];
sx q[2];
rz(-2.0217102) q[2];
sx q[2];
rz(-1.2620743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3818647) q[1];
sx q[1];
rz(-1.7654395) q[1];
sx q[1];
rz(1.0641644) q[1];
rz(-1.9862513) q[3];
sx q[3];
rz(-1.0329773) q[3];
sx q[3];
rz(-0.95738639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1404861) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(2.005119) q[2];
rz(0.040955695) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(-1.3142746) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8052335) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(-2.1144755) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(0.2480416) q[2];
sx q[2];
rz(-1.3477865) q[2];
sx q[2];
rz(-1.8852521) q[2];
rz(2.0102262) q[3];
sx q[3];
rz(-1.3888748) q[3];
sx q[3];
rz(0.5007762) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];