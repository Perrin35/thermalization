OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4797526) q[0];
sx q[0];
rz(-2.2979484) q[0];
sx q[0];
rz(2.9736829) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(-2.8462703) q[1];
sx q[1];
rz(0.056161031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9553298) q[0];
sx q[0];
rz(-1.7832527) q[0];
sx q[0];
rz(2.0583378) q[0];
rz(-pi) q[1];
rz(1.8499608) q[2];
sx q[2];
rz(-0.66510519) q[2];
sx q[2];
rz(-1.6601738) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1311156) q[1];
sx q[1];
rz(-1.3033723) q[1];
sx q[1];
rz(2.1285776) q[1];
x q[2];
rz(-1.6151186) q[3];
sx q[3];
rz(-1.3111776) q[3];
sx q[3];
rz(2.0093902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37796676) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(2.7089233) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(-0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52779657) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(-1.312785) q[0];
rz(-0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.9899433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1412927) q[0];
sx q[0];
rz(-0.89501689) q[0];
sx q[0];
rz(0.38210259) q[0];
rz(-pi) q[1];
rz(1.0863016) q[2];
sx q[2];
rz(-1.0420348) q[2];
sx q[2];
rz(-0.27258401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.093106) q[1];
sx q[1];
rz(-1.152532) q[1];
sx q[1];
rz(0.47936819) q[1];
x q[2];
rz(-1.1851951) q[3];
sx q[3];
rz(-1.10023) q[3];
sx q[3];
rz(3.1009931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(-0.22182626) q[2];
rz(-2.7644073) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-0.036852766) q[0];
rz(2.316078) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(0.056578606) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3765592) q[0];
sx q[0];
rz(-0.74900904) q[0];
sx q[0];
rz(0.88699938) q[0];
rz(-pi) q[1];
rz(0.47358863) q[2];
sx q[2];
rz(-0.28652546) q[2];
sx q[2];
rz(0.63602704) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.024836) q[1];
sx q[1];
rz(-2.7783238) q[1];
sx q[1];
rz(-2.5519752) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9583086) q[3];
sx q[3];
rz(-2.1790677) q[3];
sx q[3];
rz(-2.0351978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-0.92612129) q[2];
rz(2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91519231) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-2.3994989) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-0.4793872) q[1];
sx q[1];
rz(0.46359584) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2984943) q[0];
sx q[0];
rz(-1.8659741) q[0];
sx q[0];
rz(-0.85423268) q[0];
rz(-pi) q[1];
rz(1.6323339) q[2];
sx q[2];
rz(-1.2386285) q[2];
sx q[2];
rz(-2.5663944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0193034) q[1];
sx q[1];
rz(-0.79320723) q[1];
sx q[1];
rz(1.3293468) q[1];
x q[2];
rz(-0.47834088) q[3];
sx q[3];
rz(-1.0676749) q[3];
sx q[3];
rz(-0.92418811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(-0.049499361) q[2];
rz(3.0130623) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(-0.4822576) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(-2.1972426) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0194861) q[0];
sx q[0];
rz(-1.6158316) q[0];
sx q[0];
rz(1.3615863) q[0];
rz(-pi) q[1];
rz(-0.40446754) q[2];
sx q[2];
rz(-0.81768113) q[2];
sx q[2];
rz(-0.58194619) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.797232) q[1];
sx q[1];
rz(-1.5096944) q[1];
sx q[1];
rz(0.0021211591) q[1];
rz(-pi) q[2];
rz(0.43290187) q[3];
sx q[3];
rz(-1.3372984) q[3];
sx q[3];
rz(-1.7683065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(2.8172857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43679431) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-0.58445245) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(-3.086673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4024825) q[0];
sx q[0];
rz(-1.4882898) q[0];
sx q[0];
rz(-1.8399747) q[0];
rz(-pi) q[1];
rz(-1.6266277) q[2];
sx q[2];
rz(-0.32368127) q[2];
sx q[2];
rz(1.8813546) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7683148) q[1];
sx q[1];
rz(-1.0265961) q[1];
sx q[1];
rz(-1.1276223) q[1];
rz(-pi) q[2];
rz(0.46338007) q[3];
sx q[3];
rz(-2.1043092) q[3];
sx q[3];
rz(-2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3871258) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(2.1248655) q[2];
rz(-0.54404849) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-2.8570535) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(-0.91032666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0383354) q[0];
sx q[0];
rz(-0.084780134) q[0];
sx q[0];
rz(-2.2708562) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1548642) q[2];
sx q[2];
rz(-1.2565194) q[2];
sx q[2];
rz(-1.2736125) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6701723) q[1];
sx q[1];
rz(-0.52392611) q[1];
sx q[1];
rz(-2.7736204) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1479285) q[3];
sx q[3];
rz(-2.2054407) q[3];
sx q[3];
rz(-1.5821379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(2.2195623) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(0.63240504) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(-0.30050373) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6209517) q[0];
sx q[0];
rz(-1.2004939) q[0];
sx q[0];
rz(-2.305549) q[0];
x q[1];
rz(2.5007162) q[2];
sx q[2];
rz(-0.87101988) q[2];
sx q[2];
rz(-1.7388294) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.955519) q[1];
sx q[1];
rz(-1.3961853) q[1];
sx q[1];
rz(0.34592918) q[1];
rz(-pi) q[2];
rz(-0.31605966) q[3];
sx q[3];
rz(-2.3028767) q[3];
sx q[3];
rz(2.4162606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5552716) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(2.590495) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(-0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-0.12776275) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(-0.75884563) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1369143) q[0];
sx q[0];
rz(-1.5877962) q[0];
sx q[0];
rz(-1.1466115) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7636289) q[2];
sx q[2];
rz(-2.170993) q[2];
sx q[2];
rz(2.6892975) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2488333) q[1];
sx q[1];
rz(-1.7264688) q[1];
sx q[1];
rz(0.70181904) q[1];
rz(-pi) q[2];
rz(3.127029) q[3];
sx q[3];
rz(-2.2363538) q[3];
sx q[3];
rz(-1.927782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0163429) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(-2.6515567) q[2];
rz(-1.7193433) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(-1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(0.075335659) q[0];
rz(2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-2.5316701) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88203726) q[0];
sx q[0];
rz(-1.2872739) q[0];
sx q[0];
rz(-2.3638704) q[0];
rz(-2.7650325) q[2];
sx q[2];
rz(-1.8278404) q[2];
sx q[2];
rz(0.10274796) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6590609) q[1];
sx q[1];
rz(-2.3606803) q[1];
sx q[1];
rz(0.34300967) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0697332) q[3];
sx q[3];
rz(-0.89322972) q[3];
sx q[3];
rz(0.60615218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23218368) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(2.4278736) q[2];
rz(-2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(1.7238293) q[2];
sx q[2];
rz(-0.19822181) q[2];
sx q[2];
rz(-0.76186686) q[2];
rz(2.8888632) q[3];
sx q[3];
rz(-2.0287632) q[3];
sx q[3];
rz(-0.91026929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
