OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.1640373) q[0];
sx q[0];
rz(-1.9174175) q[0];
sx q[0];
rz(1.8608215) q[0];
rz(-0.7723074) q[1];
sx q[1];
rz(-0.63900715) q[1];
sx q[1];
rz(0.26689902) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.229202) q[0];
sx q[0];
rz(-1.4773737) q[0];
sx q[0];
rz(1.4742729) q[0];
rz(2.884567) q[2];
sx q[2];
rz(-0.89867175) q[2];
sx q[2];
rz(0.66115236) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.34632979) q[1];
sx q[1];
rz(-2.2861233) q[1];
sx q[1];
rz(-2.1500509) q[1];
rz(2.6990864) q[3];
sx q[3];
rz(-1.6638262) q[3];
sx q[3];
rz(-2.0131567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.018517) q[2];
sx q[2];
rz(-0.96394959) q[2];
sx q[2];
rz(-2.4289971) q[2];
rz(2.9739042) q[3];
sx q[3];
rz(-0.84961397) q[3];
sx q[3];
rz(0.4445506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1009104) q[0];
sx q[0];
rz(-1.3586783) q[0];
sx q[0];
rz(-1.4605301) q[0];
rz(-1.4617317) q[1];
sx q[1];
rz(-1.0667543) q[1];
sx q[1];
rz(1.1163968) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3718263) q[0];
sx q[0];
rz(-0.049749181) q[0];
sx q[0];
rz(-2.0905963) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1266534) q[2];
sx q[2];
rz(-1.4420605) q[2];
sx q[2];
rz(-1.5644846) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6608289) q[1];
sx q[1];
rz(-2.6763981) q[1];
sx q[1];
rz(0.61252266) q[1];
rz(-pi) q[2];
rz(-2.7313569) q[3];
sx q[3];
rz(-1.4465152) q[3];
sx q[3];
rz(1.6491485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1560893) q[2];
sx q[2];
rz(-0.49571338) q[2];
sx q[2];
rz(-0.25225857) q[2];
rz(1.0559399) q[3];
sx q[3];
rz(-2.2838433) q[3];
sx q[3];
rz(2.1411538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58241874) q[0];
sx q[0];
rz(-2.2370339) q[0];
sx q[0];
rz(-0.82758033) q[0];
rz(2.6170392) q[1];
sx q[1];
rz(-1.8962212) q[1];
sx q[1];
rz(-0.96447271) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37999718) q[0];
sx q[0];
rz(-0.58172136) q[0];
sx q[0];
rz(-1.6935478) q[0];
rz(-pi) q[1];
rz(-3.0542637) q[2];
sx q[2];
rz(-1.9448408) q[2];
sx q[2];
rz(0.61368521) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2169184) q[1];
sx q[1];
rz(-0.63230941) q[1];
sx q[1];
rz(-3.0576474) q[1];
rz(-pi) q[2];
rz(-1.6286704) q[3];
sx q[3];
rz(-2.0116968) q[3];
sx q[3];
rz(-0.98516516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7437637) q[2];
sx q[2];
rz(-0.22481329) q[2];
sx q[2];
rz(1.942983) q[2];
rz(-1.648929) q[3];
sx q[3];
rz(-0.97380939) q[3];
sx q[3];
rz(-0.57884136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8156133) q[0];
sx q[0];
rz(-0.15376832) q[0];
sx q[0];
rz(-1.2157259) q[0];
rz(2.1513596) q[1];
sx q[1];
rz(-0.74140397) q[1];
sx q[1];
rz(0.56891099) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5241107) q[0];
sx q[0];
rz(-2.4883399) q[0];
sx q[0];
rz(-2.0481443) q[0];
rz(-pi) q[1];
rz(-0.96069889) q[2];
sx q[2];
rz(-0.44126661) q[2];
sx q[2];
rz(0.91741761) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4601656) q[1];
sx q[1];
rz(-0.54124628) q[1];
sx q[1];
rz(-0.98249225) q[1];
x q[2];
rz(0.74796933) q[3];
sx q[3];
rz(-1.0587721) q[3];
sx q[3];
rz(-1.0486163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27900055) q[2];
sx q[2];
rz(-0.96070868) q[2];
sx q[2];
rz(0.29435364) q[2];
rz(0.66458464) q[3];
sx q[3];
rz(-2.3817101) q[3];
sx q[3];
rz(-1.4544646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6890474) q[0];
sx q[0];
rz(-1.713546) q[0];
sx q[0];
rz(-0.31387615) q[0];
rz(2.2398056) q[1];
sx q[1];
rz(-1.3548464) q[1];
sx q[1];
rz(-1.6759759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0756158) q[0];
sx q[0];
rz(-0.1217723) q[0];
sx q[0];
rz(-1.3568702) q[0];
rz(-3.0283233) q[2];
sx q[2];
rz(-1.3217032) q[2];
sx q[2];
rz(1.8367565) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.92815968) q[1];
sx q[1];
rz(-1.4959236) q[1];
sx q[1];
rz(-2.7533745) q[1];
rz(-pi) q[2];
rz(3.0907822) q[3];
sx q[3];
rz(-0.63119315) q[3];
sx q[3];
rz(-2.7630591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9471211) q[2];
sx q[2];
rz(-2.3919545) q[2];
sx q[2];
rz(-2.0571902) q[2];
rz(2.8179152) q[3];
sx q[3];
rz(-2.4249228) q[3];
sx q[3];
rz(-2.9173541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9913469) q[0];
sx q[0];
rz(-2.8960189) q[0];
sx q[0];
rz(2.7299951) q[0];
rz(0.72548524) q[1];
sx q[1];
rz(-1.8975703) q[1];
sx q[1];
rz(-1.4010319) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66564225) q[0];
sx q[0];
rz(-0.67806292) q[0];
sx q[0];
rz(2.8286788) q[0];
rz(-3.0040222) q[2];
sx q[2];
rz(-0.74512945) q[2];
sx q[2];
rz(0.18411769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.42565517) q[1];
sx q[1];
rz(-1.2285985) q[1];
sx q[1];
rz(-1.2340496) q[1];
x q[2];
rz(-2.9855491) q[3];
sx q[3];
rz(-0.45609176) q[3];
sx q[3];
rz(1.9182916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1650042) q[2];
sx q[2];
rz(-0.62825957) q[2];
sx q[2];
rz(1.4264301) q[2];
rz(0.12380883) q[3];
sx q[3];
rz(-2.0721469) q[3];
sx q[3];
rz(-2.6482705) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2774778) q[0];
sx q[0];
rz(-1.9954229) q[0];
sx q[0];
rz(1.3364828) q[0];
rz(-0.9388963) q[1];
sx q[1];
rz(-2.0195596) q[1];
sx q[1];
rz(2.8575361) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2920609) q[0];
sx q[0];
rz(-2.1725328) q[0];
sx q[0];
rz(-1.1547778) q[0];
rz(0.16347537) q[2];
sx q[2];
rz(-1.7942085) q[2];
sx q[2];
rz(2.7976409) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.56834451) q[1];
sx q[1];
rz(-1.580984) q[1];
sx q[1];
rz(-1.6022026) q[1];
rz(-0.18064336) q[3];
sx q[3];
rz(-1.4339281) q[3];
sx q[3];
rz(1.5113854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17860086) q[2];
sx q[2];
rz(-0.68135771) q[2];
sx q[2];
rz(2.2006688) q[2];
rz(0.00087794463) q[3];
sx q[3];
rz(-1.2255171) q[3];
sx q[3];
rz(-0.0860478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0963928) q[0];
sx q[0];
rz(-0.89458507) q[0];
sx q[0];
rz(-0.17564242) q[0];
rz(1.9215709) q[1];
sx q[1];
rz(-0.86985391) q[1];
sx q[1];
rz(-0.87179914) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.652013) q[0];
sx q[0];
rz(-0.30569363) q[0];
sx q[0];
rz(-2.0920805) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9398676) q[2];
sx q[2];
rz(-0.82316527) q[2];
sx q[2];
rz(1.5368423) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.84874547) q[1];
sx q[1];
rz(-2.1715144) q[1];
sx q[1];
rz(-1.7639746) q[1];
rz(-pi) q[2];
rz(1.4967887) q[3];
sx q[3];
rz(-1.0153061) q[3];
sx q[3];
rz(-3.1233226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4153727) q[2];
sx q[2];
rz(-1.8612334) q[2];
sx q[2];
rz(1.0225164) q[2];
rz(-2.7546049) q[3];
sx q[3];
rz(-1.3849247) q[3];
sx q[3];
rz(-2.3302087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33685327) q[0];
sx q[0];
rz(-1.9751208) q[0];
sx q[0];
rz(2.8785896) q[0];
rz(0.31245843) q[1];
sx q[1];
rz(-1.0954906) q[1];
sx q[1];
rz(-1.1824664) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1240886) q[0];
sx q[0];
rz(-2.5944983) q[0];
sx q[0];
rz(-1.1897683) q[0];
rz(-pi) q[1];
rz(-0.83177213) q[2];
sx q[2];
rz(-1.1397994) q[2];
sx q[2];
rz(0.43048358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7484365) q[1];
sx q[1];
rz(-1.1999722) q[1];
sx q[1];
rz(1.3575587) q[1];
rz(-2.1402604) q[3];
sx q[3];
rz(-1.1181485) q[3];
sx q[3];
rz(-2.4695726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.075228127) q[2];
sx q[2];
rz(-1.790739) q[2];
sx q[2];
rz(-2.8933375) q[2];
rz(0.17812854) q[3];
sx q[3];
rz(-1.0823366) q[3];
sx q[3];
rz(0.78290141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.42661509) q[0];
sx q[0];
rz(-2.7474032) q[0];
sx q[0];
rz(1.0618427) q[0];
rz(-1.8408403) q[1];
sx q[1];
rz(-1.5251093) q[1];
sx q[1];
rz(1.1311857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0164364) q[0];
sx q[0];
rz(-1.1068692) q[0];
sx q[0];
rz(0.50986503) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6386191) q[2];
sx q[2];
rz(-1.1736693) q[2];
sx q[2];
rz(0.056588825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85581814) q[1];
sx q[1];
rz(-0.40935358) q[1];
sx q[1];
rz(2.9079283) q[1];
rz(-pi) q[2];
rz(-1.3374109) q[3];
sx q[3];
rz(-2.1595229) q[3];
sx q[3];
rz(1.0860541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0807557) q[2];
sx q[2];
rz(-2.8751825) q[2];
sx q[2];
rz(0.60144919) q[2];
rz(-1.0414177) q[3];
sx q[3];
rz(-1.4789378) q[3];
sx q[3];
rz(1.7634332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9926485) q[0];
sx q[0];
rz(-2.5288378) q[0];
sx q[0];
rz(0.80768325) q[0];
rz(0.07769892) q[1];
sx q[1];
rz(-0.54010375) q[1];
sx q[1];
rz(1.7355951) q[1];
rz(-0.40602691) q[2];
sx q[2];
rz(-2.1320504) q[2];
sx q[2];
rz(-0.8036094) q[2];
rz(-0.79176767) q[3];
sx q[3];
rz(-1.0825915) q[3];
sx q[3];
rz(-1.3655567) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
