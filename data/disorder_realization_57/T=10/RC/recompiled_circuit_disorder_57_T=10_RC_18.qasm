OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6141619) q[0];
sx q[0];
rz(-0.80978137) q[0];
sx q[0];
rz(2.6101987) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(1.3637435) q[1];
sx q[1];
rz(10.627828) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42703585) q[0];
sx q[0];
rz(-1.3792975) q[0];
sx q[0];
rz(1.3546076) q[0];
rz(0.7026303) q[2];
sx q[2];
rz(-2.812817) q[2];
sx q[2];
rz(1.6871014) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10567084) q[1];
sx q[1];
rz(-2.8997313) q[1];
sx q[1];
rz(2.9015404) q[1];
x q[2];
rz(1.5337464) q[3];
sx q[3];
rz(-1.2780446) q[3];
sx q[3];
rz(2.1238413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27665859) q[2];
sx q[2];
rz(-2.0099137) q[2];
sx q[2];
rz(0.19908389) q[2];
rz(1.6137326) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(0.73408192) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1428225) q[0];
sx q[0];
rz(-3.0144189) q[0];
sx q[0];
rz(-0.81940991) q[0];
rz(2.8557414) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(-1.2664638) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74840141) q[0];
sx q[0];
rz(-1.9494434) q[0];
sx q[0];
rz(-0.65403954) q[0];
rz(0.23898791) q[2];
sx q[2];
rz(-1.1766489) q[2];
sx q[2];
rz(0.1220526) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5310865) q[1];
sx q[1];
rz(-0.55263457) q[1];
sx q[1];
rz(-2.3081739) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64139069) q[3];
sx q[3];
rz(-0.68123276) q[3];
sx q[3];
rz(0.89279803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1566029) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(1.161969) q[2];
rz(-3.0544288) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(-3.0676837) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53482985) q[0];
sx q[0];
rz(-2.02089) q[0];
sx q[0];
rz(2.8379922) q[0];
rz(-1.7595694) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(-0.64750013) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9377797) q[0];
sx q[0];
rz(-0.63596361) q[0];
sx q[0];
rz(0.075809191) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71543773) q[2];
sx q[2];
rz(-1.1416661) q[2];
sx q[2];
rz(-1.5769584) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4574617) q[1];
sx q[1];
rz(-2.0254454) q[1];
sx q[1];
rz(-0.40747868) q[1];
rz(-2.5472574) q[3];
sx q[3];
rz(-0.96386516) q[3];
sx q[3];
rz(0.88483525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8422164) q[2];
sx q[2];
rz(-0.78263485) q[2];
sx q[2];
rz(-1.5807318) q[2];
rz(-0.98172274) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2340387) q[0];
sx q[0];
rz(-0.83775318) q[0];
sx q[0];
rz(0.87189829) q[0];
rz(2.7543228) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(-2.8505468) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29845333) q[0];
sx q[0];
rz(-1.342134) q[0];
sx q[0];
rz(-2.162621) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2497068) q[2];
sx q[2];
rz(-1.2920657) q[2];
sx q[2];
rz(-0.89821399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0709666) q[1];
sx q[1];
rz(-0.91849209) q[1];
sx q[1];
rz(2.2233306) q[1];
rz(-0.16104161) q[3];
sx q[3];
rz(-1.0003918) q[3];
sx q[3];
rz(2.963484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7724093) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(-0.54405653) q[2];
rz(-0.50928515) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(-0.86597401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13957025) q[0];
sx q[0];
rz(-2.1152571) q[0];
sx q[0];
rz(-0.64506662) q[0];
rz(-2.5158665) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(-2.5114139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7932927) q[0];
sx q[0];
rz(-2.5508946) q[0];
sx q[0];
rz(-2.9212055) q[0];
rz(-pi) q[1];
rz(2.3216256) q[2];
sx q[2];
rz(-1.8659288) q[2];
sx q[2];
rz(-0.68857869) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9471776) q[1];
sx q[1];
rz(-1.0395323) q[1];
sx q[1];
rz(2.8814425) q[1];
rz(-pi) q[2];
rz(1.7230677) q[3];
sx q[3];
rz(-1.0183183) q[3];
sx q[3];
rz(-3.1395562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3877635) q[2];
sx q[2];
rz(-1.8176273) q[2];
sx q[2];
rz(1.7774263) q[2];
rz(1.8917313) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(-0.92145872) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0548148) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(-2.4247647) q[0];
rz(-2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(0.81370083) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2700304) q[0];
sx q[0];
rz(-2.9189778) q[0];
sx q[0];
rz(-1.2827669) q[0];
rz(-0.72197638) q[2];
sx q[2];
rz(-1.1140031) q[2];
sx q[2];
rz(1.7320088) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8061714) q[1];
sx q[1];
rz(-0.85039447) q[1];
sx q[1];
rz(-2.6909188) q[1];
rz(0.86088647) q[3];
sx q[3];
rz(-2.6885899) q[3];
sx q[3];
rz(-2.3080254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2563236) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(1.1876594) q[2];
rz(-0.93196431) q[3];
sx q[3];
rz(-2.0075802) q[3];
sx q[3];
rz(2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1535783) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(-2.4555092) q[0];
rz(-1.6185121) q[1];
sx q[1];
rz(-0.97421092) q[1];
sx q[1];
rz(-2.4783321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2880741) q[0];
sx q[0];
rz(-0.92068866) q[0];
sx q[0];
rz(-3.1115565) q[0];
x q[1];
rz(-1.7887605) q[2];
sx q[2];
rz(-2.1969165) q[2];
sx q[2];
rz(1.8916212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38695947) q[1];
sx q[1];
rz(-0.52801758) q[1];
sx q[1];
rz(-1.1623043) q[1];
rz(-2.5797964) q[3];
sx q[3];
rz(-1.1601163) q[3];
sx q[3];
rz(-2.3464349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.474581) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(-1.0117426) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(-2.570178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052977234) q[0];
sx q[0];
rz(-2.4089854) q[0];
sx q[0];
rz(-0.10738871) q[0];
rz(-2.836851) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(-0.4577786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25154197) q[0];
sx q[0];
rz(-0.72351447) q[0];
sx q[0];
rz(-2.0680122) q[0];
rz(-pi) q[1];
rz(-2.5758178) q[2];
sx q[2];
rz(-2.3104295) q[2];
sx q[2];
rz(-2.4772252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1220408) q[1];
sx q[1];
rz(-1.8035839) q[1];
sx q[1];
rz(-1.2047393) q[1];
rz(-pi) q[2];
rz(-2.5896766) q[3];
sx q[3];
rz(-1.2593927) q[3];
sx q[3];
rz(-0.50406721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7765939) q[2];
sx q[2];
rz(-1.4564617) q[2];
sx q[2];
rz(-1.4314852) q[2];
rz(-1.4252023) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(-1.8553998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15437056) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(2.1881058) q[0];
rz(1.3257239) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(2.5568533) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9033602) q[0];
sx q[0];
rz(-2.0582709) q[0];
sx q[0];
rz(-1.6977298) q[0];
x q[1];
rz(-1.2051177) q[2];
sx q[2];
rz(-1.056991) q[2];
sx q[2];
rz(-2.6363381) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1251724) q[1];
sx q[1];
rz(-2.6965045) q[1];
sx q[1];
rz(-2.9801286) q[1];
rz(-pi) q[2];
rz(-2.9986266) q[3];
sx q[3];
rz(-1.8355832) q[3];
sx q[3];
rz(2.9395482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3328302) q[2];
sx q[2];
rz(-2.3213883) q[2];
sx q[2];
rz(2.8653223) q[2];
rz(-0.55109465) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(2.7579257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.0693414) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(1.1570702) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(1.5225333) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9690543) q[0];
sx q[0];
rz(-1.2612871) q[0];
sx q[0];
rz(-2.4597416) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7340639) q[2];
sx q[2];
rz(-0.73020836) q[2];
sx q[2];
rz(0.66116316) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5787705) q[1];
sx q[1];
rz(-2.4914503) q[1];
sx q[1];
rz(-0.27362089) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9736245) q[3];
sx q[3];
rz(-1.1790923) q[3];
sx q[3];
rz(0.21564461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15554252) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(2.299451) q[2];
rz(-1.5367674) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(-0.026528927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2355272) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(-0.72262598) q[1];
sx q[1];
rz(-2.2650748) q[1];
sx q[1];
rz(1.5320019) q[1];
rz(3.1149574) q[2];
sx q[2];
rz(-1.7799108) q[2];
sx q[2];
rz(-1.5540661) q[2];
rz(-1.6668672) q[3];
sx q[3];
rz(-0.69312743) q[3];
sx q[3];
rz(3.0467924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
