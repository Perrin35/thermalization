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
rz(-2.9057246) q[0];
sx q[0];
rz(-1.9960825) q[0];
sx q[0];
rz(3.0930162) q[0];
rz(1.5161169) q[1];
sx q[1];
rz(4.2344344) q[1];
sx q[1];
rz(8.7090127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5465821) q[0];
sx q[0];
rz(-1.5413741) q[0];
sx q[0];
rz(1.3999697) q[0];
rz(-2.2015436) q[2];
sx q[2];
rz(-0.34206451) q[2];
sx q[2];
rz(-2.2243481) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6977344) q[1];
sx q[1];
rz(-0.83699534) q[1];
sx q[1];
rz(-2.5910282) q[1];
x q[2];
rz(-1.1668901) q[3];
sx q[3];
rz(-2.0527186) q[3];
sx q[3];
rz(-1.9306077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6526661) q[2];
sx q[2];
rz(-0.55997866) q[2];
sx q[2];
rz(1.0555335) q[2];
rz(2.7355898) q[3];
sx q[3];
rz(-1.8279653) q[3];
sx q[3];
rz(0.82160151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15776289) q[0];
sx q[0];
rz(-2.1021748) q[0];
sx q[0];
rz(0.55617547) q[0];
rz(1.3293386) q[1];
sx q[1];
rz(-2.5228597) q[1];
sx q[1];
rz(-3.0583196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8502959) q[0];
sx q[0];
rz(-1.8938188) q[0];
sx q[0];
rz(-0.6685386) q[0];
rz(2.5640798) q[2];
sx q[2];
rz(-1.0339875) q[2];
sx q[2];
rz(0.84371882) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.00035380414) q[1];
sx q[1];
rz(-1.3508995) q[1];
sx q[1];
rz(-2.6882437) q[1];
rz(-pi) q[2];
rz(1.0038297) q[3];
sx q[3];
rz(-1.9975845) q[3];
sx q[3];
rz(0.7184283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9437774) q[2];
sx q[2];
rz(-1.752172) q[2];
sx q[2];
rz(2.6488292) q[2];
rz(2.1150151) q[3];
sx q[3];
rz(-0.30288282) q[3];
sx q[3];
rz(1.4125642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3781085) q[0];
sx q[0];
rz(-2.5788827) q[0];
sx q[0];
rz(-0.43877959) q[0];
rz(1.6890866) q[1];
sx q[1];
rz(-2.8197598) q[1];
sx q[1];
rz(2.7486393) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9742115) q[0];
sx q[0];
rz(-2.0618453) q[0];
sx q[0];
rz(-0.031662861) q[0];
rz(-pi) q[1];
rz(1.6268756) q[2];
sx q[2];
rz(-2.6125557) q[2];
sx q[2];
rz(-0.38944405) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.445791) q[1];
sx q[1];
rz(-1.6417827) q[1];
sx q[1];
rz(1.9893212) q[1];
x q[2];
rz(-0.10458979) q[3];
sx q[3];
rz(-0.65192332) q[3];
sx q[3];
rz(-0.014217941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68985933) q[2];
sx q[2];
rz(-2.2971575) q[2];
sx q[2];
rz(-0.27352697) q[2];
rz(0.61255974) q[3];
sx q[3];
rz(-1.7275683) q[3];
sx q[3];
rz(0.85649049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3278811) q[0];
sx q[0];
rz(-2.5484945) q[0];
sx q[0];
rz(-2.3408422) q[0];
rz(2.4294991) q[1];
sx q[1];
rz(-2.021603) q[1];
sx q[1];
rz(-0.48316479) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33813223) q[0];
sx q[0];
rz(-2.3134941) q[0];
sx q[0];
rz(-1.5716193) q[0];
x q[1];
rz(-1.0739378) q[2];
sx q[2];
rz(-1.0845079) q[2];
sx q[2];
rz(-0.72975791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.871328) q[1];
sx q[1];
rz(-0.92496269) q[1];
sx q[1];
rz(-1.6069018) q[1];
rz(0.68928257) q[3];
sx q[3];
rz(-1.3069671) q[3];
sx q[3];
rz(2.0949013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.157865) q[2];
sx q[2];
rz(-2.4309776) q[2];
sx q[2];
rz(1.857081) q[2];
rz(2.6537248) q[3];
sx q[3];
rz(-1.7043461) q[3];
sx q[3];
rz(-1.6871066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81216413) q[0];
sx q[0];
rz(-0.85947961) q[0];
sx q[0];
rz(-0.46464768) q[0];
rz(-2.1931785) q[1];
sx q[1];
rz(-0.83810884) q[1];
sx q[1];
rz(-1.4533739) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8226979) q[0];
sx q[0];
rz(-1.6056716) q[0];
sx q[0];
rz(0.025392763) q[0];
rz(-pi) q[1];
rz(-1.8273962) q[2];
sx q[2];
rz(-2.294696) q[2];
sx q[2];
rz(-0.43780299) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.075045295) q[1];
sx q[1];
rz(-2.5270224) q[1];
sx q[1];
rz(-0.11654186) q[1];
rz(-pi) q[2];
rz(-1.0539758) q[3];
sx q[3];
rz(-2.2272416) q[3];
sx q[3];
rz(2.2051728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.98443085) q[2];
sx q[2];
rz(-1.0474019) q[2];
sx q[2];
rz(-2.7471527) q[2];
rz(-1.1557584) q[3];
sx q[3];
rz(-2.211536) q[3];
sx q[3];
rz(1.3084779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30094576) q[0];
sx q[0];
rz(-2.7257958) q[0];
sx q[0];
rz(2.0796602) q[0];
rz(2.1901219) q[1];
sx q[1];
rz(-0.87363344) q[1];
sx q[1];
rz(-1.5207312) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011764048) q[0];
sx q[0];
rz(-1.6239663) q[0];
sx q[0];
rz(-1.2935678) q[0];
rz(-pi) q[1];
rz(2.4967147) q[2];
sx q[2];
rz(-2.0487924) q[2];
sx q[2];
rz(-2.1410112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6013612) q[1];
sx q[1];
rz(-1.4077606) q[1];
sx q[1];
rz(-0.77107471) q[1];
rz(-0.49455182) q[3];
sx q[3];
rz(-0.44122094) q[3];
sx q[3];
rz(-1.0426723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.025853) q[2];
sx q[2];
rz(-1.9151177) q[2];
sx q[2];
rz(-2.8958877) q[2];
rz(1.687441) q[3];
sx q[3];
rz(-1.629963) q[3];
sx q[3];
rz(-0.83034849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65569735) q[0];
sx q[0];
rz(-0.73807722) q[0];
sx q[0];
rz(-0.6947211) q[0];
rz(-0.10063902) q[1];
sx q[1];
rz(-2.1058319) q[1];
sx q[1];
rz(-0.86520854) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9964906) q[0];
sx q[0];
rz(-1.1528413) q[0];
sx q[0];
rz(-1.6161902) q[0];
rz(-2.2312282) q[2];
sx q[2];
rz(-1.5558793) q[2];
sx q[2];
rz(2.2098324) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3733417) q[1];
sx q[1];
rz(-0.96186559) q[1];
sx q[1];
rz(0.18280289) q[1];
rz(-2.4860704) q[3];
sx q[3];
rz(-0.78981863) q[3];
sx q[3];
rz(1.5721377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18409099) q[2];
sx q[2];
rz(-0.75347334) q[2];
sx q[2];
rz(-0.91020477) q[2];
rz(-1.6893859) q[3];
sx q[3];
rz(-1.9096749) q[3];
sx q[3];
rz(2.6773101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24377395) q[0];
sx q[0];
rz(-2.0938566) q[0];
sx q[0];
rz(2.5183103) q[0];
rz(-2.9858164) q[1];
sx q[1];
rz(-0.77729762) q[1];
sx q[1];
rz(2.8782841) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0000738) q[0];
sx q[0];
rz(-0.21528582) q[0];
sx q[0];
rz(1.6208642) q[0];
rz(-3.0650862) q[2];
sx q[2];
rz(-1.6311809) q[2];
sx q[2];
rz(-3.1186207) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38674445) q[1];
sx q[1];
rz(-1.6476908) q[1];
sx q[1];
rz(-2.5645178) q[1];
rz(-pi) q[2];
rz(-0.78531475) q[3];
sx q[3];
rz(-1.1668596) q[3];
sx q[3];
rz(0.31111003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.89173633) q[2];
sx q[2];
rz(-2.8920434) q[2];
sx q[2];
rz(-2.2535394) q[2];
rz(0.92787162) q[3];
sx q[3];
rz(-1.1279305) q[3];
sx q[3];
rz(2.1685062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4230098) q[0];
sx q[0];
rz(-2.6319478) q[0];
sx q[0];
rz(-0.50073671) q[0];
rz(2.0210733) q[1];
sx q[1];
rz(-2.3605533) q[1];
sx q[1];
rz(-0.80148554) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716176) q[0];
sx q[0];
rz(-1.5727676) q[0];
sx q[0];
rz(-3.0706186) q[0];
rz(-2.0960484) q[2];
sx q[2];
rz(-2.0727951) q[2];
sx q[2];
rz(1.0848622) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75569713) q[1];
sx q[1];
rz(-2.3508126) q[1];
sx q[1];
rz(-0.62259953) q[1];
rz(-pi) q[2];
rz(1.3128223) q[3];
sx q[3];
rz(-1.7593062) q[3];
sx q[3];
rz(-2.9303355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4215309) q[2];
sx q[2];
rz(-1.820182) q[2];
sx q[2];
rz(0.23019543) q[2];
rz(0.0097533334) q[3];
sx q[3];
rz(-1.516927) q[3];
sx q[3];
rz(-0.15596381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941876) q[0];
sx q[0];
rz(-1.7740086) q[0];
sx q[0];
rz(0.76960027) q[0];
rz(0.96254483) q[1];
sx q[1];
rz(-1.3178408) q[1];
sx q[1];
rz(0.98709551) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5819823) q[0];
sx q[0];
rz(-1.501069) q[0];
sx q[0];
rz(2.0980673) q[0];
x q[1];
rz(-0.15322282) q[2];
sx q[2];
rz(-1.8109057) q[2];
sx q[2];
rz(2.3599412) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0278575) q[1];
sx q[1];
rz(-0.052885508) q[1];
sx q[1];
rz(-0.89370339) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77642595) q[3];
sx q[3];
rz(-1.2033495) q[3];
sx q[3];
rz(2.020379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7139682) q[2];
sx q[2];
rz(-0.93728137) q[2];
sx q[2];
rz(3.0734708) q[2];
rz(2.6918329) q[3];
sx q[3];
rz(-0.41523784) q[3];
sx q[3];
rz(1.7033887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.149067) q[0];
sx q[0];
rz(-1.5328007) q[0];
sx q[0];
rz(-1.3971064) q[0];
rz(-0.29009157) q[1];
sx q[1];
rz(-2.7128704) q[1];
sx q[1];
rz(-2.6147978) q[1];
rz(-2.5505603) q[2];
sx q[2];
rz(-2.5606511) q[2];
sx q[2];
rz(-2.64369) q[2];
rz(2.7431106) q[3];
sx q[3];
rz(-1.553411) q[3];
sx q[3];
rz(-0.51391272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
