OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12299744) q[0];
sx q[0];
rz(-1.7641492) q[0];
sx q[0];
rz(-2.7121845) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(7.9908854) q[1];
sx q[1];
rz(7.8828852) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2875007) q[0];
sx q[0];
rz(-1.8563636) q[0];
sx q[0];
rz(2.9182316) q[0];
rz(-pi) q[1];
rz(2.9978031) q[2];
sx q[2];
rz(-1.6056639) q[2];
sx q[2];
rz(-0.26511017) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0860025) q[1];
sx q[1];
rz(-1.1430972) q[1];
sx q[1];
rz(-1.5375288) q[1];
rz(0.67486169) q[3];
sx q[3];
rz(-0.45025846) q[3];
sx q[3];
rz(-0.78801149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8371007) q[2];
sx q[2];
rz(-1.203275) q[2];
sx q[2];
rz(-0.692918) q[2];
rz(2.9414226) q[3];
sx q[3];
rz(-2.9514511) q[3];
sx q[3];
rz(-1.1339124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3226586) q[0];
sx q[0];
rz(-2.942473) q[0];
sx q[0];
rz(-1.06485) q[0];
rz(2.5191567) q[1];
sx q[1];
rz(-1.348) q[1];
sx q[1];
rz(-2.650824) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1932632) q[0];
sx q[0];
rz(-1.5800467) q[0];
sx q[0];
rz(-0.023764334) q[0];
rz(2.6726637) q[2];
sx q[2];
rz(-2.6750419) q[2];
sx q[2];
rz(0.55794898) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8339757) q[1];
sx q[1];
rz(-2.3217891) q[1];
sx q[1];
rz(-1.7012672) q[1];
x q[2];
rz(-1.2964046) q[3];
sx q[3];
rz(-0.17582045) q[3];
sx q[3];
rz(2.4309513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6530767) q[2];
sx q[2];
rz(-1.4986897) q[2];
sx q[2];
rz(-0.24031362) q[2];
rz(2.6416685) q[3];
sx q[3];
rz(-0.58568716) q[3];
sx q[3];
rz(1.8723996) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2568473) q[0];
sx q[0];
rz(-2.186543) q[0];
sx q[0];
rz(-2.1599059) q[0];
rz(0.49199545) q[1];
sx q[1];
rz(-2.056608) q[1];
sx q[1];
rz(-1.5031776) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7273151) q[0];
sx q[0];
rz(-0.87782598) q[0];
sx q[0];
rz(0.63590617) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0122224) q[2];
sx q[2];
rz(-0.40310848) q[2];
sx q[2];
rz(2.108824) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99276421) q[1];
sx q[1];
rz(-1.9566613) q[1];
sx q[1];
rz(0.35986114) q[1];
rz(-pi) q[2];
rz(-3.1206661) q[3];
sx q[3];
rz(-1.4731579) q[3];
sx q[3];
rz(0.00033631246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.65695277) q[2];
sx q[2];
rz(-2.6159365) q[2];
sx q[2];
rz(0.12118113) q[2];
rz(2.1335404) q[3];
sx q[3];
rz(-2.0262148) q[3];
sx q[3];
rz(2.0602267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0693531) q[0];
sx q[0];
rz(-0.00088748137) q[0];
sx q[0];
rz(0.51373154) q[0];
rz(2.1836102) q[1];
sx q[1];
rz(-1.1424516) q[1];
sx q[1];
rz(-1.4428008) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4289249) q[0];
sx q[0];
rz(-0.44433549) q[0];
sx q[0];
rz(2.7508468) q[0];
x q[1];
rz(-0.36607217) q[2];
sx q[2];
rz(-2.5806081) q[2];
sx q[2];
rz(1.6268886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0569257) q[1];
sx q[1];
rz(-0.62090331) q[1];
sx q[1];
rz(-1.4292745) q[1];
x q[2];
rz(0.12421457) q[3];
sx q[3];
rz(-0.66541568) q[3];
sx q[3];
rz(-1.6197325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5865667) q[2];
sx q[2];
rz(-0.96379605) q[2];
sx q[2];
rz(-1.2653992) q[2];
rz(1.966656) q[3];
sx q[3];
rz(-0.73052162) q[3];
sx q[3];
rz(-1.1635715) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88605276) q[0];
sx q[0];
rz(-0.20296725) q[0];
sx q[0];
rz(-1.8170005) q[0];
rz(-0.96877226) q[1];
sx q[1];
rz(-1.6561534) q[1];
sx q[1];
rz(-1.0383777) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7234112) q[0];
sx q[0];
rz(-1.9514284) q[0];
sx q[0];
rz(-1.7443954) q[0];
x q[1];
rz(-2.3275787) q[2];
sx q[2];
rz(-0.90384877) q[2];
sx q[2];
rz(2.2499354) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9818791) q[1];
sx q[1];
rz(-0.37283373) q[1];
sx q[1];
rz(-0.11319686) q[1];
rz(-1.0594924) q[3];
sx q[3];
rz(-0.97701525) q[3];
sx q[3];
rz(2.4932414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6394627) q[2];
sx q[2];
rz(-2.837193) q[2];
sx q[2];
rz(0.89140618) q[2];
rz(-1.8799479) q[3];
sx q[3];
rz(-1.5104048) q[3];
sx q[3];
rz(2.4752899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0618184) q[0];
sx q[0];
rz(-0.49538716) q[0];
sx q[0];
rz(-2.1451982) q[0];
rz(0.90006104) q[1];
sx q[1];
rz(-0.78949094) q[1];
sx q[1];
rz(0.17734227) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9437761) q[0];
sx q[0];
rz(-1.7223486) q[0];
sx q[0];
rz(0.00023784296) q[0];
x q[1];
rz(3.1204865) q[2];
sx q[2];
rz(-0.93230196) q[2];
sx q[2];
rz(0.47551949) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2832433) q[1];
sx q[1];
rz(-0.56811404) q[1];
sx q[1];
rz(-1.0198221) q[1];
rz(-pi) q[2];
rz(0.097058724) q[3];
sx q[3];
rz(-2.0776111) q[3];
sx q[3];
rz(2.7926796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8511054) q[2];
sx q[2];
rz(-1.1815716) q[2];
sx q[2];
rz(-0.30430749) q[2];
rz(-2.815222) q[3];
sx q[3];
rz(-2.5984952) q[3];
sx q[3];
rz(2.2296026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0684763) q[0];
sx q[0];
rz(-2.731972) q[0];
sx q[0];
rz(-0.9325183) q[0];
rz(-0.069123507) q[1];
sx q[1];
rz(-0.2176452) q[1];
sx q[1];
rz(-1.0401915) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0737586) q[0];
sx q[0];
rz(-2.375575) q[0];
sx q[0];
rz(-3.104408) q[0];
rz(-pi) q[1];
rz(2.1625453) q[2];
sx q[2];
rz(-0.30610105) q[2];
sx q[2];
rz(-3.0476168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5835147) q[1];
sx q[1];
rz(-2.3574074) q[1];
sx q[1];
rz(2.9072443) q[1];
rz(-pi) q[2];
rz(0.69871083) q[3];
sx q[3];
rz(-2.0255528) q[3];
sx q[3];
rz(-0.48331279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13504623) q[2];
sx q[2];
rz(-0.73837215) q[2];
sx q[2];
rz(-2.4086003) q[2];
rz(2.0586355) q[3];
sx q[3];
rz(-1.0561918) q[3];
sx q[3];
rz(1.0068896) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55423823) q[0];
sx q[0];
rz(-3.0576958) q[0];
sx q[0];
rz(-2.6485637) q[0];
rz(-0.21656491) q[1];
sx q[1];
rz(-2.000587) q[1];
sx q[1];
rz(-1.0239748) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39505491) q[0];
sx q[0];
rz(-0.97486541) q[0];
sx q[0];
rz(1.5851969) q[0];
x q[1];
rz(-2.8220213) q[2];
sx q[2];
rz(-0.34863483) q[2];
sx q[2];
rz(2.2120668) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9453911) q[1];
sx q[1];
rz(-1.7111254) q[1];
sx q[1];
rz(2.7811471) q[1];
rz(-pi) q[2];
rz(2.2129503) q[3];
sx q[3];
rz(-1.3285884) q[3];
sx q[3];
rz(0.95706576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0852802) q[2];
sx q[2];
rz(-1.5093426) q[2];
sx q[2];
rz(0.67576605) q[2];
rz(-2.7285649) q[3];
sx q[3];
rz(-1.690381) q[3];
sx q[3];
rz(-0.54615027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6290879) q[0];
sx q[0];
rz(-0.65852037) q[0];
sx q[0];
rz(3.107048) q[0];
rz(-0.054556219) q[1];
sx q[1];
rz(-1.9787534) q[1];
sx q[1];
rz(-2.3202855) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4865882) q[0];
sx q[0];
rz(-1.6704217) q[0];
sx q[0];
rz(2.5667046) q[0];
rz(1.5615145) q[2];
sx q[2];
rz(-1.4820274) q[2];
sx q[2];
rz(0.87997681) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2149722) q[1];
sx q[1];
rz(-1.4550303) q[1];
sx q[1];
rz(2.9909219) q[1];
rz(-pi) q[2];
rz(0.48252941) q[3];
sx q[3];
rz(-0.97597968) q[3];
sx q[3];
rz(-1.9687259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.36046946) q[2];
sx q[2];
rz(-2.1253773) q[2];
sx q[2];
rz(0.13295573) q[2];
rz(-0.41108701) q[3];
sx q[3];
rz(-1.6229595) q[3];
sx q[3];
rz(-1.0857922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046831176) q[0];
sx q[0];
rz(-0.60698858) q[0];
sx q[0];
rz(1.2588311) q[0];
rz(1.5065441) q[1];
sx q[1];
rz(-2.4848487) q[1];
sx q[1];
rz(-2.3283995) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5482322) q[0];
sx q[0];
rz(-1.5357067) q[0];
sx q[0];
rz(0.078087383) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7178463) q[2];
sx q[2];
rz(-2.4644901) q[2];
sx q[2];
rz(0.33237095) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5845994) q[1];
sx q[1];
rz(-2.0210588) q[1];
sx q[1];
rz(-2.492639) q[1];
rz(1.7045295) q[3];
sx q[3];
rz(-1.1220891) q[3];
sx q[3];
rz(2.1577842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64858156) q[2];
sx q[2];
rz(-1.1375256) q[2];
sx q[2];
rz(1.0515593) q[2];
rz(2.1227396) q[3];
sx q[3];
rz(-2.2994883) q[3];
sx q[3];
rz(-2.4837608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9857585) q[0];
sx q[0];
rz(-0.65407615) q[0];
sx q[0];
rz(0.88055897) q[0];
rz(-1.3195994) q[1];
sx q[1];
rz(-1.9965912) q[1];
sx q[1];
rz(0.051864787) q[1];
rz(-0.32991275) q[2];
sx q[2];
rz(-2.920426) q[2];
sx q[2];
rz(-2.5830808) q[2];
rz(1.6307835) q[3];
sx q[3];
rz(-1.3740448) q[3];
sx q[3];
rz(-3.1136953) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
