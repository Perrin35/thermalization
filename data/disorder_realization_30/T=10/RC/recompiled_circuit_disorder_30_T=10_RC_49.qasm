OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(-0.2645275) q[0];
sx q[0];
rz(-0.39443031) q[0];
rz(-3.1354304) q[1];
sx q[1];
rz(-2.8013464) q[1];
sx q[1];
rz(-1.9415829) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8631247) q[0];
sx q[0];
rz(-1.5651363) q[0];
sx q[0];
rz(-1.6186884) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94304396) q[2];
sx q[2];
rz(-2.606751) q[2];
sx q[2];
rz(-0.56378555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29772273) q[1];
sx q[1];
rz(-2.5281318) q[1];
sx q[1];
rz(-0.86617275) q[1];
rz(-2.7634611) q[3];
sx q[3];
rz(-1.5298801) q[3];
sx q[3];
rz(3.0331628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(0.28960323) q[2];
rz(-0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(0.084331766) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.3551691) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2750759) q[0];
sx q[0];
rz(-1.6356902) q[0];
sx q[0];
rz(2.1823723) q[0];
x q[1];
rz(-0.64135735) q[2];
sx q[2];
rz(-0.20198447) q[2];
sx q[2];
rz(-0.55822492) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4210216) q[1];
sx q[1];
rz(-2.6877626) q[1];
sx q[1];
rz(1.4245016) q[1];
rz(-pi) q[2];
rz(0.077640688) q[3];
sx q[3];
rz(-1.3878126) q[3];
sx q[3];
rz(-1.3444572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8460059) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(-0.53768349) q[2];
rz(0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4780592) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-0.74869853) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-2.0764988) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0707866) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(1.5887567) q[0];
rz(-pi) q[1];
rz(0.89945729) q[2];
sx q[2];
rz(-0.81842917) q[2];
sx q[2];
rz(-2.1607272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6667337) q[1];
sx q[1];
rz(-1.6068646) q[1];
sx q[1];
rz(2.1875847) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8912656) q[3];
sx q[3];
rz(-0.86391376) q[3];
sx q[3];
rz(1.7771429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.76434) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(2.4242145) q[2];
rz(-2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(-0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82729572) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(2.8919019) q[0];
rz(-2.1266134) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(-0.011118523) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5549094) q[0];
sx q[0];
rz(-1.1197829) q[0];
sx q[0];
rz(-2.6888072) q[0];
x q[1];
rz(-0.87984933) q[2];
sx q[2];
rz(-0.86266154) q[2];
sx q[2];
rz(-0.81530064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3938155) q[1];
sx q[1];
rz(-1.897246) q[1];
sx q[1];
rz(-1.0541037) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3783781) q[3];
sx q[3];
rz(-1.2580066) q[3];
sx q[3];
rz(0.49945143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.49542385) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(2.8584976) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0304612) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.3269075) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0622612) q[0];
sx q[0];
rz(-0.766303) q[0];
sx q[0];
rz(2.9192231) q[0];
rz(-2.3438498) q[2];
sx q[2];
rz(-1.1290871) q[2];
sx q[2];
rz(-1.1898578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8718308) q[1];
sx q[1];
rz(-2.7360536) q[1];
sx q[1];
rz(2.521442) q[1];
rz(0.29019659) q[3];
sx q[3];
rz(-1.962933) q[3];
sx q[3];
rz(0.41662595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(-1.8959321) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44928837) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(0.6814878) q[0];
rz(2.9340414) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(1.1157657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078243144) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(1.0043762) q[0];
rz(-2.7115466) q[2];
sx q[2];
rz(-2.1950245) q[2];
sx q[2];
rz(-2.4732694) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11568497) q[1];
sx q[1];
rz(-2.5528862) q[1];
sx q[1];
rz(3.1097417) q[1];
rz(0.22478215) q[3];
sx q[3];
rz(-0.68192476) q[3];
sx q[3];
rz(-3.0608321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-2.7872655) q[2];
rz(-2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0091244) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(-1.1122423) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(2.5792714) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6816662) q[0];
sx q[0];
rz(-2.473218) q[0];
sx q[0];
rz(1.5947123) q[0];
rz(-pi) q[1];
rz(-1.4698896) q[2];
sx q[2];
rz(-1.2404053) q[2];
sx q[2];
rz(-2.908387) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.454969) q[1];
sx q[1];
rz(-1.7640055) q[1];
sx q[1];
rz(1.4229694) q[1];
rz(1.698877) q[3];
sx q[3];
rz(-0.60045419) q[3];
sx q[3];
rz(-1.8734991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(-2.8015461) q[2];
rz(-0.21751054) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927004) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-2.7451519) q[0];
rz(-3.0026644) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(1.6202392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5312708) q[0];
sx q[0];
rz(-1.5025508) q[0];
sx q[0];
rz(-2.9444429) q[0];
rz(2.0422158) q[2];
sx q[2];
rz(-1.9328914) q[2];
sx q[2];
rz(-0.21381703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13008598) q[1];
sx q[1];
rz(-0.79400051) q[1];
sx q[1];
rz(-1.4873051) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2500854) q[3];
sx q[3];
rz(-1.4818947) q[3];
sx q[3];
rz(2.8241983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56269318) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(-0.29433027) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(0.94611478) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(2.8709581) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3753189) q[0];
sx q[0];
rz(-0.047485504) q[0];
sx q[0];
rz(-2.4256698) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6793849) q[2];
sx q[2];
rz(-2.2773829) q[2];
sx q[2];
rz(2.8414937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2808025) q[1];
sx q[1];
rz(-1.5250912) q[1];
sx q[1];
rz(-1.6153107) q[1];
rz(2.4326914) q[3];
sx q[3];
rz(-2.0376251) q[3];
sx q[3];
rz(-2.4469568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82751194) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(2.6861526) q[2];
rz(0.81196249) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(-0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51046002) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(0.23070681) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(-2.646692) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56787581) q[0];
sx q[0];
rz(-1.084869) q[0];
sx q[0];
rz(0.042851187) q[0];
rz(-1.5779183) q[2];
sx q[2];
rz(-1.8986423) q[2];
sx q[2];
rz(0.42051007) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1315688) q[1];
sx q[1];
rz(-1.459517) q[1];
sx q[1];
rz(-0.93025031) q[1];
rz(-pi) q[2];
rz(-3.0462618) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(1.8180083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13359244) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(-2.949529) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(-1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1223758) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(2.3090251) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(2.9359948) q[2];
sx q[2];
rz(-2.1046706) q[2];
sx q[2];
rz(-0.54866366) q[2];
rz(-1.3151863) q[3];
sx q[3];
rz(-1.3127463) q[3];
sx q[3];
rz(-1.128872) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
