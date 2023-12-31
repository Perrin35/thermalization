OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7944613) q[0];
sx q[0];
rz(-2.1262655) q[0];
sx q[0];
rz(-0.46749687) q[0];
rz(2.6214018) q[1];
sx q[1];
rz(-1.7953035) q[1];
sx q[1];
rz(-2.4612114) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2517393) q[0];
sx q[0];
rz(-1.7008874) q[0];
sx q[0];
rz(-0.84530172) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19619588) q[2];
sx q[2];
rz(-1.2805403) q[2];
sx q[2];
rz(2.6930075) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3375652) q[1];
sx q[1];
rz(-2.0731888) q[1];
sx q[1];
rz(0.2663836) q[1];
rz(-pi) q[2];
rz(-0.33403553) q[3];
sx q[3];
rz(-1.738027) q[3];
sx q[3];
rz(-0.86080307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9900069) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(0.72922373) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-1.0013642) q[0];
rz(-2.9691866) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(-0.52406812) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66602089) q[0];
sx q[0];
rz(-2.4040939) q[0];
sx q[0];
rz(-1.1230418) q[0];
rz(2.7153154) q[2];
sx q[2];
rz(-0.91765109) q[2];
sx q[2];
rz(3.1003568) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14428917) q[1];
sx q[1];
rz(-1.3376437) q[1];
sx q[1];
rz(-0.38645978) q[1];
rz(-pi) q[2];
rz(2.6045813) q[3];
sx q[3];
rz(-1.4376663) q[3];
sx q[3];
rz(0.46768026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9841763) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(1.3400419) q[2];
rz(0.79483461) q[3];
sx q[3];
rz(-1.1398311) q[3];
sx q[3];
rz(-0.036858233) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3443417) q[0];
sx q[0];
rz(-0.69673711) q[0];
sx q[0];
rz(2.5168193) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(-2.952081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0441168) q[0];
sx q[0];
rz(-2.2783845) q[0];
sx q[0];
rz(1.2114899) q[0];
rz(-0.37522845) q[2];
sx q[2];
rz(-1.4725176) q[2];
sx q[2];
rz(2.3972942) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5901958) q[1];
sx q[1];
rz(-1.3314178) q[1];
sx q[1];
rz(-1.5134642) q[1];
rz(-pi) q[2];
rz(-1.8241006) q[3];
sx q[3];
rz(-2.5251303) q[3];
sx q[3];
rz(-2.9951819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1069964) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(-1.754388) q[2];
rz(-0.43131367) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(1.0881902) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65961924) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(-1.5455998) q[0];
rz(2.003147) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(1.0901573) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50463146) q[0];
sx q[0];
rz(-1.22389) q[0];
sx q[0];
rz(-2.7142801) q[0];
rz(-pi) q[1];
rz(2.3218669) q[2];
sx q[2];
rz(-1.956454) q[2];
sx q[2];
rz(-2.0476598) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5174597) q[1];
sx q[1];
rz(-2.1412666) q[1];
sx q[1];
rz(-1.073451) q[1];
rz(-pi) q[2];
x q[2];
rz(1.777321) q[3];
sx q[3];
rz(-1.9897451) q[3];
sx q[3];
rz(-1.2057613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1426992) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(0.85582716) q[2];
rz(1.9479729) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(-2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(-2.9918616) q[0];
rz(-2.1504452) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.9715462) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4677508) q[0];
sx q[0];
rz(-1.8918599) q[0];
sx q[0];
rz(2.2348316) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1779551) q[2];
sx q[2];
rz(-2.422214) q[2];
sx q[2];
rz(2.9658085) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7331446) q[1];
sx q[1];
rz(-2.0093579) q[1];
sx q[1];
rz(2.0088197) q[1];
rz(-pi) q[2];
rz(0.31801362) q[3];
sx q[3];
rz(-0.56088305) q[3];
sx q[3];
rz(-0.47919264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(-0.13723792) q[2];
rz(1.8042701) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(-1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.18734922) q[0];
sx q[0];
rz(-1.3784778) q[0];
sx q[0];
rz(-1.7911918) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(2.4687016) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85149375) q[0];
sx q[0];
rz(-2.5585472) q[0];
sx q[0];
rz(-1.3821938) q[0];
rz(-pi) q[1];
x q[1];
rz(0.024725155) q[2];
sx q[2];
rz(-2.284986) q[2];
sx q[2];
rz(-1.6753472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2313378) q[1];
sx q[1];
rz(-2.5006223) q[1];
sx q[1];
rz(-0.35229589) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9665885) q[3];
sx q[3];
rz(-0.41737469) q[3];
sx q[3];
rz(-2.008703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3488397) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(1.7815636) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(0.24766651) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9687987) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(-2.0794179) q[0];
rz(-2.0299714) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(1.7395082) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8334478) q[0];
sx q[0];
rz(-1.2045367) q[0];
sx q[0];
rz(1.863198) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3767396) q[2];
sx q[2];
rz(-1.3991038) q[2];
sx q[2];
rz(2.9219251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78571037) q[1];
sx q[1];
rz(-1.7240228) q[1];
sx q[1];
rz(-0.36550891) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13774638) q[3];
sx q[3];
rz(-2.1546954) q[3];
sx q[3];
rz(-3.0018842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7193675) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(2.4613703) q[2];
rz(-2.7133572) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(1.7817106) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5087886) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(-1.5493786) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.1523749) q[1];
sx q[1];
rz(-2.6002398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4562948) q[0];
sx q[0];
rz(-1.0066427) q[0];
sx q[0];
rz(-1.5575404) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5112553) q[2];
sx q[2];
rz(-2.6001843) q[2];
sx q[2];
rz(2.2556717) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.036877) q[1];
sx q[1];
rz(-0.86537213) q[1];
sx q[1];
rz(-2.5887262) q[1];
rz(-0.3803216) q[3];
sx q[3];
rz(-1.5725279) q[3];
sx q[3];
rz(2.1960432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0970739) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-2.3642335) q[2];
rz(-0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(-1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38354307) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(1.0154356) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(1.6962956) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24628595) q[0];
sx q[0];
rz(-0.9406957) q[0];
sx q[0];
rz(2.2340328) q[0];
rz(2.5601013) q[2];
sx q[2];
rz(-2.1610689) q[2];
sx q[2];
rz(1.3587388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.054113764) q[1];
sx q[1];
rz(-2.2009654) q[1];
sx q[1];
rz(1.5900882) q[1];
x q[2];
rz(2.2441838) q[3];
sx q[3];
rz(-2.071278) q[3];
sx q[3];
rz(0.15411479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4426667) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(-0.60738579) q[2];
rz(-1.4390885) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8059175) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(0.6138531) q[0];
rz(-1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(-2.7493431) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046767226) q[0];
sx q[0];
rz(-1.9691159) q[0];
sx q[0];
rz(2.4368068) q[0];
rz(-pi) q[1];
rz(0.62551542) q[2];
sx q[2];
rz(-1.0872456) q[2];
sx q[2];
rz(-0.12550437) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1019198) q[1];
sx q[1];
rz(-2.3452248) q[1];
sx q[1];
rz(2.2069195) q[1];
rz(0.20499968) q[3];
sx q[3];
rz(-1.5424171) q[3];
sx q[3];
rz(1.320822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0739416) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(0.60047853) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4338715) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(-2.1451163) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(1.2148576) q[2];
sx q[2];
rz(-0.83649737) q[2];
sx q[2];
rz(2.5964824) q[2];
rz(0.89041238) q[3];
sx q[3];
rz(-1.1287516) q[3];
sx q[3];
rz(0.26656084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
