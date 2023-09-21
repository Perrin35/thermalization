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
rz(2.8770652) q[0];
sx q[0];
rz(9.8192083) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(-1.2000097) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7317176) q[0];
sx q[0];
rz(-0.048225064) q[0];
sx q[0];
rz(1.4531141) q[0];
rz(-0.33485246) q[2];
sx q[2];
rz(-1.9960253) q[2];
sx q[2];
rz(1.2644757) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66558054) q[1];
sx q[1];
rz(-1.952938) q[1];
sx q[1];
rz(1.0784472) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6148189) q[3];
sx q[3];
rz(-1.1929973) q[3];
sx q[3];
rz(1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(0.28960323) q[2];
rz(0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(3.0818821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(-0.34399024) q[0];
rz(-0.084331766) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.7864236) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25027572) q[0];
sx q[0];
rz(-0.96069562) q[0];
sx q[0];
rz(0.079205714) q[0];
x q[1];
rz(-1.6927035) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(0.09300692) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.018108798) q[1];
sx q[1];
rz(-1.6347486) q[1];
sx q[1];
rz(-1.121184) q[1];
x q[2];
rz(-0.077640688) q[3];
sx q[3];
rz(-1.3878126) q[3];
sx q[3];
rz(-1.7971355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8460059) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(-2.6039092) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(-1.0104377) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780592) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-1.0650939) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5168034) q[0];
sx q[0];
rz(-1.5644801) q[0];
sx q[0];
rz(1.9301231) q[0];
rz(-pi) q[1];
rz(-0.87419072) q[2];
sx q[2];
rz(-2.0421931) q[2];
sx q[2];
rz(-0.092560571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1467495) q[1];
sx q[1];
rz(-0.61770505) q[1];
sx q[1];
rz(-1.5084933) q[1];
rz(-pi) q[2];
rz(2.2933526) q[3];
sx q[3];
rz(-1.3812997) q[3];
sx q[3];
rz(3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.76434) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(-0.71737814) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(-0.29754105) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(2.8919019) q[0];
rz(1.0149792) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(0.011118523) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4269037) q[0];
sx q[0];
rz(-0.62781292) q[0];
sx q[0];
rz(2.3054302) q[0];
rz(-0.83799329) q[2];
sx q[2];
rz(-1.0655155) q[2];
sx q[2];
rz(0.2620286) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3938155) q[1];
sx q[1];
rz(-1.2443466) q[1];
sx q[1];
rz(1.0541037) q[1];
rz(-pi) q[2];
rz(1.3783781) q[3];
sx q[3];
rz(-1.8835861) q[3];
sx q[3];
rz(0.49945143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6461688) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(-2.8584976) q[2];
rz(-2.4781573) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0304612) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(2.8097613) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.3269075) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4715695) q[0];
sx q[0];
rz(-1.7243392) q[0];
sx q[0];
rz(2.3877386) q[0];
rz(-0.97557108) q[2];
sx q[2];
rz(-2.274548) q[2];
sx q[2];
rz(3.1095568) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8618968) q[1];
sx q[1];
rz(-1.8021291) q[1];
sx q[1];
rz(2.8054603) q[1];
rz(-2.1760686) q[3];
sx q[3];
rz(-2*pi/13) q[3];
sx q[3];
rz(2.0616639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.999324) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(1.8959321) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44928837) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(-2.4601049) q[0];
rz(-2.9340414) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(2.025827) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.241248) q[0];
sx q[0];
rz(-1.9755409) q[0];
sx q[0];
rz(-0.83165283) q[0];
rz(2.7115466) q[2];
sx q[2];
rz(-0.94656813) q[2];
sx q[2];
rz(0.66832322) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.15397729) q[1];
sx q[1];
rz(-2.1591641) q[1];
sx q[1];
rz(-1.5920559) q[1];
rz(-pi) q[2];
rz(-0.22478215) q[3];
sx q[3];
rz(-0.68192476) q[3];
sx q[3];
rz(3.0608321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21268022) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(-2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0091244) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(-2.5792714) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0494941) q[0];
sx q[0];
rz(-1.5856165) q[0];
sx q[0];
rz(-2.2390319) q[0];
rz(-pi) q[1];
rz(-1.6717031) q[2];
sx q[2];
rz(-1.9011874) q[2];
sx q[2];
rz(0.23320564) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.454969) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(-1.7186233) q[1];
rz(-pi) q[2];
rz(2.167422) q[3];
sx q[3];
rz(-1.6430292) q[3];
sx q[3];
rz(0.40856397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0397296) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(-2.8015461) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44889221) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(2.7451519) q[0];
rz(3.0026644) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(1.5213535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0884468) q[0];
sx q[0];
rz(-1.3741115) q[0];
sx q[0];
rz(1.501207) q[0];
rz(-0.40201681) q[2];
sx q[2];
rz(-2.0094299) q[2];
sx q[2];
rz(-1.6059665) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.24890451) q[1];
sx q[1];
rz(-2.3612594) q[1];
sx q[1];
rz(-0.08463879) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11407125) q[3];
sx q[3];
rz(-0.89469203) q[3];
sx q[3];
rz(1.8166208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5788995) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(-2.8472624) q[2];
rz(-1.1307905) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(0.35933581) q[0];
rz(2.1954779) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(0.27063453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7662738) q[0];
sx q[0];
rz(-0.047485504) q[0];
sx q[0];
rz(-2.4256698) q[0];
x q[1];
rz(-2.3324899) q[2];
sx q[2];
rz(-1.9168233) q[2];
sx q[2];
rz(1.5835294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0535069) q[1];
sx q[1];
rz(-0.063789531) q[1];
sx q[1];
rz(2.3699058) q[1];
rz(-pi) q[2];
rz(2.4828033) q[3];
sx q[3];
rz(-0.82595347) q[3];
sx q[3];
rz(-2.7487019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(-0.81196249) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6311326) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(2.4023138) q[0];
rz(-2.9108858) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(0.49490067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5737168) q[0];
sx q[0];
rz(-2.0567237) q[0];
sx q[0];
rz(3.0987415) q[0];
rz(-pi) q[1];
rz(1.5636744) q[2];
sx q[2];
rz(-1.2429503) q[2];
sx q[2];
rz(-0.42051007) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6634076) q[1];
sx q[1];
rz(-0.93485281) q[1];
sx q[1];
rz(0.138476) q[1];
rz(-0.095330843) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(-1.8180083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13359244) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(2.949529) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(2.1081934) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-2.1140425) q[2];
sx q[2];
rz(-1.3941358) q[2];
sx q[2];
rz(0.91640581) q[2];
rz(2.8752747) q[3];
sx q[3];
rz(-1.3238293) q[3];
sx q[3];
rz(-2.7662591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];