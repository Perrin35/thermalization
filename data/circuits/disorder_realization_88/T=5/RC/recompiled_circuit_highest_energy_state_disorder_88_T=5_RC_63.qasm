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
rz(0.54230827) q[0];
sx q[0];
rz(-0.13442726) q[0];
sx q[0];
rz(2.0943213) q[0];
rz(-0.32416999) q[1];
sx q[1];
rz(-0.14961641) q[1];
sx q[1];
rz(-2.1967998) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1143998) q[0];
sx q[0];
rz(-1.8437705) q[0];
sx q[0];
rz(1.9045715) q[0];
rz(-pi) q[1];
rz(0.23350291) q[2];
sx q[2];
rz(-2.2036607) q[2];
sx q[2];
rz(-1.7720745) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85192363) q[1];
sx q[1];
rz(-1.1186558) q[1];
sx q[1];
rz(-0.55638066) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0202623) q[3];
sx q[3];
rz(-1.0912195) q[3];
sx q[3];
rz(-0.12313719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2904539) q[2];
sx q[2];
rz(-1.4522499) q[2];
sx q[2];
rz(-1.5770844) q[2];
rz(-0.56646281) q[3];
sx q[3];
rz(-2.5201859) q[3];
sx q[3];
rz(1.2982781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8241149) q[0];
sx q[0];
rz(-0.7190187) q[0];
sx q[0];
rz(-0.7199921) q[0];
rz(-0.27436817) q[1];
sx q[1];
rz(-0.73148483) q[1];
sx q[1];
rz(-2.5879587) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8047129) q[0];
sx q[0];
rz(-1.9229445) q[0];
sx q[0];
rz(3.0180305) q[0];
rz(0.64525147) q[2];
sx q[2];
rz(-0.49075365) q[2];
sx q[2];
rz(0.09749271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.95018847) q[1];
sx q[1];
rz(-0.85655115) q[1];
sx q[1];
rz(2.6493303) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6194044) q[3];
sx q[3];
rz(-2.0117674) q[3];
sx q[3];
rz(-2.9984634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1729892) q[2];
sx q[2];
rz(-1.9393238) q[2];
sx q[2];
rz(-2.7776862) q[2];
rz(-3.0299752) q[3];
sx q[3];
rz(-0.93897396) q[3];
sx q[3];
rz(-2.3811471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.97742057) q[0];
sx q[0];
rz(-0.82472473) q[0];
sx q[0];
rz(0.0053996276) q[0];
rz(-0.81575704) q[1];
sx q[1];
rz(-2.0033671) q[1];
sx q[1];
rz(-2.7556748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9078131) q[0];
sx q[0];
rz(-0.44615567) q[0];
sx q[0];
rz(-2.4491007) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0089985) q[2];
sx q[2];
rz(-1.8799861) q[2];
sx q[2];
rz(-0.55436347) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.94066835) q[1];
sx q[1];
rz(-1.4438119) q[1];
sx q[1];
rz(-1.1881571) q[1];
rz(-pi) q[2];
rz(3.0073062) q[3];
sx q[3];
rz(-1.1597753) q[3];
sx q[3];
rz(-2.3255796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5526814) q[2];
sx q[2];
rz(-2.3083355) q[2];
sx q[2];
rz(-1.7595278) q[2];
rz(-0.16119334) q[3];
sx q[3];
rz(-1.7101945) q[3];
sx q[3];
rz(-2.5126357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2564119) q[0];
sx q[0];
rz(-1.2762524) q[0];
sx q[0];
rz(-1.3389583) q[0];
rz(-1.3171875) q[1];
sx q[1];
rz(-2.1664797) q[1];
sx q[1];
rz(2.8417974) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2071964) q[0];
sx q[0];
rz(-2.3863433) q[0];
sx q[0];
rz(2.62225) q[0];
rz(2.2251769) q[2];
sx q[2];
rz(-0.70087003) q[2];
sx q[2];
rz(0.13417164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2721709) q[1];
sx q[1];
rz(-1.7079087) q[1];
sx q[1];
rz(-1.2161944) q[1];
x q[2];
rz(1.3922774) q[3];
sx q[3];
rz(-2.4981294) q[3];
sx q[3];
rz(-1.3993128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.51912159) q[2];
sx q[2];
rz(-0.93834472) q[2];
sx q[2];
rz(-2.8224714) q[2];
rz(3.0506813) q[3];
sx q[3];
rz(-0.14486434) q[3];
sx q[3];
rz(1.3566141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262481) q[0];
sx q[0];
rz(-2.321796) q[0];
sx q[0];
rz(3.1392642) q[0];
rz(1.5709411) q[1];
sx q[1];
rz(-0.80051533) q[1];
sx q[1];
rz(-1.194582) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2187506) q[0];
sx q[0];
rz(-2.9463769) q[0];
sx q[0];
rz(2.4890635) q[0];
rz(2.2885913) q[2];
sx q[2];
rz(-1.5403131) q[2];
sx q[2];
rz(-0.76678813) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4264226) q[1];
sx q[1];
rz(-1.914901) q[1];
sx q[1];
rz(-0.10870966) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2890599) q[3];
sx q[3];
rz(-2.1093371) q[3];
sx q[3];
rz(2.0592214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2583367) q[2];
sx q[2];
rz(-0.77815762) q[2];
sx q[2];
rz(3.0641595) q[2];
rz(2.1954913) q[3];
sx q[3];
rz(-2.6587722) q[3];
sx q[3];
rz(-0.4536804) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2133863) q[0];
sx q[0];
rz(-3.0551857) q[0];
sx q[0];
rz(1.7267831) q[0];
rz(0.66697031) q[1];
sx q[1];
rz(-1.7844776) q[1];
sx q[1];
rz(0.40474969) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0694097) q[0];
sx q[0];
rz(-1.2889922) q[0];
sx q[0];
rz(-0.3148766) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91333977) q[2];
sx q[2];
rz(-1.7796675) q[2];
sx q[2];
rz(0.058519017) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.09913832) q[1];
sx q[1];
rz(-2.6746866) q[1];
sx q[1];
rz(1.7068638) q[1];
rz(-2.1920106) q[3];
sx q[3];
rz(-0.57692617) q[3];
sx q[3];
rz(-2.2872137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3177967) q[2];
sx q[2];
rz(-2.3522289) q[2];
sx q[2];
rz(0.44343597) q[2];
rz(0.41380841) q[3];
sx q[3];
rz(-2.8518854) q[3];
sx q[3];
rz(-2.2558291) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7498748) q[0];
sx q[0];
rz(-1.3626008) q[0];
sx q[0];
rz(2.474127) q[0];
rz(-0.04774566) q[1];
sx q[1];
rz(-2.0011438) q[1];
sx q[1];
rz(2.023229) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76343173) q[0];
sx q[0];
rz(-1.1677647) q[0];
sx q[0];
rz(2.5600938) q[0];
x q[1];
rz(-0.99230687) q[2];
sx q[2];
rz(-1.6725146) q[2];
sx q[2];
rz(-2.5755239) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1350007) q[1];
sx q[1];
rz(-1.1524832) q[1];
sx q[1];
rz(2.5043284) q[1];
x q[2];
rz(-0.067906109) q[3];
sx q[3];
rz(-1.0722245) q[3];
sx q[3];
rz(-0.25085051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.095801) q[2];
sx q[2];
rz(-1.2597151) q[2];
sx q[2];
rz(3.1328787) q[2];
rz(-1.2039315) q[3];
sx q[3];
rz(-2.7976076) q[3];
sx q[3];
rz(-3.1194527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1190204) q[0];
sx q[0];
rz(-0.38563269) q[0];
sx q[0];
rz(-2.2517396) q[0];
rz(-1.285137) q[1];
sx q[1];
rz(-1.0533918) q[1];
sx q[1];
rz(-0.13592517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0014627731) q[0];
sx q[0];
rz(-2.9076932) q[0];
sx q[0];
rz(-2.9359096) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0338342) q[2];
sx q[2];
rz(-1.0091678) q[2];
sx q[2];
rz(1.6800576) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7116075) q[1];
sx q[1];
rz(-1.2384336) q[1];
sx q[1];
rz(1.7794987) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9985524) q[3];
sx q[3];
rz(-1.2238127) q[3];
sx q[3];
rz(0.59511371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4385684) q[2];
sx q[2];
rz(-0.88006222) q[2];
sx q[2];
rz(1.0620037) q[2];
rz(0.92987531) q[3];
sx q[3];
rz(-2.3558741) q[3];
sx q[3];
rz(0.78833956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7787665) q[0];
sx q[0];
rz(-2.7500948) q[0];
sx q[0];
rz(-2.7407001) q[0];
rz(2.3800384) q[1];
sx q[1];
rz(-1.4925894) q[1];
sx q[1];
rz(-2.3525995) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74247201) q[0];
sx q[0];
rz(-1.6982887) q[0];
sx q[0];
rz(-1.5546868) q[0];
x q[1];
rz(-2.9862271) q[2];
sx q[2];
rz(-2.2399358) q[2];
sx q[2];
rz(-2.5071627) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4915062) q[1];
sx q[1];
rz(-1.6900926) q[1];
sx q[1];
rz(-2.5507798) q[1];
rz(-1.5212359) q[3];
sx q[3];
rz(-1.0465099) q[3];
sx q[3];
rz(-1.4567284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5679428) q[2];
sx q[2];
rz(-1.2929792) q[2];
sx q[2];
rz(-0.72081494) q[2];
rz(1.4912262) q[3];
sx q[3];
rz(-2.3590915) q[3];
sx q[3];
rz(1.3097552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1184621) q[0];
sx q[0];
rz(-0.11496249) q[0];
sx q[0];
rz(-2.3638828) q[0];
rz(-0.65981162) q[1];
sx q[1];
rz(-2.2751364) q[1];
sx q[1];
rz(-0.39001098) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4068425) q[0];
sx q[0];
rz(-1.5336541) q[0];
sx q[0];
rz(-1.4358621) q[0];
rz(2.397178) q[2];
sx q[2];
rz(-1.0420024) q[2];
sx q[2];
rz(-2.2802558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.240814) q[1];
sx q[1];
rz(-1.1692746) q[1];
sx q[1];
rz(-2.8177849) q[1];
rz(-2.1999851) q[3];
sx q[3];
rz(-2.3857085) q[3];
sx q[3];
rz(2.1516277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4593792) q[2];
sx q[2];
rz(-0.43872681) q[2];
sx q[2];
rz(2.6518346) q[2];
rz(-0.30719906) q[3];
sx q[3];
rz(-2.2178853) q[3];
sx q[3];
rz(1.3908305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3601892) q[0];
sx q[0];
rz(-1.6506945) q[0];
sx q[0];
rz(-1.6682464) q[0];
rz(2.0669943) q[1];
sx q[1];
rz(-0.85292024) q[1];
sx q[1];
rz(-1.9180752) q[1];
rz(0.39911453) q[2];
sx q[2];
rz(-1.21798) q[2];
sx q[2];
rz(-1.1026364) q[2];
rz(-2.1297602) q[3];
sx q[3];
rz(-1.0444416) q[3];
sx q[3];
rz(-2.039644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
