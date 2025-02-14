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
rz(1.1433831) q[0];
sx q[0];
rz(-1.5902061) q[0];
sx q[0];
rz(-2.8144612) q[0];
rz(1.7786572) q[1];
sx q[1];
rz(-2.5442446) q[1];
sx q[1];
rz(-0.30244952) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0344447) q[0];
sx q[0];
rz(-1.4054789) q[0];
sx q[0];
rz(-1.3470696) q[0];
rz(-pi) q[1];
rz(-2.2313868) q[2];
sx q[2];
rz(-1.1288092) q[2];
sx q[2];
rz(3.0750753) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0781344) q[1];
sx q[1];
rz(-1.7094104) q[1];
sx q[1];
rz(1.8891508) q[1];
rz(-0.47306319) q[3];
sx q[3];
rz(-1.6636208) q[3];
sx q[3];
rz(-1.0582093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3737619) q[2];
sx q[2];
rz(-0.094002873) q[2];
sx q[2];
rz(-2.8685699) q[2];
rz(0.36429575) q[3];
sx q[3];
rz(-1.7184869) q[3];
sx q[3];
rz(-0.020615904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6881707) q[0];
sx q[0];
rz(-2.0437129) q[0];
sx q[0];
rz(1.4349487) q[0];
rz(2.9393328) q[1];
sx q[1];
rz(-2.5113998) q[1];
sx q[1];
rz(-0.30466255) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6759535) q[0];
sx q[0];
rz(-2.3871564) q[0];
sx q[0];
rz(-1.937768) q[0];
rz(-0.96370187) q[2];
sx q[2];
rz(-0.4956565) q[2];
sx q[2];
rz(-0.021481422) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.92485904) q[1];
sx q[1];
rz(-1.2302515) q[1];
sx q[1];
rz(1.7664643) q[1];
rz(2.3000003) q[3];
sx q[3];
rz(-2.9805956) q[3];
sx q[3];
rz(2.5567233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8649586) q[2];
sx q[2];
rz(-2.7101597) q[2];
sx q[2];
rz(-1.5383833) q[2];
rz(2.3254584) q[3];
sx q[3];
rz(-2.3147801) q[3];
sx q[3];
rz(-0.89474595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64033878) q[0];
sx q[0];
rz(-2.8479939) q[0];
sx q[0];
rz(-1.6295992) q[0];
rz(-1.8005796) q[1];
sx q[1];
rz(-0.87313849) q[1];
sx q[1];
rz(-2.8626633) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19353774) q[0];
sx q[0];
rz(-0.16866651) q[0];
sx q[0];
rz(-2.4381766) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90143335) q[2];
sx q[2];
rz(-2.3877904) q[2];
sx q[2];
rz(1.19687) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2532369) q[1];
sx q[1];
rz(-0.86730343) q[1];
sx q[1];
rz(2.7363267) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5666714) q[3];
sx q[3];
rz(-0.71625159) q[3];
sx q[3];
rz(-0.28820693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5035847) q[2];
sx q[2];
rz(-1.4853442) q[2];
sx q[2];
rz(-0.33835641) q[2];
rz(1.7342957) q[3];
sx q[3];
rz(-1.1275848) q[3];
sx q[3];
rz(-2.1710763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19891837) q[0];
sx q[0];
rz(-3.119454) q[0];
sx q[0];
rz(2.4778147) q[0];
rz(-2.5860419) q[1];
sx q[1];
rz(-1.5063565) q[1];
sx q[1];
rz(-1.7864071) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595918) q[0];
sx q[0];
rz(-2.2170236) q[0];
sx q[0];
rz(-1.5311956) q[0];
rz(-pi) q[1];
x q[1];
rz(0.023472114) q[2];
sx q[2];
rz(-1.9292574) q[2];
sx q[2];
rz(-1.1662999) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.042204) q[1];
sx q[1];
rz(-0.70800938) q[1];
sx q[1];
rz(-0.96570496) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.066800634) q[3];
sx q[3];
rz(-1.1054174) q[3];
sx q[3];
rz(-2.3540879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3163471) q[2];
sx q[2];
rz(-2.1630042) q[2];
sx q[2];
rz(-0.52581954) q[2];
rz(-0.60872269) q[3];
sx q[3];
rz(-2.5059097) q[3];
sx q[3];
rz(-0.70613247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.71259251) q[0];
sx q[0];
rz(-1.9569995) q[0];
sx q[0];
rz(-2.1767966) q[0];
rz(0.4110128) q[1];
sx q[1];
rz(-2.1844468) q[1];
sx q[1];
rz(2.029665) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85025609) q[0];
sx q[0];
rz(-1.5865069) q[0];
sx q[0];
rz(-1.5474623) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70198029) q[2];
sx q[2];
rz(-1.9147083) q[2];
sx q[2];
rz(-1.3683343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7644123) q[1];
sx q[1];
rz(-2.0736338) q[1];
sx q[1];
rz(-0.92288121) q[1];
rz(-pi) q[2];
rz(0.49647496) q[3];
sx q[3];
rz(-0.30204812) q[3];
sx q[3];
rz(-0.99629096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8370221) q[2];
sx q[2];
rz(-1.8111753) q[2];
sx q[2];
rz(-1.042761) q[2];
rz(1.1834831) q[3];
sx q[3];
rz(-2.2730998) q[3];
sx q[3];
rz(-0.97064251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17827621) q[0];
sx q[0];
rz(-1.2238418) q[0];
sx q[0];
rz(-2.8455612) q[0];
rz(3.1405247) q[1];
sx q[1];
rz(-1.3384534) q[1];
sx q[1];
rz(1.0329049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72275513) q[0];
sx q[0];
rz(-1.5745409) q[0];
sx q[0];
rz(-3.1360955) q[0];
rz(-pi) q[1];
rz(3.0216713) q[2];
sx q[2];
rz(-1.9397221) q[2];
sx q[2];
rz(-3.1110087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.008625) q[1];
sx q[1];
rz(-1.6264249) q[1];
sx q[1];
rz(2.8939496) q[1];
x q[2];
rz(-1.9719129) q[3];
sx q[3];
rz(-0.84904957) q[3];
sx q[3];
rz(-2.7665319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.380015) q[2];
sx q[2];
rz(-1.6753316) q[2];
sx q[2];
rz(-1.5029017) q[2];
rz(2.5401529) q[3];
sx q[3];
rz(-2.1419339) q[3];
sx q[3];
rz(-0.39826605) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10694207) q[0];
sx q[0];
rz(-1.5489738) q[0];
sx q[0];
rz(2.4079127) q[0];
rz(0.94379464) q[1];
sx q[1];
rz(-1.5422041) q[1];
sx q[1];
rz(2.5234047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9551463) q[0];
sx q[0];
rz(-1.6206564) q[0];
sx q[0];
rz(1.3912892) q[0];
rz(-1.5812065) q[2];
sx q[2];
rz(-0.73042831) q[2];
sx q[2];
rz(1.324162) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.427278) q[1];
sx q[1];
rz(-2.4763417) q[1];
sx q[1];
rz(1.6131496) q[1];
x q[2];
rz(-0.47241601) q[3];
sx q[3];
rz(-0.69184408) q[3];
sx q[3];
rz(-1.7418246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1980629) q[2];
sx q[2];
rz(-1.1897831) q[2];
sx q[2];
rz(2.5243536) q[2];
rz(1.9700358) q[3];
sx q[3];
rz(-0.37552437) q[3];
sx q[3];
rz(2.530781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010855762) q[0];
sx q[0];
rz(-2.3541088) q[0];
sx q[0];
rz(0.47295824) q[0];
rz(-1.6821776) q[1];
sx q[1];
rz(-1.4143896) q[1];
sx q[1];
rz(0.032729538) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41389909) q[0];
sx q[0];
rz(-1.5581456) q[0];
sx q[0];
rz(0.032300579) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47686994) q[2];
sx q[2];
rz(-0.82399659) q[2];
sx q[2];
rz(-1.0007953) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1188726) q[1];
sx q[1];
rz(-1.1274844) q[1];
sx q[1];
rz(2.9037649) q[1];
x q[2];
rz(-3.1113137) q[3];
sx q[3];
rz(-1.4631728) q[3];
sx q[3];
rz(-0.054747907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1590283) q[2];
sx q[2];
rz(-0.50243598) q[2];
sx q[2];
rz(1.095613) q[2];
rz(-0.76121965) q[3];
sx q[3];
rz(-1.8571564) q[3];
sx q[3];
rz(-2.718954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6765167) q[0];
sx q[0];
rz(-3.0618771) q[0];
sx q[0];
rz(-2.4327143) q[0];
rz(0.31157663) q[1];
sx q[1];
rz(-1.0919002) q[1];
sx q[1];
rz(0.30837217) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.209359) q[0];
sx q[0];
rz(-1.7361802) q[0];
sx q[0];
rz(2.7424521) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4089085) q[2];
sx q[2];
rz(-0.53712979) q[2];
sx q[2];
rz(-0.86974517) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.93058298) q[1];
sx q[1];
rz(-2.4786754) q[1];
sx q[1];
rz(-0.40056132) q[1];
rz(-pi) q[2];
rz(-2.3544054) q[3];
sx q[3];
rz(-1.7554253) q[3];
sx q[3];
rz(-0.61239374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.10820216) q[2];
sx q[2];
rz(-1.4806925) q[2];
sx q[2];
rz(2.9208276) q[2];
rz(1.0734142) q[3];
sx q[3];
rz(-2.150841) q[3];
sx q[3];
rz(2.3555135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49755001) q[0];
sx q[0];
rz(-2.2798517) q[0];
sx q[0];
rz(-2.1047237) q[0];
rz(-1.6715624) q[1];
sx q[1];
rz(-1.0187047) q[1];
sx q[1];
rz(1.1933614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84813335) q[0];
sx q[0];
rz(-2.6984697) q[0];
sx q[0];
rz(3.0713017) q[0];
rz(2.4429351) q[2];
sx q[2];
rz(-1.9113052) q[2];
sx q[2];
rz(-0.27891544) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9876839) q[1];
sx q[1];
rz(-1.7584956) q[1];
sx q[1];
rz(0.69324268) q[1];
rz(0.12659215) q[3];
sx q[3];
rz(-1.8014354) q[3];
sx q[3];
rz(-2.0960208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6679473) q[2];
sx q[2];
rz(-1.8569943) q[2];
sx q[2];
rz(0.13772193) q[2];
rz(1.600945) q[3];
sx q[3];
rz(-2.9100304) q[3];
sx q[3];
rz(0.93129492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4731428) q[0];
sx q[0];
rz(-1.4746329) q[0];
sx q[0];
rz(2.0869577) q[0];
rz(2.6558381) q[1];
sx q[1];
rz(-1.2243441) q[1];
sx q[1];
rz(-1.4996554) q[1];
rz(-0.68436868) q[2];
sx q[2];
rz(-1.3791313) q[2];
sx q[2];
rz(-2.287938) q[2];
rz(-2.6089062) q[3];
sx q[3];
rz(-0.83704994) q[3];
sx q[3];
rz(-0.34216135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
