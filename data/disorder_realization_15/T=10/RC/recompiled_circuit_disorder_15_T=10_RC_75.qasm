OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8437334) q[0];
sx q[0];
rz(-0.61367404) q[0];
sx q[0];
rz(0.71917978) q[0];
rz(1.367388) q[1];
sx q[1];
rz(-0.24582882) q[1];
sx q[1];
rz(2.153102) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0702695) q[0];
sx q[0];
rz(-1.9914314) q[0];
sx q[0];
rz(3.0778528) q[0];
rz(-0.16015957) q[2];
sx q[2];
rz(-0.9848435) q[2];
sx q[2];
rz(1.1908659) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4157384) q[1];
sx q[1];
rz(-2.038318) q[1];
sx q[1];
rz(1.0456677) q[1];
x q[2];
rz(-0.037976102) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(-1.9248885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6926379) q[2];
sx q[2];
rz(-1.8779034) q[2];
sx q[2];
rz(0.56837481) q[2];
rz(-1.189399) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(2.9590759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73873591) q[0];
sx q[0];
rz(-2.0500654) q[0];
sx q[0];
rz(-2.7785595) q[0];
rz(-2.1733213) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(-1.8889069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7583) q[0];
sx q[0];
rz(-0.18342605) q[0];
sx q[0];
rz(-2.6276876) q[0];
rz(-pi) q[1];
rz(2.5622796) q[2];
sx q[2];
rz(-0.99537731) q[2];
sx q[2];
rz(-0.79798165) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15236552) q[1];
sx q[1];
rz(-0.75132912) q[1];
sx q[1];
rz(-0.82528021) q[1];
rz(-0.18685347) q[3];
sx q[3];
rz(-1.1396176) q[3];
sx q[3];
rz(-0.29160515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12053717) q[2];
sx q[2];
rz(-1.3201822) q[2];
sx q[2];
rz(2.7978314) q[2];
rz(-0.3343285) q[3];
sx q[3];
rz(-2.7407586) q[3];
sx q[3];
rz(-2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47135982) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(3.1345471) q[0];
rz(0.37653157) q[1];
sx q[1];
rz(-0.9286325) q[1];
sx q[1];
rz(2.4287756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8757652) q[0];
sx q[0];
rz(-1.7808502) q[0];
sx q[0];
rz(1.7957627) q[0];
rz(-1.1926646) q[2];
sx q[2];
rz(-1.7308123) q[2];
sx q[2];
rz(0.68607054) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1944113) q[1];
sx q[1];
rz(-0.57465034) q[1];
sx q[1];
rz(1.9504207) q[1];
rz(-pi) q[2];
rz(1.4180693) q[3];
sx q[3];
rz(-2.6543791) q[3];
sx q[3];
rz(-1.0917851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3537447) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(-2.4776069) q[2];
rz(1.9021696) q[3];
sx q[3];
rz(-2.911471) q[3];
sx q[3];
rz(-0.91528875) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5666714) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(-2.5675039) q[0];
rz(-0.29218778) q[1];
sx q[1];
rz(-2.8657587) q[1];
sx q[1];
rz(-2.2132197) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5260122) q[0];
sx q[0];
rz(-0.86475879) q[0];
sx q[0];
rz(-3.1403149) q[0];
rz(-1.9119772) q[2];
sx q[2];
rz(-0.8114292) q[2];
sx q[2];
rz(0.84412837) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0005562) q[1];
sx q[1];
rz(-0.49266854) q[1];
sx q[1];
rz(1.0286742) q[1];
x q[2];
rz(-2.6392322) q[3];
sx q[3];
rz(-0.97548786) q[3];
sx q[3];
rz(0.82752284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.092768) q[2];
sx q[2];
rz(-2.6269045) q[2];
sx q[2];
rz(1.9862004) q[2];
rz(-2.998735) q[3];
sx q[3];
rz(-0.5345878) q[3];
sx q[3];
rz(2.3891881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0578385) q[0];
sx q[0];
rz(-0.46828073) q[0];
sx q[0];
rz(-0.074247867) q[0];
rz(1.6429365) q[1];
sx q[1];
rz(-1.5131806) q[1];
sx q[1];
rz(2.1898851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4518406) q[0];
sx q[0];
rz(-1.4983488) q[0];
sx q[0];
rz(-2.3206582) q[0];
rz(-0.98624595) q[2];
sx q[2];
rz(-2.2918321) q[2];
sx q[2];
rz(0.69794387) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4068027) q[1];
sx q[1];
rz(-0.64195913) q[1];
sx q[1];
rz(1.2924679) q[1];
rz(-2.1577155) q[3];
sx q[3];
rz(-0.65423274) q[3];
sx q[3];
rz(0.65140843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4962861) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(-1.9892233) q[2];
rz(-2.0955829) q[3];
sx q[3];
rz(-0.73863107) q[3];
sx q[3];
rz(3.1402821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.1482658) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(-2.6519725) q[0];
rz(2.1173677) q[1];
sx q[1];
rz(-0.63413292) q[1];
sx q[1];
rz(1.5354059) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4682966) q[0];
sx q[0];
rz(-1.522994) q[0];
sx q[0];
rz(0.099138069) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68131955) q[2];
sx q[2];
rz(-1.0208924) q[2];
sx q[2];
rz(0.51598179) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3052169) q[1];
sx q[1];
rz(-1.4123385) q[1];
sx q[1];
rz(1.8578908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9165886) q[3];
sx q[3];
rz(-0.64151728) q[3];
sx q[3];
rz(-2.0819506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9635222) q[2];
sx q[2];
rz(-0.81729752) q[2];
sx q[2];
rz(-2.9658588) q[2];
rz(-2.0641616) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(-0.0065461672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068129383) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(-1.8716795) q[0];
rz(-2.9636256) q[1];
sx q[1];
rz(-2.3033419) q[1];
sx q[1];
rz(1.7756745) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7570092) q[0];
sx q[0];
rz(-2.2582044) q[0];
sx q[0];
rz(-2.0100726) q[0];
rz(-1.3283417) q[2];
sx q[2];
rz(-1.2521267) q[2];
sx q[2];
rz(2.6048425) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.3733702) q[1];
sx q[1];
rz(-1.7473514) q[1];
sx q[1];
rz(-2.8819487) q[1];
rz(-pi) q[2];
rz(2.9317022) q[3];
sx q[3];
rz(-0.53901796) q[3];
sx q[3];
rz(0.15912661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0697249) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(0.89795566) q[2];
rz(0.14363025) q[3];
sx q[3];
rz(-1.8840021) q[3];
sx q[3];
rz(2.7974421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1865858) q[0];
sx q[0];
rz(-2.2877559) q[0];
sx q[0];
rz(-0.32178497) q[0];
rz(0.92542648) q[1];
sx q[1];
rz(-1.7089475) q[1];
sx q[1];
rz(0.61703533) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0641545) q[0];
sx q[0];
rz(-1.371521) q[0];
sx q[0];
rz(-1.827924) q[0];
x q[1];
rz(-3.1084193) q[2];
sx q[2];
rz(-1.355194) q[2];
sx q[2];
rz(2.6278091) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.033213381) q[1];
sx q[1];
rz(-1.9404267) q[1];
sx q[1];
rz(-1.0934248) q[1];
rz(-pi) q[2];
rz(-1.5259605) q[3];
sx q[3];
rz(-2.3196844) q[3];
sx q[3];
rz(-2.3007948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59163219) q[2];
sx q[2];
rz(-2.154921) q[2];
sx q[2];
rz(-3.0275596) q[2];
rz(-2.7791801) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(1.9053649) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6645901) q[0];
sx q[0];
rz(-2.6053071) q[0];
sx q[0];
rz(-1.7846918) q[0];
rz(-0.82018954) q[1];
sx q[1];
rz(-2.7892022) q[1];
sx q[1];
rz(1.483451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8720855) q[0];
sx q[0];
rz(-1.5604696) q[0];
sx q[0];
rz(-1.6540065) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16295095) q[2];
sx q[2];
rz(-0.59646791) q[2];
sx q[2];
rz(-0.60566723) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7328651) q[1];
sx q[1];
rz(-2.7795526) q[1];
sx q[1];
rz(-1.15508) q[1];
rz(-pi) q[2];
rz(0.55012196) q[3];
sx q[3];
rz(-1.9024444) q[3];
sx q[3];
rz(-1.5105997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1411529) q[2];
sx q[2];
rz(-0.65594643) q[2];
sx q[2];
rz(1.2443939) q[2];
rz(2.7164298) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(1.6311133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4193831) q[0];
sx q[0];
rz(-0.49383759) q[0];
sx q[0];
rz(-0.37049946) q[0];
rz(2.0314979) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(-1.8006905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140708) q[0];
sx q[0];
rz(-2.7033269) q[0];
sx q[0];
rz(-2.1379495) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9843763) q[2];
sx q[2];
rz(-2.4210857) q[2];
sx q[2];
rz(0.23888982) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6262628) q[1];
sx q[1];
rz(-1.9326107) q[1];
sx q[1];
rz(2.8742909) q[1];
x q[2];
rz(-2.9756536) q[3];
sx q[3];
rz(-3.0186742) q[3];
sx q[3];
rz(2.3636706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5122539) q[2];
sx q[2];
rz(-0.24015716) q[2];
sx q[2];
rz(-1.5226927) q[2];
rz(-0.86087888) q[3];
sx q[3];
rz(-1.4168134) q[3];
sx q[3];
rz(3.1057152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9027949) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(0.32342708) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(-1.515083) q[2];
sx q[2];
rz(-2.734676) q[2];
sx q[2];
rz(2.717072) q[2];
rz(-0.38701804) q[3];
sx q[3];
rz(-1.1292463) q[3];
sx q[3];
rz(-2.675727) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
