OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34831369) q[0];
sx q[0];
rz(8.0170595) q[0];
sx q[0];
rz(9.4693139) q[0];
rz(2.0308004) q[1];
sx q[1];
rz(-1.9171311) q[1];
sx q[1];
rz(0.78723025) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6985748) q[0];
sx q[0];
rz(-2.0788143) q[0];
sx q[0];
rz(-2.1430911) q[0];
x q[1];
rz(-1.2455315) q[2];
sx q[2];
rz(-1.533154) q[2];
sx q[2];
rz(1.4311439) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2582142) q[1];
sx q[1];
rz(-1.6401263) q[1];
sx q[1];
rz(-0.99154559) q[1];
x q[2];
rz(0.23870339) q[3];
sx q[3];
rz(-1.2269964) q[3];
sx q[3];
rz(1.3581004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74701509) q[2];
sx q[2];
rz(-1.0801103) q[2];
sx q[2];
rz(2.3874095) q[2];
rz(1.4150103) q[3];
sx q[3];
rz(-2.6455046) q[3];
sx q[3];
rz(2.8126341) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0873782) q[0];
sx q[0];
rz(-0.35254756) q[0];
sx q[0];
rz(1.3341599) q[0];
rz(1.2913903) q[1];
sx q[1];
rz(-2.1033557) q[1];
sx q[1];
rz(2.2659567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.401225) q[0];
sx q[0];
rz(-2.2384907) q[0];
sx q[0];
rz(0.69913787) q[0];
rz(-1.0742168) q[2];
sx q[2];
rz(-1.5160217) q[2];
sx q[2];
rz(-0.77067256) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8671682) q[1];
sx q[1];
rz(-2.6058785) q[1];
sx q[1];
rz(-2.1650326) q[1];
rz(-pi) q[2];
rz(0.61296566) q[3];
sx q[3];
rz(-0.35109441) q[3];
sx q[3];
rz(0.28051091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4448173) q[2];
sx q[2];
rz(-0.70187086) q[2];
sx q[2];
rz(-2.2689421) q[2];
rz(2.3526092) q[3];
sx q[3];
rz(-2.240447) q[3];
sx q[3];
rz(-2.9130329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2149684) q[0];
sx q[0];
rz(-1.7890395) q[0];
sx q[0];
rz(0.0080000814) q[0];
rz(2.3232715) q[1];
sx q[1];
rz(-2.3921831) q[1];
sx q[1];
rz(-1.9047033) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.711126) q[0];
sx q[0];
rz(-2.240181) q[0];
sx q[0];
rz(-0.41394625) q[0];
rz(2.0047999) q[2];
sx q[2];
rz(-1.889359) q[2];
sx q[2];
rz(0.54036507) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8244918) q[1];
sx q[1];
rz(-2.503509) q[1];
sx q[1];
rz(0.53148766) q[1];
rz(-pi) q[2];
x q[2];
rz(1.453492) q[3];
sx q[3];
rz(-1.8059732) q[3];
sx q[3];
rz(-1.4405574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7311953) q[2];
sx q[2];
rz(-3.0072913) q[2];
sx q[2];
rz(0.94089874) q[2];
rz(-1.1040374) q[3];
sx q[3];
rz(-1.7887812) q[3];
sx q[3];
rz(-1.1436536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10821548) q[0];
sx q[0];
rz(-0.5431076) q[0];
sx q[0];
rz(-2.9040842) q[0];
rz(0.6595276) q[1];
sx q[1];
rz(-2.567465) q[1];
sx q[1];
rz(-0.017017078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1579817) q[0];
sx q[0];
rz(-1.5226304) q[0];
sx q[0];
rz(1.7397797) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3624362) q[2];
sx q[2];
rz(-2.0063498) q[2];
sx q[2];
rz(-0.58247936) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5524872) q[1];
sx q[1];
rz(-1.0081588) q[1];
sx q[1];
rz(1.4472812) q[1];
x q[2];
rz(-2.3364725) q[3];
sx q[3];
rz(-1.0110572) q[3];
sx q[3];
rz(-0.92529682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.235405) q[2];
sx q[2];
rz(-1.9402639) q[2];
sx q[2];
rz(-0.91442937) q[2];
rz(-0.36214456) q[3];
sx q[3];
rz(-2.3320964) q[3];
sx q[3];
rz(2.5563498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.888716) q[0];
sx q[0];
rz(-1.7676366) q[0];
sx q[0];
rz(3.0294898) q[0];
rz(2.0762699) q[1];
sx q[1];
rz(-2.0707097) q[1];
sx q[1];
rz(-2.0723453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99226588) q[0];
sx q[0];
rz(-2.7130824) q[0];
sx q[0];
rz(0.70121588) q[0];
rz(-pi) q[1];
rz(-0.71938558) q[2];
sx q[2];
rz(-0.66443887) q[2];
sx q[2];
rz(-1.7094572) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.62110626) q[1];
sx q[1];
rz(-2.6992528) q[1];
sx q[1];
rz(0.81653313) q[1];
x q[2];
rz(0.97158708) q[3];
sx q[3];
rz(-1.6808482) q[3];
sx q[3];
rz(-0.19796619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4273044) q[2];
sx q[2];
rz(-0.71791831) q[2];
sx q[2];
rz(-0.53421268) q[2];
rz(-0.30321768) q[3];
sx q[3];
rz(-1.5568308) q[3];
sx q[3];
rz(1.6259954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7233906) q[0];
sx q[0];
rz(-1.284282) q[0];
sx q[0];
rz(2.2499625) q[0];
rz(-1.760969) q[1];
sx q[1];
rz(-1.2029519) q[1];
sx q[1];
rz(0.92323971) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5864201) q[0];
sx q[0];
rz(-1.5186383) q[0];
sx q[0];
rz(-3.0432493) q[0];
rz(-pi) q[1];
rz(-0.60439749) q[2];
sx q[2];
rz(-2.6816419) q[2];
sx q[2];
rz(-0.15397554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90369836) q[1];
sx q[1];
rz(-0.35939068) q[1];
sx q[1];
rz(0.87025799) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9736863) q[3];
sx q[3];
rz(-2.1729762) q[3];
sx q[3];
rz(2.9704579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.096006958) q[2];
sx q[2];
rz(-0.43232375) q[2];
sx q[2];
rz(-1.879479) q[2];
rz(-1.9803842) q[3];
sx q[3];
rz(-1.5491756) q[3];
sx q[3];
rz(1.8956634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4534509) q[0];
sx q[0];
rz(-2.3095486) q[0];
sx q[0];
rz(-1.6081109) q[0];
rz(0.43117943) q[1];
sx q[1];
rz(-2.2240413) q[1];
sx q[1];
rz(1.4088438) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36700248) q[0];
sx q[0];
rz(-0.95159114) q[0];
sx q[0];
rz(-1.9253233) q[0];
x q[1];
rz(-1.7869912) q[2];
sx q[2];
rz(-0.39297418) q[2];
sx q[2];
rz(-2.4525688) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0558989) q[1];
sx q[1];
rz(-1.4454464) q[1];
sx q[1];
rz(2.9841444) q[1];
rz(3.0173637) q[3];
sx q[3];
rz(-1.0719187) q[3];
sx q[3];
rz(1.4829829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93671736) q[2];
sx q[2];
rz(-2.016158) q[2];
sx q[2];
rz(-1.3196866) q[2];
rz(0.034505757) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(-1.7721133) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7337604) q[0];
sx q[0];
rz(-0.5540846) q[0];
sx q[0];
rz(1.2169417) q[0];
rz(-0.92562428) q[1];
sx q[1];
rz(-1.3652912) q[1];
sx q[1];
rz(-1.2006753) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1797829) q[0];
sx q[0];
rz(-2.2968074) q[0];
sx q[0];
rz(-2.2711193) q[0];
rz(-pi) q[1];
rz(1.9782135) q[2];
sx q[2];
rz(-1.0098741) q[2];
sx q[2];
rz(1.7307868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0266704) q[1];
sx q[1];
rz(-1.5296796) q[1];
sx q[1];
rz(2.369472) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86494495) q[3];
sx q[3];
rz(-2.6078913) q[3];
sx q[3];
rz(0.362901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.919148) q[2];
sx q[2];
rz(-1.5347975) q[2];
sx q[2];
rz(-2.0929125) q[2];
rz(-2.9564296) q[3];
sx q[3];
rz(-1.1566297) q[3];
sx q[3];
rz(-3.0345501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.01934) q[0];
sx q[0];
rz(-2.4801319) q[0];
sx q[0];
rz(0.81073236) q[0];
rz(-2.4108389) q[1];
sx q[1];
rz(-1.0565051) q[1];
sx q[1];
rz(-0.96053851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98314032) q[0];
sx q[0];
rz(-3.1380655) q[0];
sx q[0];
rz(2.8615224) q[0];
rz(-pi) q[1];
rz(-0.87709092) q[2];
sx q[2];
rz(-0.6232647) q[2];
sx q[2];
rz(0.40568144) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5170578) q[1];
sx q[1];
rz(-2.0987928) q[1];
sx q[1];
rz(0.30450423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0960363) q[3];
sx q[3];
rz(-0.48431319) q[3];
sx q[3];
rz(-2.7698295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9591799) q[2];
sx q[2];
rz(-1.7824087) q[2];
sx q[2];
rz(-1.5394999) q[2];
rz(-2.0942073) q[3];
sx q[3];
rz(-1.3771219) q[3];
sx q[3];
rz(1.876095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3613116) q[0];
sx q[0];
rz(-2.2762716) q[0];
sx q[0];
rz(-2.5541232) q[0];
rz(-2.7777708) q[1];
sx q[1];
rz(-1.1452585) q[1];
sx q[1];
rz(0.90715539) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22940561) q[0];
sx q[0];
rz(-0.36021458) q[0];
sx q[0];
rz(0.2304669) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1829218) q[2];
sx q[2];
rz(-2.0128314) q[2];
sx q[2];
rz(-1.7427874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75395155) q[1];
sx q[1];
rz(-0.50261897) q[1];
sx q[1];
rz(2.6941264) q[1];
rz(0.085368319) q[3];
sx q[3];
rz(-1.2078346) q[3];
sx q[3];
rz(-0.1051108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.63984799) q[2];
sx q[2];
rz(-0.088968337) q[2];
sx q[2];
rz(-0.42195827) q[2];
rz(0.0095327775) q[3];
sx q[3];
rz(-1.5671174) q[3];
sx q[3];
rz(-0.52614051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137153) q[0];
sx q[0];
rz(-1.8462702) q[0];
sx q[0];
rz(0.12311002) q[0];
rz(0.70115024) q[1];
sx q[1];
rz(-1.9155365) q[1];
sx q[1];
rz(1.6446) q[1];
rz(-2.7098165) q[2];
sx q[2];
rz(-2.3974621) q[2];
sx q[2];
rz(1.5324788) q[2];
rz(1.1604068) q[3];
sx q[3];
rz(-1.8262564) q[3];
sx q[3];
rz(-1.5953596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
