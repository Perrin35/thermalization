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
rz(-1.4077185) q[0];
sx q[0];
rz(-0.044535927) q[0];
rz(-1.1107923) q[1];
sx q[1];
rz(-1.2244616) q[1];
sx q[1];
rz(-0.78723025) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7101332) q[0];
sx q[0];
rz(-2.0636807) q[0];
sx q[0];
rz(-2.5566275) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8960612) q[2];
sx q[2];
rz(-1.6084387) q[2];
sx q[2];
rz(-1.4311439) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7234955) q[1];
sx q[1];
rz(-2.5586794) q[1];
sx q[1];
rz(1.6969796) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98712943) q[3];
sx q[3];
rz(-0.41582022) q[3];
sx q[3];
rz(0.73279954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.74701509) q[2];
sx q[2];
rz(-2.0614823) q[2];
sx q[2];
rz(2.3874095) q[2];
rz(-1.7265823) q[3];
sx q[3];
rz(-2.6455046) q[3];
sx q[3];
rz(2.8126341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-1.0873782) q[0];
sx q[0];
rz(-0.35254756) q[0];
sx q[0];
rz(-1.3341599) q[0];
rz(1.2913903) q[1];
sx q[1];
rz(-1.038237) q[1];
sx q[1];
rz(0.87563595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8311616) q[0];
sx q[0];
rz(-1.0410032) q[0];
sx q[0];
rz(2.3710662) q[0];
x q[1];
rz(-0.062280999) q[2];
sx q[2];
rz(-2.0665633) q[2];
sx q[2];
rz(-0.77046662) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.39159472) q[1];
sx q[1];
rz(-1.1340967) q[1];
sx q[1];
rz(-2.8207247) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8505136) q[3];
sx q[3];
rz(-1.3716231) q[3];
sx q[3];
rz(1.8739623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69677532) q[2];
sx q[2];
rz(-2.4397218) q[2];
sx q[2];
rz(-2.2689421) q[2];
rz(0.78898346) q[3];
sx q[3];
rz(-2.240447) q[3];
sx q[3];
rz(-0.22855973) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92662421) q[0];
sx q[0];
rz(-1.7890395) q[0];
sx q[0];
rz(3.1335926) q[0];
rz(0.81832111) q[1];
sx q[1];
rz(-0.74940959) q[1];
sx q[1];
rz(-1.9047033) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0465572) q[0];
sx q[0];
rz(-0.76991428) q[0];
sx q[0];
rz(-1.1004992) q[0];
rz(2.7929467) q[2];
sx q[2];
rz(-1.9816035) q[2];
sx q[2];
rz(-0.88627671) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8244918) q[1];
sx q[1];
rz(-0.6380837) q[1];
sx q[1];
rz(-0.53148766) q[1];
rz(-pi) q[2];
rz(1.6881006) q[3];
sx q[3];
rz(-1.3356195) q[3];
sx q[3];
rz(1.7010353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41039738) q[2];
sx q[2];
rz(-3.0072913) q[2];
sx q[2];
rz(-0.94089874) q[2];
rz(2.0375552) q[3];
sx q[3];
rz(-1.3528115) q[3];
sx q[3];
rz(1.1436536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0333772) q[0];
sx q[0];
rz(-0.5431076) q[0];
sx q[0];
rz(-0.23750842) q[0];
rz(-0.6595276) q[1];
sx q[1];
rz(-0.57412761) q[1];
sx q[1];
rz(-0.017017078) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.983611) q[0];
sx q[0];
rz(-1.5226304) q[0];
sx q[0];
rz(1.401813) q[0];
x q[1];
rz(-2.723316) q[2];
sx q[2];
rz(-0.47992963) q[2];
sx q[2];
rz(-2.0945702) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5524872) q[1];
sx q[1];
rz(-2.1334339) q[1];
sx q[1];
rz(1.4472812) q[1];
rz(-0.83570273) q[3];
sx q[3];
rz(-0.91360211) q[3];
sx q[3];
rz(3.0007255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90618769) q[2];
sx q[2];
rz(-1.9402639) q[2];
sx q[2];
rz(-2.2271633) q[2];
rz(0.36214456) q[3];
sx q[3];
rz(-2.3320964) q[3];
sx q[3];
rz(-2.5563498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.888716) q[0];
sx q[0];
rz(-1.7676366) q[0];
sx q[0];
rz(-0.1121029) q[0];
rz(-1.0653227) q[1];
sx q[1];
rz(-2.0707097) q[1];
sx q[1];
rz(1.0692474) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0651848) q[0];
sx q[0];
rz(-1.2994081) q[0];
sx q[0];
rz(2.8057764) q[0];
x q[1];
rz(-2.4222071) q[2];
sx q[2];
rz(-2.4771538) q[2];
sx q[2];
rz(1.4321355) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4878975) q[1];
sx q[1];
rz(-1.2733165) q[1];
sx q[1];
rz(1.2384227) q[1];
x q[2];
rz(-2.1700056) q[3];
sx q[3];
rz(-1.6808482) q[3];
sx q[3];
rz(2.9436265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4273044) q[2];
sx q[2];
rz(-2.4236743) q[2];
sx q[2];
rz(0.53421268) q[2];
rz(-0.30321768) q[3];
sx q[3];
rz(-1.5568308) q[3];
sx q[3];
rz(1.6259954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41820207) q[0];
sx q[0];
rz(-1.284282) q[0];
sx q[0];
rz(-2.2499625) q[0];
rz(-1.3806237) q[1];
sx q[1];
rz(-1.9386407) q[1];
sx q[1];
rz(0.92323971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5017446) q[0];
sx q[0];
rz(-3.0303133) q[0];
sx q[0];
rz(-2.6529045) q[0];
rz(2.7545287) q[2];
sx q[2];
rz(-1.3157857) q[2];
sx q[2];
rz(-0.86266359) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90369836) q[1];
sx q[1];
rz(-0.35939068) q[1];
sx q[1];
rz(2.2713347) q[1];
x q[2];
rz(-0.51839101) q[3];
sx q[3];
rz(-0.71037358) q[3];
sx q[3];
rz(-0.47391665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0455857) q[2];
sx q[2];
rz(-0.43232375) q[2];
sx q[2];
rz(-1.879479) q[2];
rz(1.9803842) q[3];
sx q[3];
rz(-1.5491756) q[3];
sx q[3];
rz(1.2459292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4534509) q[0];
sx q[0];
rz(-2.3095486) q[0];
sx q[0];
rz(-1.6081109) q[0];
rz(-0.43117943) q[1];
sx q[1];
rz(-0.9175514) q[1];
sx q[1];
rz(-1.7327488) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4154177) q[0];
sx q[0];
rz(-1.2841932) q[0];
sx q[0];
rz(0.64985259) q[0];
rz(-pi) q[1];
rz(3.0529019) q[2];
sx q[2];
rz(-1.9541395) q[2];
sx q[2];
rz(0.92244043) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.085693749) q[1];
sx q[1];
rz(-1.6961462) q[1];
sx q[1];
rz(2.9841444) q[1];
rz(-0.124229) q[3];
sx q[3];
rz(-1.0719187) q[3];
sx q[3];
rz(1.4829829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93671736) q[2];
sx q[2];
rz(-2.016158) q[2];
sx q[2];
rz(-1.3196866) q[2];
rz(0.034505757) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(1.3694793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7337604) q[0];
sx q[0];
rz(-0.5540846) q[0];
sx q[0];
rz(-1.2169417) q[0];
rz(0.92562428) q[1];
sx q[1];
rz(-1.3652912) q[1];
sx q[1];
rz(-1.9409174) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.901163) q[0];
sx q[0];
rz(-1.0679185) q[0];
sx q[0];
rz(-0.85977413) q[0];
rz(-0.56270069) q[2];
sx q[2];
rz(-0.68018736) q[2];
sx q[2];
rz(-2.0923751) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4958805) q[1];
sx q[1];
rz(-0.79949841) q[1];
sx q[1];
rz(-1.6281566) q[1];
x q[2];
rz(0.36603277) q[3];
sx q[3];
rz(-1.1732374) q[3];
sx q[3];
rz(1.1432858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.919148) q[2];
sx q[2];
rz(-1.6067952) q[2];
sx q[2];
rz(-1.0486802) q[2];
rz(0.18516304) q[3];
sx q[3];
rz(-1.1566297) q[3];
sx q[3];
rz(-3.0345501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.01934) q[0];
sx q[0];
rz(-0.66146079) q[0];
sx q[0];
rz(-0.81073236) q[0];
rz(-2.4108389) q[1];
sx q[1];
rz(-2.0850875) q[1];
sx q[1];
rz(-2.1810541) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98314032) q[0];
sx q[0];
rz(-3.1380655) q[0];
sx q[0];
rz(2.8615224) q[0];
x q[1];
rz(2.075718) q[2];
sx q[2];
rz(-1.9532579) q[2];
sx q[2];
rz(1.7590211) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3523622) q[1];
sx q[1];
rz(-1.8327729) q[1];
sx q[1];
rz(2.1195223) q[1];
x q[2];
rz(2.8836684) q[3];
sx q[3];
rz(-1.1561794) q[3];
sx q[3];
rz(2.1900511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1824128) q[2];
sx q[2];
rz(-1.7824087) q[2];
sx q[2];
rz(-1.5394999) q[2];
rz(-1.0473853) q[3];
sx q[3];
rz(-1.7644707) q[3];
sx q[3];
rz(-1.2654977) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78028107) q[0];
sx q[0];
rz(-0.86532101) q[0];
sx q[0];
rz(-2.5541232) q[0];
rz(-2.7777708) q[1];
sx q[1];
rz(-1.9963341) q[1];
sx q[1];
rz(-0.90715539) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22940561) q[0];
sx q[0];
rz(-0.36021458) q[0];
sx q[0];
rz(-0.2304669) q[0];
rz(-1.9586708) q[2];
sx q[2];
rz(-2.0128314) q[2];
sx q[2];
rz(1.3988053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7228667) q[1];
sx q[1];
rz(-1.360824) q[1];
sx q[1];
rz(-0.46011113) q[1];
rz(1.9349707) q[3];
sx q[3];
rz(-1.6505906) q[3];
sx q[3];
rz(1.4353115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5017447) q[2];
sx q[2];
rz(-0.088968337) q[2];
sx q[2];
rz(-0.42195827) q[2];
rz(-3.1320599) q[3];
sx q[3];
rz(-1.5671174) q[3];
sx q[3];
rz(2.6154521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52787732) q[0];
sx q[0];
rz(-1.8462702) q[0];
sx q[0];
rz(0.12311002) q[0];
rz(2.4404424) q[1];
sx q[1];
rz(-1.2260561) q[1];
sx q[1];
rz(-1.4969926) q[1];
rz(-2.4451679) q[2];
sx q[2];
rz(-1.8581894) q[2];
sx q[2];
rz(-2.8530864) q[2];
rz(-0.27746874) q[3];
sx q[3];
rz(-1.1744842) q[3];
sx q[3];
rz(-0.13406772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
