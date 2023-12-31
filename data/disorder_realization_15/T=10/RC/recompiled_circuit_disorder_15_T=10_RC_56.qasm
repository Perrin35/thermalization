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
rz(-1.7742046) q[1];
sx q[1];
rz(-2.8957638) q[1];
sx q[1];
rz(-2.153102) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6160625) q[0];
sx q[0];
rz(-1.5126192) q[0];
sx q[0];
rz(1.9921897) q[0];
rz(-0.97889401) q[2];
sx q[2];
rz(-1.4375293) q[2];
sx q[2];
rz(2.8507581) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.58952744) q[1];
sx q[1];
rz(-1.1067992) q[1];
sx q[1];
rz(0.5281756) q[1];
rz(-3.1036166) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(1.9248885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6926379) q[2];
sx q[2];
rz(-1.8779034) q[2];
sx q[2];
rz(-2.5732178) q[2];
rz(1.189399) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(-2.9590759) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73873591) q[0];
sx q[0];
rz(-1.0915272) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(-0.96827132) q[1];
sx q[1];
rz(-2.4666511) q[1];
sx q[1];
rz(-1.8889069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2371212) q[0];
sx q[0];
rz(-1.7303109) q[0];
sx q[0];
rz(-1.6617387) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86979903) q[2];
sx q[2];
rz(-0.79248488) q[2];
sx q[2];
rz(3.062641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15236552) q[1];
sx q[1];
rz(-0.75132912) q[1];
sx q[1];
rz(-0.82528021) q[1];
x q[2];
rz(2.0086803) q[3];
sx q[3];
rz(-1.4012194) q[3];
sx q[3];
rz(1.358043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0210555) q[2];
sx q[2];
rz(-1.3201822) q[2];
sx q[2];
rz(2.7978314) q[2];
rz(-0.3343285) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(-0.46935579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6702328) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(-0.0070455889) q[0];
rz(0.37653157) q[1];
sx q[1];
rz(-2.2129602) q[1];
sx q[1];
rz(0.71281707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8757652) q[0];
sx q[0];
rz(-1.7808502) q[0];
sx q[0];
rz(1.34583) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9489281) q[2];
sx q[2];
rz(-1.4107804) q[2];
sx q[2];
rz(-0.68607054) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.30046001) q[1];
sx q[1];
rz(-1.3679879) q[1];
sx q[1];
rz(-2.1122785) q[1];
x q[2];
rz(-3.0611638) q[3];
sx q[3];
rz(-2.0518528) q[3];
sx q[3];
rz(-1.2642494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.787848) q[2];
sx q[2];
rz(-1.5793844) q[2];
sx q[2];
rz(-0.66398579) q[2];
rz(-1.239423) q[3];
sx q[3];
rz(-2.911471) q[3];
sx q[3];
rz(-0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5749213) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(-2.5675039) q[0];
rz(-2.8494049) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(-2.2132197) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5260122) q[0];
sx q[0];
rz(-0.86475879) q[0];
sx q[0];
rz(-3.1403149) q[0];
rz(0.33889126) q[2];
sx q[2];
rz(-2.3232984) q[2];
sx q[2];
rz(-1.3202867) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2239383) q[1];
sx q[1];
rz(-1.3242711) q[1];
sx q[1];
rz(-1.1397821) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5023605) q[3];
sx q[3];
rz(-2.1661048) q[3];
sx q[3];
rz(2.3140698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0488247) q[2];
sx q[2];
rz(-2.6269045) q[2];
sx q[2];
rz(-1.1553923) q[2];
rz(2.998735) q[3];
sx q[3];
rz(-2.6070049) q[3];
sx q[3];
rz(2.3891881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0578385) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(0.074247867) q[0];
rz(-1.6429365) q[1];
sx q[1];
rz(-1.5131806) q[1];
sx q[1];
rz(-2.1898851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1002016) q[0];
sx q[0];
rz(-0.75267422) q[0];
sx q[0];
rz(1.4647096) q[0];
rz(-pi) q[1];
rz(0.98624595) q[2];
sx q[2];
rz(-2.2918321) q[2];
sx q[2];
rz(-0.69794387) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.73478991) q[1];
sx q[1];
rz(-2.4996335) q[1];
sx q[1];
rz(-1.2924679) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0025025) q[3];
sx q[3];
rz(-1.9145402) q[3];
sx q[3];
rz(2.7078201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4962861) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(-1.9892233) q[2];
rz(-2.0955829) q[3];
sx q[3];
rz(-2.4029616) q[3];
sx q[3];
rz(0.0013105198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1482658) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(-2.6519725) q[0];
rz(-2.1173677) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(1.5354059) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5912356) q[0];
sx q[0];
rz(-3.0315657) q[0];
sx q[0];
rz(-0.45022924) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68131955) q[2];
sx q[2];
rz(-2.1207003) q[2];
sx q[2];
rz(2.6256109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8294576) q[1];
sx q[1];
rz(-1.8541938) q[1];
sx q[1];
rz(2.9764922) q[1];
rz(-pi) q[2];
rz(-2.9165886) q[3];
sx q[3];
rz(-2.5000754) q[3];
sx q[3];
rz(-2.0819506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9635222) q[2];
sx q[2];
rz(-0.81729752) q[2];
sx q[2];
rz(-2.9658588) q[2];
rz(1.0774311) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(3.1350465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734633) q[0];
sx q[0];
rz(-0.21502762) q[0];
sx q[0];
rz(1.2699132) q[0];
rz(2.9636256) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(1.7756745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2451536) q[0];
sx q[0];
rz(-1.23587) q[0];
sx q[0];
rz(-0.73672898) q[0];
rz(1.3283417) q[2];
sx q[2];
rz(-1.2521267) q[2];
sx q[2];
rz(0.53675011) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7682225) q[1];
sx q[1];
rz(-1.3942413) q[1];
sx q[1];
rz(-2.8819487) q[1];
rz(-pi) q[2];
rz(-2.6122983) q[3];
sx q[3];
rz(-1.4636453) q[3];
sx q[3];
rz(1.2308434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0718677) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(-2.243637) q[2];
rz(-0.14363025) q[3];
sx q[3];
rz(-1.8840021) q[3];
sx q[3];
rz(-2.7974421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9550069) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(-2.8198077) q[0];
rz(-2.2161662) q[1];
sx q[1];
rz(-1.7089475) q[1];
sx q[1];
rz(0.61703533) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077438146) q[0];
sx q[0];
rz(-1.7700717) q[0];
sx q[0];
rz(-1.827924) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7210984) q[2];
sx q[2];
rz(-0.21810025) q[2];
sx q[2];
rz(0.66767603) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7223282) q[1];
sx q[1];
rz(-1.1280578) q[1];
sx q[1];
rz(-0.41132136) q[1];
rz(-0.048181941) q[3];
sx q[3];
rz(-0.74996862) q[3];
sx q[3];
rz(-0.90660209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59163219) q[2];
sx q[2];
rz(-2.154921) q[2];
sx q[2];
rz(0.11403306) q[2];
rz(2.7791801) q[3];
sx q[3];
rz(-2.896538) q[3];
sx q[3];
rz(-1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47700259) q[0];
sx q[0];
rz(-2.6053071) q[0];
sx q[0];
rz(-1.7846918) q[0];
rz(-0.82018954) q[1];
sx q[1];
rz(-2.7892022) q[1];
sx q[1];
rz(1.483451) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8720855) q[0];
sx q[0];
rz(-1.5604696) q[0];
sx q[0];
rz(1.4875862) q[0];
rz(-pi) q[1];
rz(1.6805029) q[2];
sx q[2];
rz(-0.98329558) q[2];
sx q[2];
rz(-0.40949958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7328651) q[1];
sx q[1];
rz(-0.36204007) q[1];
sx q[1];
rz(-1.15508) q[1];
x q[2];
rz(0.58247329) q[3];
sx q[3];
rz(-0.63342047) q[3];
sx q[3];
rz(0.42788351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1411529) q[2];
sx q[2];
rz(-0.65594643) q[2];
sx q[2];
rz(-1.8971987) q[2];
rz(0.42516285) q[3];
sx q[3];
rz(-0.77912283) q[3];
sx q[3];
rz(-1.5104793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.72220951) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(-0.37049946) q[0];
rz(-2.0314979) q[1];
sx q[1];
rz(-2.4659174) q[1];
sx q[1];
rz(-1.8006905) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0144893) q[0];
sx q[0];
rz(-1.2047486) q[0];
sx q[0];
rz(-0.24665477) q[0];
rz(-pi) q[1];
rz(-1.4341899) q[2];
sx q[2];
rz(-2.280526) q[2];
sx q[2];
rz(0.030985706) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9677876) q[1];
sx q[1];
rz(-2.6952744) q[1];
sx q[1];
rz(0.96149573) q[1];
x q[2];
rz(0.16593905) q[3];
sx q[3];
rz(-3.0186742) q[3];
sx q[3];
rz(-0.77792203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.6293388) q[2];
sx q[2];
rz(-0.24015716) q[2];
sx q[2];
rz(1.6188999) q[2];
rz(-0.86087888) q[3];
sx q[3];
rz(-1.4168134) q[3];
sx q[3];
rz(-0.035877429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9027949) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(-2.8181656) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(1.9771489) q[2];
sx q[2];
rz(-1.5928369) q[2];
sx q[2];
rz(-1.9441446) q[2];
rz(1.0988416) q[3];
sx q[3];
rz(-1.9190211) q[3];
sx q[3];
rz(2.2091051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
