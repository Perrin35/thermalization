OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(1.4899878) q[0];
sx q[0];
rz(8.4943354) q[0];
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(1.1073444) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9288899) q[0];
sx q[0];
rz(-1.2455997) q[0];
sx q[0];
rz(0.70744275) q[0];
rz(-2.1030175) q[2];
sx q[2];
rz(-2.3468809) q[2];
sx q[2];
rz(2.6796535) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.63065527) q[1];
sx q[1];
rz(-1.8488548) q[1];
sx q[1];
rz(3.0250938) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7928365) q[3];
sx q[3];
rz(-0.96491279) q[3];
sx q[3];
rz(-2.6989258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.063623108) q[2];
sx q[2];
rz(-2.412553) q[2];
sx q[2];
rz(1.3280274) q[2];
rz(-2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(-0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6593453) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(2.2609718) q[0];
rz(1.2940787) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(2.3243288) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6200703) q[0];
sx q[0];
rz(-1.9154473) q[0];
sx q[0];
rz(-0.14981139) q[0];
x q[1];
rz(0.21821071) q[2];
sx q[2];
rz(-2.0514224) q[2];
sx q[2];
rz(-2.6618119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.17727795) q[1];
sx q[1];
rz(-0.69523584) q[1];
sx q[1];
rz(-2.7362105) q[1];
x q[2];
rz(-1.0975295) q[3];
sx q[3];
rz(-1.739199) q[3];
sx q[3];
rz(-0.67697224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91784224) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(-2.7775653) q[2];
rz(-0.98637995) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-0.36689511) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(2.1752775) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(1.6945217) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9054748) q[0];
sx q[0];
rz(-0.21572278) q[0];
sx q[0];
rz(-3.059379) q[0];
x q[1];
rz(1.2800531) q[2];
sx q[2];
rz(-1.650052) q[2];
sx q[2];
rz(2.5153164) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6009532) q[1];
sx q[1];
rz(-2.1689231) q[1];
sx q[1];
rz(-1.1126493) q[1];
rz(-2.8964642) q[3];
sx q[3];
rz(-1.9760625) q[3];
sx q[3];
rz(0.95642904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43859279) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(-1.1052216) q[2];
rz(-2.3953719) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(2.1658649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352017) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(1.4659457) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-2.7526061) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3272414) q[0];
sx q[0];
rz(-1.8530493) q[0];
sx q[0];
rz(1.809657) q[0];
rz(-pi) q[1];
rz(0.19548266) q[2];
sx q[2];
rz(-2.0494378) q[2];
sx q[2];
rz(-3.0340956) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71967857) q[1];
sx q[1];
rz(-0.72243566) q[1];
sx q[1];
rz(0.6694442) q[1];
rz(-2.9205434) q[3];
sx q[3];
rz(-2.3813558) q[3];
sx q[3];
rz(1.5628712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-2.6848865) q[2];
rz(-1.6263973) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(0.28234282) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221955) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(-0.36002457) q[0];
rz(0.64741627) q[1];
sx q[1];
rz(-1.5292239) q[1];
sx q[1];
rz(2.6470851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9279328) q[0];
sx q[0];
rz(-1.4466009) q[0];
sx q[0];
rz(2.859982) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4353384) q[2];
sx q[2];
rz(-1.6098621) q[2];
sx q[2];
rz(-1.0328229) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9322885) q[1];
sx q[1];
rz(-2.5087025) q[1];
sx q[1];
rz(-3.0401405) q[1];
rz(-pi) q[2];
rz(2.6974929) q[3];
sx q[3];
rz(-1.5937623) q[3];
sx q[3];
rz(2.3820153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4313724) q[2];
sx q[2];
rz(-2.4857095) q[2];
sx q[2];
rz(2.8430856) q[2];
rz(-2.8295637) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(2.503094) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9962149) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(2.9456855) q[0];
rz(-0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.9063937) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0876448) q[0];
sx q[0];
rz(-1.2637648) q[0];
sx q[0];
rz(-0.35015492) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7487095) q[2];
sx q[2];
rz(-2.3977931) q[2];
sx q[2];
rz(1.5768029) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0022618) q[1];
sx q[1];
rz(-1.1552703) q[1];
sx q[1];
rz(-2.5050487) q[1];
rz(2.9052832) q[3];
sx q[3];
rz(-1.6522539) q[3];
sx q[3];
rz(2.8636275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49332508) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(-2.7098999) q[2];
rz(1.4124983) q[3];
sx q[3];
rz(-2.4192211) q[3];
sx q[3];
rz(-3.0055962) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7504808) q[0];
sx q[0];
rz(-1.1688122) q[0];
sx q[0];
rz(3.0294763) q[0];
rz(-0.21513367) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(-2.0281866) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33681413) q[0];
sx q[0];
rz(-1.830173) q[0];
sx q[0];
rz(-0.39774261) q[0];
x q[1];
rz(2.8261975) q[2];
sx q[2];
rz(-1.3239667) q[2];
sx q[2];
rz(-2.18404) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50892936) q[1];
sx q[1];
rz(-0.76970657) q[1];
sx q[1];
rz(-0.44429227) q[1];
rz(-pi) q[2];
rz(1.9951622) q[3];
sx q[3];
rz(-0.48454912) q[3];
sx q[3];
rz(-0.68819203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1404184) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(0.12602885) q[2];
rz(-2.0942988) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(2.4333911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.4814608) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(1.2449645) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(1.9326899) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93223244) q[0];
sx q[0];
rz(-0.27420843) q[0];
sx q[0];
rz(1.2742395) q[0];
rz(-pi) q[1];
rz(2.4822794) q[2];
sx q[2];
rz(-2.1954143) q[2];
sx q[2];
rz(-2.5059932) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.066684494) q[1];
sx q[1];
rz(-1.9004596) q[1];
sx q[1];
rz(2.3782905) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2743837) q[3];
sx q[3];
rz(-0.65215462) q[3];
sx q[3];
rz(1.6782827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7631491) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(1.3195999) q[2];
rz(2.54946) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(-0.035141703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9086583) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(-1.7804902) q[0];
rz(1.9305485) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(0.39168721) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93145934) q[0];
sx q[0];
rz(-1.1362695) q[0];
sx q[0];
rz(-1.5772485) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7462256) q[2];
sx q[2];
rz(-2.6569416) q[2];
sx q[2];
rz(-1.6389099) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0605676) q[1];
sx q[1];
rz(-2.276366) q[1];
sx q[1];
rz(2.702436) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6158995) q[3];
sx q[3];
rz(-1.5075397) q[3];
sx q[3];
rz(-0.69378187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6212375) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(-1.3367782) q[2];
rz(2.3796066) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8005463) q[0];
sx q[0];
rz(-0.30277345) q[0];
sx q[0];
rz(0.57089943) q[0];
rz(1.4292498) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(-0.16194078) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83308342) q[0];
sx q[0];
rz(-2.1662795) q[0];
sx q[0];
rz(2.7916629) q[0];
rz(0.0091153092) q[2];
sx q[2];
rz(-1.3822123) q[2];
sx q[2];
rz(-1.3537223) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3028114) q[1];
sx q[1];
rz(-1.4421842) q[1];
sx q[1];
rz(-1.7437115) q[1];
rz(-pi) q[2];
rz(-1.6041683) q[3];
sx q[3];
rz(-1.0479234) q[3];
sx q[3];
rz(0.066730412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1841715) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(1.7133678) q[2];
rz(-1.1994294) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(-1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4409055) q[0];
sx q[0];
rz(-2.6661243) q[0];
sx q[0];
rz(2.0211924) q[0];
rz(-1.3700925) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(3.0225171) q[2];
sx q[2];
rz(-1.4243813) q[2];
sx q[2];
rz(1.4204155) q[2];
rz(0.45392848) q[3];
sx q[3];
rz(-0.47149999) q[3];
sx q[3];
rz(-0.75054689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
