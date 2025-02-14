OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0082173) q[0];
sx q[0];
rz(-1.2674588) q[0];
sx q[0];
rz(-0.01292364) q[0];
rz(-2.456993) q[1];
sx q[1];
rz(-0.79889387) q[1];
sx q[1];
rz(-1.0577143) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2051281) q[0];
sx q[0];
rz(-1.3442486) q[0];
sx q[0];
rz(0.675987) q[0];
rz(-pi) q[1];
rz(0.77785141) q[2];
sx q[2];
rz(-1.9490644) q[2];
sx q[2];
rz(3.0562378) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5200708) q[1];
sx q[1];
rz(-2.3316262) q[1];
sx q[1];
rz(2.6317627) q[1];
x q[2];
rz(-1.4679883) q[3];
sx q[3];
rz(-2.7479726) q[3];
sx q[3];
rz(3.032544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.58618033) q[2];
sx q[2];
rz(-1.5988028) q[2];
sx q[2];
rz(0.093322873) q[2];
rz(3.1210476) q[3];
sx q[3];
rz(-2.8829657) q[3];
sx q[3];
rz(-1.7790022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697295) q[0];
sx q[0];
rz(-1.7427895) q[0];
sx q[0];
rz(-0.82759696) q[0];
rz(-0.96356511) q[1];
sx q[1];
rz(-1.5255442) q[1];
sx q[1];
rz(2.4049984) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44880906) q[0];
sx q[0];
rz(-2.8570456) q[0];
sx q[0];
rz(1.8440767) q[0];
rz(0.57666333) q[2];
sx q[2];
rz(-1.0500337) q[2];
sx q[2];
rz(1.7885838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1087649) q[1];
sx q[1];
rz(-1.2695644) q[1];
sx q[1];
rz(-0.21183045) q[1];
rz(0.16421825) q[3];
sx q[3];
rz(-3.0053557) q[3];
sx q[3];
rz(0.0020023684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57614148) q[2];
sx q[2];
rz(-0.57500035) q[2];
sx q[2];
rz(-0.74419332) q[2];
rz(0.60892504) q[3];
sx q[3];
rz(-2.3600793) q[3];
sx q[3];
rz(-1.6412546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610483) q[0];
sx q[0];
rz(-0.86549509) q[0];
sx q[0];
rz(1.3737099) q[0];
rz(-2.3420077) q[1];
sx q[1];
rz(-2.1345963) q[1];
sx q[1];
rz(-2.0094357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4904259) q[0];
sx q[0];
rz(-1.4803783) q[0];
sx q[0];
rz(1.7866485) q[0];
x q[1];
rz(2.8707809) q[2];
sx q[2];
rz(-0.083148227) q[2];
sx q[2];
rz(-1.6390334) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3721322) q[1];
sx q[1];
rz(-0.67382183) q[1];
sx q[1];
rz(-0.12427434) q[1];
rz(-pi) q[2];
rz(-1.6145094) q[3];
sx q[3];
rz(-2.1792886) q[3];
sx q[3];
rz(-0.19245806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75480294) q[2];
sx q[2];
rz(-2.5567882) q[2];
sx q[2];
rz(-1.7542138) q[2];
rz(-2.7925708) q[3];
sx q[3];
rz(-1.6882221) q[3];
sx q[3];
rz(-2.6121228) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.741852) q[0];
sx q[0];
rz(-2.1735503) q[0];
sx q[0];
rz(2.7834748) q[0];
rz(1.0707567) q[1];
sx q[1];
rz(-2.4929969) q[1];
sx q[1];
rz(1.7020285) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62240206) q[0];
sx q[0];
rz(-2.3163412) q[0];
sx q[0];
rz(-2.78016) q[0];
x q[1];
rz(-2.1522572) q[2];
sx q[2];
rz(-2.3609567) q[2];
sx q[2];
rz(2.6065741) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.44612158) q[1];
sx q[1];
rz(-2.2955986) q[1];
sx q[1];
rz(-2.3430969) q[1];
rz(2.374745) q[3];
sx q[3];
rz(-1.7550751) q[3];
sx q[3];
rz(-0.16367463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7122571) q[2];
sx q[2];
rz(-1.9108994) q[2];
sx q[2];
rz(-2.5168391) q[2];
rz(-1.5077) q[3];
sx q[3];
rz(-2.3816536) q[3];
sx q[3];
rz(-1.1225351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37609142) q[0];
sx q[0];
rz(-2.2608345) q[0];
sx q[0];
rz(1.1464024) q[0];
rz(1.435185) q[1];
sx q[1];
rz(-0.51992661) q[1];
sx q[1];
rz(-0.82040876) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4389184) q[0];
sx q[0];
rz(-1.6587388) q[0];
sx q[0];
rz(1.8431208) q[0];
rz(-1.7960288) q[2];
sx q[2];
rz(-2.4625255) q[2];
sx q[2];
rz(-2.1175543) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8294404) q[1];
sx q[1];
rz(-1.7746801) q[1];
sx q[1];
rz(1.3711434) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6040398) q[3];
sx q[3];
rz(-0.98131949) q[3];
sx q[3];
rz(-1.4840028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72941214) q[2];
sx q[2];
rz(-2.086144) q[2];
sx q[2];
rz(-2.6030276) q[2];
rz(-1.6719079) q[3];
sx q[3];
rz(-1.6939751) q[3];
sx q[3];
rz(-3.0830429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84022123) q[0];
sx q[0];
rz(-3.0026307) q[0];
sx q[0];
rz(-3.050991) q[0];
rz(0.36965707) q[1];
sx q[1];
rz(-1.548111) q[1];
sx q[1];
rz(-2.7634117) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9210631) q[0];
sx q[0];
rz(-0.53595966) q[0];
sx q[0];
rz(2.6820763) q[0];
x q[1];
rz(3.0514977) q[2];
sx q[2];
rz(-1.6125154) q[2];
sx q[2];
rz(2.1509474) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.52247574) q[1];
sx q[1];
rz(-1.4475087) q[1];
sx q[1];
rz(-2.6465764) q[1];
rz(-pi) q[2];
x q[2];
rz(1.380081) q[3];
sx q[3];
rz(-1.6487953) q[3];
sx q[3];
rz(-1.6140779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52046627) q[2];
sx q[2];
rz(-1.6660606) q[2];
sx q[2];
rz(0.10144083) q[2];
rz(1.8488047) q[3];
sx q[3];
rz(-0.83270508) q[3];
sx q[3];
rz(-1.4147991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8349649) q[0];
sx q[0];
rz(-1.2766301) q[0];
sx q[0];
rz(-2.2391338) q[0];
rz(-0.2392256) q[1];
sx q[1];
rz(-1.3419515) q[1];
sx q[1];
rz(-0.36144027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1162565) q[0];
sx q[0];
rz(-2.6317208) q[0];
sx q[0];
rz(-1.340853) q[0];
rz(-pi) q[1];
rz(1.5148411) q[2];
sx q[2];
rz(-1.6978886) q[2];
sx q[2];
rz(-0.54886078) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9779404) q[1];
sx q[1];
rz(-1.1457902) q[1];
sx q[1];
rz(-0.24042701) q[1];
rz(-pi) q[2];
rz(-1.4124198) q[3];
sx q[3];
rz(-1.90622) q[3];
sx q[3];
rz(0.80663825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.07301894) q[2];
sx q[2];
rz(-1.34015) q[2];
sx q[2];
rz(-0.87693357) q[2];
rz(-2.1584611) q[3];
sx q[3];
rz(-1.5777595) q[3];
sx q[3];
rz(0.94299281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8125732) q[0];
sx q[0];
rz(-2.946377) q[0];
sx q[0];
rz(-0.83874291) q[0];
rz(-2.4328531) q[1];
sx q[1];
rz(-2.510431) q[1];
sx q[1];
rz(0.90604025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7278465) q[0];
sx q[0];
rz(-2.4165396) q[0];
sx q[0];
rz(-0.031812761) q[0];
x q[1];
rz(-1.0308517) q[2];
sx q[2];
rz(-2.738852) q[2];
sx q[2];
rz(2.7583721) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3248277) q[1];
sx q[1];
rz(-2.205619) q[1];
sx q[1];
rz(-2.2544508) q[1];
rz(-2.8894807) q[3];
sx q[3];
rz(-2.1641755) q[3];
sx q[3];
rz(-0.93170792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16443843) q[2];
sx q[2];
rz(-0.70455569) q[2];
sx q[2];
rz(-1.1158811) q[2];
rz(-0.62853938) q[3];
sx q[3];
rz(-1.081859) q[3];
sx q[3];
rz(-2.5270497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81166613) q[0];
sx q[0];
rz(-2.3217432) q[0];
sx q[0];
rz(0.16127583) q[0];
rz(-2.2992086) q[1];
sx q[1];
rz(-1.7120275) q[1];
sx q[1];
rz(-0.17002034) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9901491) q[0];
sx q[0];
rz(-0.020016622) q[0];
sx q[0];
rz(2.722094) q[0];
x q[1];
rz(1.4030946) q[2];
sx q[2];
rz(-1.1360886) q[2];
sx q[2];
rz(-1.1230043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61115757) q[1];
sx q[1];
rz(-0.47947219) q[1];
sx q[1];
rz(-1.1392639) q[1];
rz(-pi) q[2];
rz(0.028548553) q[3];
sx q[3];
rz(-1.1880837) q[3];
sx q[3];
rz(3.1392424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9876447) q[2];
sx q[2];
rz(-1.597007) q[2];
sx q[2];
rz(1.4863996) q[2];
rz(-0.060128309) q[3];
sx q[3];
rz(-2.3190053) q[3];
sx q[3];
rz(-1.4415584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63115591) q[0];
sx q[0];
rz(-0.031351723) q[0];
sx q[0];
rz(1.4309058) q[0];
rz(0.36262861) q[1];
sx q[1];
rz(-1.5321833) q[1];
sx q[1];
rz(-1.3154202) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5820438) q[0];
sx q[0];
rz(-1.2799731) q[0];
sx q[0];
rz(1.2469588) q[0];
rz(-pi) q[1];
rz(-1.1926629) q[2];
sx q[2];
rz(-1.6844498) q[2];
sx q[2];
rz(-0.092157539) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9104256) q[1];
sx q[1];
rz(-2.1520414) q[1];
sx q[1];
rz(1.2509996) q[1];
rz(-pi) q[2];
rz(0.86851991) q[3];
sx q[3];
rz(-1.9316626) q[3];
sx q[3];
rz(-1.8250193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0843087) q[2];
sx q[2];
rz(-1.5949275) q[2];
sx q[2];
rz(-0.16051897) q[2];
rz(-0.48216835) q[3];
sx q[3];
rz(-2.4653258) q[3];
sx q[3];
rz(-0.67960656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1228444) q[0];
sx q[0];
rz(-1.7457122) q[0];
sx q[0];
rz(2.2298298) q[0];
rz(1.979076) q[1];
sx q[1];
rz(-1.5717506) q[1];
sx q[1];
rz(2.7350978) q[1];
rz(1.6655501) q[2];
sx q[2];
rz(-1.4838525) q[2];
sx q[2];
rz(-2.9529689) q[2];
rz(-1.0033458) q[3];
sx q[3];
rz(-0.88450817) q[3];
sx q[3];
rz(-1.3664288) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
