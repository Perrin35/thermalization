OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5443213) q[0];
sx q[0];
rz(-2.7379524) q[0];
sx q[0];
rz(0.37024745) q[0];
rz(-2.9397842) q[1];
sx q[1];
rz(-1.8887853) q[1];
sx q[1];
rz(-1.7226146) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0050632) q[0];
sx q[0];
rz(-0.29924527) q[0];
sx q[0];
rz(2.5314999) q[0];
rz(-pi) q[1];
rz(-0.51214829) q[2];
sx q[2];
rz(-1.3379659) q[2];
sx q[2];
rz(1.4621967) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4048684) q[1];
sx q[1];
rz(-2.1202592) q[1];
sx q[1];
rz(-3.0800372) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8098104) q[3];
sx q[3];
rz(-2.3506769) q[3];
sx q[3];
rz(3.0856109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3966763) q[2];
sx q[2];
rz(-1.4531206) q[2];
sx q[2];
rz(2.7837226) q[2];
rz(0.19168028) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1382004) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(-3.0733118) q[0];
rz(2.1349019) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(2.4904747) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2383645) q[0];
sx q[0];
rz(-2.4698386) q[0];
sx q[0];
rz(-0.3268468) q[0];
rz(-pi) q[1];
x q[1];
rz(3.100349) q[2];
sx q[2];
rz(-2.6154499) q[2];
sx q[2];
rz(-2.0522842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29775961) q[1];
sx q[1];
rz(-2.4503772) q[1];
sx q[1];
rz(0.83696951) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3531559) q[3];
sx q[3];
rz(-1.0895035) q[3];
sx q[3];
rz(2.7841115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6887001) q[2];
sx q[2];
rz(-0.15658997) q[2];
sx q[2];
rz(-0.90332705) q[2];
rz(2.3870758) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84905255) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(-0.49355155) q[0];
rz(-1.4986528) q[1];
sx q[1];
rz(-0.40619266) q[1];
sx q[1];
rz(1.0292056) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.878359) q[0];
sx q[0];
rz(-1.9089111) q[0];
sx q[0];
rz(2.6792206) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2649293) q[2];
sx q[2];
rz(-2.6138517) q[2];
sx q[2];
rz(1.9805679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5268847) q[1];
sx q[1];
rz(-1.7803368) q[1];
sx q[1];
rz(0.31063147) q[1];
rz(0.1251827) q[3];
sx q[3];
rz(-1.0297965) q[3];
sx q[3];
rz(2.560905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.8445245) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9349174) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(-2.7329965) q[0];
rz(1.4273377) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(0.88159195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6245759) q[0];
sx q[0];
rz(-1.7424184) q[0];
sx q[0];
rz(-2.4044328) q[0];
rz(-pi) q[1];
x q[1];
rz(2.300823) q[2];
sx q[2];
rz(-0.86849125) q[2];
sx q[2];
rz(-2.4384987) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3410014) q[1];
sx q[1];
rz(-2.2644271) q[1];
sx q[1];
rz(2.4800406) q[1];
x q[2];
rz(-1.6374171) q[3];
sx q[3];
rz(-1.1542218) q[3];
sx q[3];
rz(-1.6601603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0980229) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(1.1126474) q[2];
rz(2.5860795) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1495789) q[0];
sx q[0];
rz(-1.9861789) q[0];
sx q[0];
rz(1.3647112) q[0];
rz(0.31750202) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(-0.11725765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8506354) q[0];
sx q[0];
rz(-2.1521749) q[0];
sx q[0];
rz(1.4145538) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7848894) q[2];
sx q[2];
rz(-2.1630686) q[2];
sx q[2];
rz(1.5562197) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9980364) q[1];
sx q[1];
rz(-1.7245502) q[1];
sx q[1];
rz(-2.9424332) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4555143) q[3];
sx q[3];
rz(-2.5513463) q[3];
sx q[3];
rz(2.3595927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34510288) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.7618746) q[2];
rz(1.1896677) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(0.074507944) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9538486) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(-2.4353819) q[0];
rz(-2.0300991) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(3.0336753) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43150049) q[0];
sx q[0];
rz(-1.2762696) q[0];
sx q[0];
rz(-2.4157903) q[0];
rz(1.582167) q[2];
sx q[2];
rz(-1.6108542) q[2];
sx q[2];
rz(-0.96206059) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95277112) q[1];
sx q[1];
rz(-1.2880039) q[1];
sx q[1];
rz(-2.4464843) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48072731) q[3];
sx q[3];
rz(-0.93989621) q[3];
sx q[3];
rz(0.37068278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(2.0325913) q[2];
rz(1.2935982) q[3];
sx q[3];
rz(-1.3556017) q[3];
sx q[3];
rz(0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7727707) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(3.1266881) q[0];
rz(0.42124721) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(2.3419535) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0470456) q[0];
sx q[0];
rz(-2.4106815) q[0];
sx q[0];
rz(-0.38696179) q[0];
rz(-1.854935) q[2];
sx q[2];
rz(-1.26789) q[2];
sx q[2];
rz(0.45380935) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0529424) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(0.36693962) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1361254) q[3];
sx q[3];
rz(-1.6871095) q[3];
sx q[3];
rz(3.1073991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5853167) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(-0.9220534) q[2];
rz(1.8317892) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(-0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.17469445) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(-0.2640557) q[0];
rz(1.4008201) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(-2.82428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70586328) q[0];
sx q[0];
rz(-1.1840491) q[0];
sx q[0];
rz(-3.0630037) q[0];
x q[1];
rz(2.9810901) q[2];
sx q[2];
rz(-1.9151033) q[2];
sx q[2];
rz(-2.3366994) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4999769) q[1];
sx q[1];
rz(-2.9950905) q[1];
sx q[1];
rz(0.67540692) q[1];
rz(2.92225) q[3];
sx q[3];
rz(-2.4845124) q[3];
sx q[3];
rz(2.5430162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7060966) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(-1.502011) q[2];
rz(0.25029415) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(-1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.89649993) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(-1.7171575) q[0];
rz(0.4793438) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(-0.11553484) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53287017) q[0];
sx q[0];
rz(-1.0057797) q[0];
sx q[0];
rz(-0.67824058) q[0];
rz(-pi) q[1];
rz(-2.8250474) q[2];
sx q[2];
rz(-2.4015421) q[2];
sx q[2];
rz(-0.22063247) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.83976165) q[1];
sx q[1];
rz(-1.8592195) q[1];
sx q[1];
rz(2.0088058) q[1];
rz(2.4453836) q[3];
sx q[3];
rz(-1.6632102) q[3];
sx q[3];
rz(-0.19314167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4775548) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(-1.4257365) q[2];
rz(-1.4005631) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(-0.28276309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(-0.26915959) q[0];
rz(1.0956988) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(-1.7620618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36294286) q[0];
sx q[0];
rz(-1.4654219) q[0];
sx q[0];
rz(-1.4044579) q[0];
x q[1];
rz(-1.1534259) q[2];
sx q[2];
rz(-1.6438612) q[2];
sx q[2];
rz(-3.0550115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3456612) q[1];
sx q[1];
rz(-0.17594166) q[1];
sx q[1];
rz(-1.3962586) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74897154) q[3];
sx q[3];
rz(-2.7016771) q[3];
sx q[3];
rz(-0.16187748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0925838) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(1.6258378) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9243069) q[0];
sx q[0];
rz(-1.3047682) q[0];
sx q[0];
rz(-2.7346942) q[0];
rz(2.6869607) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(-0.17299962) q[2];
sx q[2];
rz(-2.0377918) q[2];
sx q[2];
rz(1.3755058) q[2];
rz(-0.061108246) q[3];
sx q[3];
rz(-1.8295049) q[3];
sx q[3];
rz(-2.6873333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
