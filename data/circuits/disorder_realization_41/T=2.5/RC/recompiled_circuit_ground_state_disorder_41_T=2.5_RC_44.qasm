OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1908258) q[0];
sx q[0];
rz(-1.289225) q[0];
sx q[0];
rz(3.0486795) q[0];
rz(3.0463123) q[1];
sx q[1];
rz(-2.4089101) q[1];
sx q[1];
rz(-1.9021775) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5372972) q[0];
sx q[0];
rz(-1.3193466) q[0];
sx q[0];
rz(-0.23668134) q[0];
x q[1];
rz(-2.5074717) q[2];
sx q[2];
rz(-2.8180052) q[2];
sx q[2];
rz(0.35558082) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.69529205) q[1];
sx q[1];
rz(-1.9009034) q[1];
sx q[1];
rz(-1.5395201) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69783437) q[3];
sx q[3];
rz(-2.4437332) q[3];
sx q[3];
rz(0.35240155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4702845) q[2];
sx q[2];
rz(-1.4802063) q[2];
sx q[2];
rz(-1.7926463) q[2];
rz(0.74364439) q[3];
sx q[3];
rz(-2.8642004) q[3];
sx q[3];
rz(0.17717895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.141356) q[0];
sx q[0];
rz(-1.5970705) q[0];
sx q[0];
rz(-2.8479688) q[0];
rz(-1.0307505) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(-2.7986599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9259277) q[0];
sx q[0];
rz(-2.5747888) q[0];
sx q[0];
rz(1.180992) q[0];
rz(-1.198771) q[2];
sx q[2];
rz(-1.4497533) q[2];
sx q[2];
rz(0.72265676) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.98951116) q[1];
sx q[1];
rz(-0.73484269) q[1];
sx q[1];
rz(2.8192725) q[1];
x q[2];
rz(1.2323805) q[3];
sx q[3];
rz(-0.65753257) q[3];
sx q[3];
rz(-2.7595208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.12884101) q[2];
sx q[2];
rz(-1.3920471) q[2];
sx q[2];
rz(0.7592321) q[2];
rz(-0.020708474) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(-2.7032963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6193806) q[0];
sx q[0];
rz(-0.601957) q[0];
sx q[0];
rz(-0.0044599175) q[0];
rz(-0.64747539) q[1];
sx q[1];
rz(-0.59440333) q[1];
sx q[1];
rz(2.3562145) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9849094) q[0];
sx q[0];
rz(-1.7298609) q[0];
sx q[0];
rz(-1.6522264) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9663229) q[2];
sx q[2];
rz(-1.2614377) q[2];
sx q[2];
rz(2.4966405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33357802) q[1];
sx q[1];
rz(-1.6393264) q[1];
sx q[1];
rz(-1.0295111) q[1];
rz(-pi) q[2];
rz(-2.249751) q[3];
sx q[3];
rz(-2.1168609) q[3];
sx q[3];
rz(-0.31072703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.37041) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(1.7837589) q[2];
rz(0.34058288) q[3];
sx q[3];
rz(-0.95652306) q[3];
sx q[3];
rz(0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99712813) q[0];
sx q[0];
rz(-1.0839394) q[0];
sx q[0];
rz(-0.88515627) q[0];
rz(-0.74686933) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(1.0964099) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80712705) q[0];
sx q[0];
rz(-1.2594481) q[0];
sx q[0];
rz(1.2047355) q[0];
x q[1];
rz(2.0463767) q[2];
sx q[2];
rz(-3.0747483) q[2];
sx q[2];
rz(0.16193709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.29624048) q[1];
sx q[1];
rz(-1.8295604) q[1];
sx q[1];
rz(0.47611632) q[1];
x q[2];
rz(-0.73170264) q[3];
sx q[3];
rz(-1.1771132) q[3];
sx q[3];
rz(-1.7997774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5235644) q[2];
sx q[2];
rz(-0.21054331) q[2];
sx q[2];
rz(3.1169685) q[2];
rz(-2.5060182) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(-0.88328254) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4516975) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(2.6746993) q[0];
rz(-0.11058552) q[1];
sx q[1];
rz(-0.46375912) q[1];
sx q[1];
rz(1.0708403) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8798053) q[0];
sx q[0];
rz(-2.6222485) q[0];
sx q[0];
rz(-2.065361) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51009615) q[2];
sx q[2];
rz(-0.3613216) q[2];
sx q[2];
rz(0.23621836) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.22510281) q[1];
sx q[1];
rz(-1.1943294) q[1];
sx q[1];
rz(2.589499) q[1];
rz(-pi) q[2];
rz(2.4399906) q[3];
sx q[3];
rz(-1.0402586) q[3];
sx q[3];
rz(0.70487937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0232627) q[2];
sx q[2];
rz(-1.3619962) q[2];
sx q[2];
rz(-2.2892717) q[2];
rz(0.037633745) q[3];
sx q[3];
rz(-1.9084385) q[3];
sx q[3];
rz(0.5948624) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1325876) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(1.2888541) q[0];
rz(1.6291078) q[1];
sx q[1];
rz(-1.1237203) q[1];
sx q[1];
rz(1.1423133) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61296755) q[0];
sx q[0];
rz(-0.36699793) q[0];
sx q[0];
rz(-2.8662843) q[0];
x q[1];
rz(-1.283848) q[2];
sx q[2];
rz(-0.21636886) q[2];
sx q[2];
rz(-1.8277825) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9276322) q[1];
sx q[1];
rz(-2.7134656) q[1];
sx q[1];
rz(0.63251782) q[1];
x q[2];
rz(0.31123881) q[3];
sx q[3];
rz(-1.4025619) q[3];
sx q[3];
rz(-0.5420891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8196572) q[2];
sx q[2];
rz(-1.8737996) q[2];
sx q[2];
rz(0.12040559) q[2];
rz(1.1558007) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(-0.22404484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33893809) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(-1.6876203) q[0];
rz(1.5441719) q[1];
sx q[1];
rz(-1.495196) q[1];
sx q[1];
rz(2.6905751) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7190349) q[0];
sx q[0];
rz(-1.8717878) q[0];
sx q[0];
rz(1.6331571) q[0];
rz(-1.0181997) q[2];
sx q[2];
rz(-0.42358735) q[2];
sx q[2];
rz(0.46679631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4828795) q[1];
sx q[1];
rz(-2.6900243) q[1];
sx q[1];
rz(2.8612479) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5805832) q[3];
sx q[3];
rz(-2.3811445) q[3];
sx q[3];
rz(-1.1475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94577998) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(-1.3607402) q[2];
rz(-1.0055379) q[3];
sx q[3];
rz(-1.1755627) q[3];
sx q[3];
rz(-1.3895234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176395) q[0];
sx q[0];
rz(-0.12549505) q[0];
sx q[0];
rz(-2.4355167) q[0];
rz(-1.5438682) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(-2.2854038) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5498915) q[0];
sx q[0];
rz(-3.1208) q[0];
sx q[0];
rz(-1.9338305) q[0];
rz(-pi) q[1];
rz(-0.82243408) q[2];
sx q[2];
rz(-1.8959799) q[2];
sx q[2];
rz(2.6564244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6996307) q[1];
sx q[1];
rz(-0.61304997) q[1];
sx q[1];
rz(2.7837201) q[1];
x q[2];
rz(1.3230349) q[3];
sx q[3];
rz(-0.50163402) q[3];
sx q[3];
rz(2.786123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84264821) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(-0.16743463) q[2];
rz(1.5873448) q[3];
sx q[3];
rz(-0.7730248) q[3];
sx q[3];
rz(-2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3779959) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(2.7698621) q[0];
rz(1.4338214) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(-0.30002123) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0473012) q[0];
sx q[0];
rz(-2.8326748) q[0];
sx q[0];
rz(-2.7500217) q[0];
rz(3.07934) q[2];
sx q[2];
rz(-2.5635898) q[2];
sx q[2];
rz(2.5684772) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2310861) q[1];
sx q[1];
rz(-0.66201895) q[1];
sx q[1];
rz(2.9390017) q[1];
x q[2];
rz(-2.2737502) q[3];
sx q[3];
rz(-1.4757089) q[3];
sx q[3];
rz(2.8114708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6649449) q[2];
sx q[2];
rz(-0.46755329) q[2];
sx q[2];
rz(2.7139968) q[2];
rz(-1.556373) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(-0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4204191) q[0];
sx q[0];
rz(-2.6009646) q[0];
sx q[0];
rz(-0.46947259) q[0];
rz(0.92533127) q[1];
sx q[1];
rz(-1.0877437) q[1];
sx q[1];
rz(0.023177711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3566475) q[0];
sx q[0];
rz(-2.0933505) q[0];
sx q[0];
rz(-1.3808668) q[0];
rz(-pi) q[1];
rz(-0.59892861) q[2];
sx q[2];
rz(-0.99236503) q[2];
sx q[2];
rz(0.79858649) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7851023) q[1];
sx q[1];
rz(-1.7447326) q[1];
sx q[1];
rz(2.9746303) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7816824) q[3];
sx q[3];
rz(-0.74112219) q[3];
sx q[3];
rz(-1.9884584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2716486) q[2];
sx q[2];
rz(-1.2056489) q[2];
sx q[2];
rz(2.6507157) q[2];
rz(2.0513746) q[3];
sx q[3];
rz(-2.1722983) q[3];
sx q[3];
rz(0.31392613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8563817) q[0];
sx q[0];
rz(-2.2638392) q[0];
sx q[0];
rz(1.0765156) q[0];
rz(-2.8648227) q[1];
sx q[1];
rz(-2.3634187) q[1];
sx q[1];
rz(1.707911) q[1];
rz(-1.2645012) q[2];
sx q[2];
rz(-0.90521348) q[2];
sx q[2];
rz(0.95059849) q[2];
rz(-1.5908949) q[3];
sx q[3];
rz(-2.1051959) q[3];
sx q[3];
rz(0.013230562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
