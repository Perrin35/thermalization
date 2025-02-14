OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95076686) q[0];
sx q[0];
rz(-1.8523676) q[0];
sx q[0];
rz(0.092913203) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(2.4089101) q[1];
sx q[1];
rz(10.664193) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1151506) q[0];
sx q[0];
rz(-1.3416933) q[0];
sx q[0];
rz(-1.8291446) q[0];
x q[1];
rz(0.6341209) q[2];
sx q[2];
rz(-2.8180052) q[2];
sx q[2];
rz(0.35558082) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2762294) q[1];
sx q[1];
rz(-1.6003834) q[1];
sx q[1];
rz(-0.33025708) q[1];
rz(-pi) q[2];
rz(0.69783437) q[3];
sx q[3];
rz(-2.4437332) q[3];
sx q[3];
rz(2.7891911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6713082) q[2];
sx q[2];
rz(-1.6613864) q[2];
sx q[2];
rz(-1.3489464) q[2];
rz(-0.74364439) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(-2.9644137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0002366) q[0];
sx q[0];
rz(-1.5970705) q[0];
sx q[0];
rz(0.29362383) q[0];
rz(-1.0307505) q[1];
sx q[1];
rz(-1.3615969) q[1];
sx q[1];
rz(-0.34293276) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6688419) q[0];
sx q[0];
rz(-2.0905417) q[0];
sx q[0];
rz(2.9042871) q[0];
rz(-pi) q[1];
rz(3.0117576) q[2];
sx q[2];
rz(-1.9399683) q[2];
sx q[2];
rz(-2.2463727) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.575042) q[1];
sx q[1];
rz(-2.2600265) q[1];
sx q[1];
rz(1.291996) q[1];
rz(-pi) q[2];
rz(2.2003056) q[3];
sx q[3];
rz(-1.7751179) q[3];
sx q[3];
rz(2.2245537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0127516) q[2];
sx q[2];
rz(-1.3920471) q[2];
sx q[2];
rz(-2.3823605) q[2];
rz(-0.020708474) q[3];
sx q[3];
rz(-1.8755251) q[3];
sx q[3];
rz(-0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52221209) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(3.1371327) q[0];
rz(-0.64747539) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(-2.3562145) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9849094) q[0];
sx q[0];
rz(-1.4117318) q[0];
sx q[0];
rz(-1.6522264) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17526971) q[2];
sx q[2];
rz(-1.2614377) q[2];
sx q[2];
rz(2.4966405) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2783616) q[1];
sx q[1];
rz(-1.0309217) q[1];
sx q[1];
rz(-0.079915483) q[1];
x q[2];
rz(-2.249751) q[3];
sx q[3];
rz(-2.1168609) q[3];
sx q[3];
rz(2.8308656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7711827) q[2];
sx q[2];
rz(-0.9202756) q[2];
sx q[2];
rz(1.7837589) q[2];
rz(0.34058288) q[3];
sx q[3];
rz(-0.95652306) q[3];
sx q[3];
rz(-2.3497439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99712813) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(-2.2564364) q[0];
rz(2.3947233) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(-2.0451827) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0521419) q[0];
sx q[0];
rz(-0.47596395) q[0];
sx q[0];
rz(2.3030998) q[0];
rz(0.030641067) q[2];
sx q[2];
rz(-1.6302135) q[2];
sx q[2];
rz(-0.63842809) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9982352) q[1];
sx q[1];
rz(-2.0298185) q[1];
sx q[1];
rz(1.281339) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0616333) q[3];
sx q[3];
rz(-0.90583767) q[3];
sx q[3];
rz(3.0388415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5235644) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(0.024624126) q[2];
rz(2.5060182) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(0.88328254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4516975) q[0];
sx q[0];
rz(-0.30245936) q[0];
sx q[0];
rz(0.46689335) q[0];
rz(-3.0310071) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(1.0708403) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8798053) q[0];
sx q[0];
rz(-2.6222485) q[0];
sx q[0];
rz(-1.0762317) q[0];
x q[1];
rz(-1.3883287) q[2];
sx q[2];
rz(-1.8844205) q[2];
sx q[2];
rz(2.8387866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9164898) q[1];
sx q[1];
rz(-1.1943294) q[1];
sx q[1];
rz(2.589499) q[1];
x q[2];
rz(-0.91584622) q[3];
sx q[3];
rz(-2.1612242) q[3];
sx q[3];
rz(-2.6797323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0232627) q[2];
sx q[2];
rz(-1.7795965) q[2];
sx q[2];
rz(2.2892717) q[2];
rz(-0.037633745) q[3];
sx q[3];
rz(-1.9084385) q[3];
sx q[3];
rz(-0.5948624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0090050176) q[0];
sx q[0];
rz(-1.9629033) q[0];
sx q[0];
rz(-1.8527385) q[0];
rz(-1.6291078) q[1];
sx q[1];
rz(-1.1237203) q[1];
sx q[1];
rz(-1.1423133) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5286251) q[0];
sx q[0];
rz(-0.36699793) q[0];
sx q[0];
rz(-0.27530833) q[0];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89216512) q[1];
sx q[1];
rz(-1.2293503) q[1];
sx q[1];
rz(1.83431) q[1];
x q[2];
rz(2.8303538) q[3];
sx q[3];
rz(-1.7390307) q[3];
sx q[3];
rz(2.5995035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(-0.12040559) q[2];
rz(-1.985792) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8026546) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(1.4539723) q[0];
rz(1.5974207) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(2.6905751) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9748443) q[0];
sx q[0];
rz(-1.5112425) q[0];
sx q[0];
rz(0.3015428) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1233929) q[2];
sx q[2];
rz(-2.7180053) q[2];
sx q[2];
rz(-2.6747963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4828795) q[1];
sx q[1];
rz(-2.6900243) q[1];
sx q[1];
rz(0.2803448) q[1];
rz(-pi) q[2];
rz(-1.5805832) q[3];
sx q[3];
rz(-0.76044816) q[3];
sx q[3];
rz(1.1475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1958127) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(-1.3607402) q[2];
rz(1.0055379) q[3];
sx q[3];
rz(-1.1755627) q[3];
sx q[3];
rz(1.3895234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176395) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(-0.70607591) q[0];
rz(-1.5977244) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(2.2854038) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5498915) q[0];
sx q[0];
rz(-0.020792637) q[0];
sx q[0];
rz(-1.2077622) q[0];
rz(-pi) q[1];
rz(0.43120877) q[2];
sx q[2];
rz(-2.2715306) q[2];
sx q[2];
rz(-1.3740115) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4256712) q[1];
sx q[1];
rz(-1.7737264) q[1];
sx q[1];
rz(0.58260609) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13369067) q[3];
sx q[3];
rz(-1.0858337) q[3];
sx q[3];
rz(-3.0669989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2989444) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(2.974158) q[2];
rz(-1.5873448) q[3];
sx q[3];
rz(-0.7730248) q[3];
sx q[3];
rz(2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3779959) q[0];
sx q[0];
rz(-1.6908228) q[0];
sx q[0];
rz(-0.3717306) q[0];
rz(-1.4338214) q[1];
sx q[1];
rz(-1.6038409) q[1];
sx q[1];
rz(2.8415714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2903295) q[0];
sx q[0];
rz(-1.6870878) q[0];
sx q[0];
rz(-2.854748) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5302363) q[2];
sx q[2];
rz(-2.1475361) q[2];
sx q[2];
rz(2.4941913) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4858497) q[1];
sx q[1];
rz(-2.2169737) q[1];
sx q[1];
rz(-1.4152526) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12436538) q[3];
sx q[3];
rz(-0.87166407) q[3];
sx q[3];
rz(1.82064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6649449) q[2];
sx q[2];
rz(-0.46755329) q[2];
sx q[2];
rz(0.42759582) q[2];
rz(-1.5852196) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(-2.3868886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4204191) q[0];
sx q[0];
rz(-2.6009646) q[0];
sx q[0];
rz(-2.6721201) q[0];
rz(-0.92533127) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(-3.1184149) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7243225) q[0];
sx q[0];
rz(-2.5886154) q[0];
sx q[0];
rz(-2.824845) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.542664) q[2];
sx q[2];
rz(-2.1492276) q[2];
sx q[2];
rz(-2.3430062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1851481) q[1];
sx q[1];
rz(-1.4063764) q[1];
sx q[1];
rz(-1.3944574) q[1];
rz(-0.84081991) q[3];
sx q[3];
rz(-1.7125907) q[3];
sx q[3];
rz(-0.26103448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8699441) q[2];
sx q[2];
rz(-1.9359438) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8563817) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(2.8648227) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(-2.4527373) q[2];
sx q[2];
rz(-1.8102472) q[2];
sx q[2];
rz(2.3285338) q[2];
rz(-2.6071045) q[3];
sx q[3];
rz(-1.5880924) q[3];
sx q[3];
rz(-1.5678035) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
