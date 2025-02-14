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
rz(-0.095280401) q[1];
sx q[1];
rz(2.4089101) q[1];
sx q[1];
rz(10.664193) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026442083) q[0];
sx q[0];
rz(-1.7998994) q[0];
sx q[0];
rz(1.312448) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26387604) q[2];
sx q[2];
rz(-1.7603163) q[2];
sx q[2];
rz(-0.60631982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3500786) q[1];
sx q[1];
rz(-2.8100612) q[1];
sx q[1];
rz(-0.091011957) q[1];
rz(1.0765355) q[3];
sx q[3];
rz(-2.0856033) q[3];
sx q[3];
rz(0.47806397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6713082) q[2];
sx q[2];
rz(-1.6613864) q[2];
sx q[2];
rz(-1.3489464) q[2];
rz(2.3979483) q[3];
sx q[3];
rz(-2.8642004) q[3];
sx q[3];
rz(-0.17717895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0002366) q[0];
sx q[0];
rz(-1.5445222) q[0];
sx q[0];
rz(0.29362383) q[0];
rz(1.0307505) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(-0.34293276) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6688419) q[0];
sx q[0];
rz(-1.051051) q[0];
sx q[0];
rz(0.23730554) q[0];
rz(-pi) q[1];
rz(1.198771) q[2];
sx q[2];
rz(-1.4497533) q[2];
sx q[2];
rz(2.4189359) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.575042) q[1];
sx q[1];
rz(-2.2600265) q[1];
sx q[1];
rz(-1.8495967) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2323805) q[3];
sx q[3];
rz(-0.65753257) q[3];
sx q[3];
rz(2.7595208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0127516) q[2];
sx q[2];
rz(-1.7495456) q[2];
sx q[2];
rz(-0.7592321) q[2];
rz(-3.1208842) q[3];
sx q[3];
rz(-1.8755251) q[3];
sx q[3];
rz(0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6193806) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(3.1371327) q[0];
rz(-2.4941173) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(2.3562145) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1566832) q[0];
sx q[0];
rz(-1.4117318) q[0];
sx q[0];
rz(1.6522264) q[0];
x q[1];
rz(0.17526971) q[2];
sx q[2];
rz(-1.2614377) q[2];
sx q[2];
rz(0.64495211) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
x q[2];
rz(2.4786199) q[3];
sx q[3];
rz(-2.1372652) q[3];
sx q[3];
rz(-0.86323767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.37041) q[2];
sx q[2];
rz(-0.9202756) q[2];
sx q[2];
rz(1.7837589) q[2];
rz(-0.34058288) q[3];
sx q[3];
rz(-2.1850696) q[3];
sx q[3];
rz(0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99712813) q[0];
sx q[0];
rz(-1.0839394) q[0];
sx q[0];
rz(2.2564364) q[0];
rz(0.74686933) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(-1.0964099) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.261026) q[0];
sx q[0];
rz(-1.9184904) q[0];
sx q[0];
rz(0.33190042) q[0];
x q[1];
rz(-2.0463767) q[2];
sx q[2];
rz(-0.066844373) q[2];
sx q[2];
rz(-2.9796556) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8453522) q[1];
sx q[1];
rz(-1.8295604) q[1];
sx q[1];
rz(2.6654763) q[1];
rz(-pi) q[2];
rz(2.5853852) q[3];
sx q[3];
rz(-2.3283151) q[3];
sx q[3];
rz(-2.5089056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5235644) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(0.024624126) q[2];
rz(-0.63557449) q[3];
sx q[3];
rz(-0.98819757) q[3];
sx q[3];
rz(2.2583101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6898952) q[0];
sx q[0];
rz(-0.30245936) q[0];
sx q[0];
rz(-0.46689335) q[0];
rz(3.0310071) q[1];
sx q[1];
rz(-0.46375912) q[1];
sx q[1];
rz(1.0708403) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2617874) q[0];
sx q[0];
rz(-0.51934411) q[0];
sx q[0];
rz(-1.0762317) q[0];
x q[1];
rz(2.8230225) q[2];
sx q[2];
rz(-1.397322) q[2];
sx q[2];
rz(-1.3248548) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5683978) q[1];
sx q[1];
rz(-1.0612951) q[1];
sx q[1];
rz(-2.0054818) q[1];
x q[2];
rz(0.91584622) q[3];
sx q[3];
rz(-0.98036843) q[3];
sx q[3];
rz(-2.6797323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.11833) q[2];
sx q[2];
rz(-1.3619962) q[2];
sx q[2];
rz(-0.85232097) q[2];
rz(3.1039589) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(-2.5467303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1325876) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(-1.8527385) q[0];
rz(-1.6291078) q[1];
sx q[1];
rz(-1.1237203) q[1];
sx q[1];
rz(-1.1423133) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4415671) q[0];
sx q[0];
rz(-1.4730994) q[0];
sx q[0];
rz(-0.35432451) q[0];
rz(-pi) q[1];
rz(-3.0794607) q[2];
sx q[2];
rz(-1.3634063) q[2];
sx q[2];
rz(-1.0203938) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.89216512) q[1];
sx q[1];
rz(-1.9122423) q[1];
sx q[1];
rz(-1.3072827) q[1];
x q[2];
rz(-1.7473502) q[3];
sx q[3];
rz(-1.8774967) q[3];
sx q[3];
rz(-0.97489417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8196572) q[2];
sx q[2];
rz(-1.8737996) q[2];
sx q[2];
rz(-3.0211871) q[2];
rz(-1.1558007) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(-2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33893809) q[0];
sx q[0];
rz(-1.3405565) q[0];
sx q[0];
rz(1.4539723) q[0];
rz(1.5441719) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(-2.6905751) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9748443) q[0];
sx q[0];
rz(-1.5112425) q[0];
sx q[0];
rz(2.8400499) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1233929) q[2];
sx q[2];
rz(-0.42358735) q[2];
sx q[2];
rz(0.46679631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6587131) q[1];
sx q[1];
rz(-2.6900243) q[1];
sx q[1];
rz(-0.2803448) q[1];
rz(-pi) q[2];
rz(-0.0093098442) q[3];
sx q[3];
rz(-2.3311989) q[3];
sx q[3];
rz(-1.9805465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94577998) q[2];
sx q[2];
rz(-1.7273629) q[2];
sx q[2];
rz(-1.7808524) q[2];
rz(1.0055379) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(-1.3895234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176395) q[0];
sx q[0];
rz(-0.12549505) q[0];
sx q[0];
rz(-2.4355167) q[0];
rz(1.5977244) q[1];
sx q[1];
rz(-1.4223301) q[1];
sx q[1];
rz(-0.85618883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2285952) q[0];
sx q[0];
rz(-1.5513591) q[0];
sx q[0];
rz(-3.134208) q[0];
x q[1];
rz(2.7103839) q[2];
sx q[2];
rz(-2.2715306) q[2];
sx q[2];
rz(1.3740115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.44196196) q[1];
sx q[1];
rz(-2.5285427) q[1];
sx q[1];
rz(0.35787257) q[1];
rz(-pi) q[2];
rz(-1.3230349) q[3];
sx q[3];
rz(-2.6399586) q[3];
sx q[3];
rz(-0.35546965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2989444) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(-2.974158) q[2];
rz(1.5873448) q[3];
sx q[3];
rz(-0.7730248) q[3];
sx q[3];
rz(1.0003482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7635968) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(-0.3717306) q[0];
rz(-1.7077712) q[1];
sx q[1];
rz(-1.6038409) q[1];
sx q[1];
rz(-2.8415714) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8269129) q[0];
sx q[0];
rz(-1.2859435) q[0];
sx q[0];
rz(1.4495984) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5644767) q[2];
sx q[2];
rz(-1.5367998) q[2];
sx q[2];
rz(-2.196072) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4858497) q[1];
sx q[1];
rz(-0.92461899) q[1];
sx q[1];
rz(1.7263401) q[1];
rz(1.4243218) q[3];
sx q[3];
rz(-0.70826847) q[3];
sx q[3];
rz(-2.012501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47664777) q[2];
sx q[2];
rz(-2.6740394) q[2];
sx q[2];
rz(-0.42759582) q[2];
rz(1.5852196) q[3];
sx q[3];
rz(-1.1170324) q[3];
sx q[3];
rz(-2.3868886) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4204191) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(0.46947259) q[0];
rz(-0.92533127) q[1];
sx q[1];
rz(-1.0877437) q[1];
sx q[1];
rz(-0.023177711) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7243225) q[0];
sx q[0];
rz(-0.5529772) q[0];
sx q[0];
rz(0.31674762) q[0];
rz(-pi) q[1];
rz(0.9018578) q[2];
sx q[2];
rz(-2.0624071) q[2];
sx q[2];
rz(2.7265446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7851023) q[1];
sx q[1];
rz(-1.7447326) q[1];
sx q[1];
rz(-2.9746303) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9523207) q[3];
sx q[3];
rz(-0.84976518) q[3];
sx q[3];
rz(-1.4355575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2716486) q[2];
sx q[2];
rz(-1.9359438) q[2];
sx q[2];
rz(-0.49087697) q[2];
rz(-2.0513746) q[3];
sx q[3];
rz(-2.1722983) q[3];
sx q[3];
rz(2.8276665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8563817) q[0];
sx q[0];
rz(-2.2638392) q[0];
sx q[0];
rz(1.0765156) q[0];
rz(-0.27676997) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(0.3666975) q[2];
sx q[2];
rz(-2.4187805) q[2];
sx q[2];
rz(0.4772966) q[2];
rz(1.5506977) q[3];
sx q[3];
rz(-2.1051959) q[3];
sx q[3];
rz(0.013230562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
