OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1462829) q[0];
sx q[0];
rz(0.028385552) q[0];
sx q[0];
rz(10.382494) q[0];
rz(-2.7081642) q[1];
sx q[1];
rz(-1.6038916) q[1];
sx q[1];
rz(0.62827194) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7564297) q[0];
sx q[0];
rz(-1.9397501) q[0];
sx q[0];
rz(1.4403309) q[0];
rz(-pi) q[1];
rz(-0.84282173) q[2];
sx q[2];
rz(-2.5113238) q[2];
sx q[2];
rz(0.81319076) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0121668) q[1];
sx q[1];
rz(-0.91246997) q[1];
sx q[1];
rz(2.7321187) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7946934) q[3];
sx q[3];
rz(-2.7198987) q[3];
sx q[3];
rz(-1.2661915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6750703) q[2];
sx q[2];
rz(-1.0202946) q[2];
sx q[2];
rz(2.7543219) q[2];
rz(1.9668503) q[3];
sx q[3];
rz(-1.6956804) q[3];
sx q[3];
rz(-2.2430879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3627477) q[0];
sx q[0];
rz(-1.7411106) q[0];
sx q[0];
rz(-1.8722906) q[0];
rz(-1.6954039) q[1];
sx q[1];
rz(-1.8480999) q[1];
sx q[1];
rz(-2.2206025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22494754) q[0];
sx q[0];
rz(-0.9592255) q[0];
sx q[0];
rz(-2.7782604) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3878421) q[2];
sx q[2];
rz(-0.91380807) q[2];
sx q[2];
rz(-2.8859101) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9780636) q[1];
sx q[1];
rz(-2.7394751) q[1];
sx q[1];
rz(2.0744214) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5778353) q[3];
sx q[3];
rz(-1.3448997) q[3];
sx q[3];
rz(-1.5510538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9903119) q[2];
sx q[2];
rz(-1.1996148) q[2];
sx q[2];
rz(-2.3178103) q[2];
rz(1.524205) q[3];
sx q[3];
rz(-1.5533605) q[3];
sx q[3];
rz(1.2113781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9059432) q[0];
sx q[0];
rz(-2.2610569) q[0];
sx q[0];
rz(-2.549951) q[0];
rz(-2.1982684) q[1];
sx q[1];
rz(-1.0572409) q[1];
sx q[1];
rz(-1.0234157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9647261) q[0];
sx q[0];
rz(-2.9826678) q[0];
sx q[0];
rz(-2.311241) q[0];
rz(-2.1322971) q[2];
sx q[2];
rz(-0.63234419) q[2];
sx q[2];
rz(0.48582669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.652214) q[1];
sx q[1];
rz(-1.4783597) q[1];
sx q[1];
rz(-3.1104065) q[1];
rz(-pi) q[2];
rz(-2.4571506) q[3];
sx q[3];
rz(-1.8307643) q[3];
sx q[3];
rz(2.3459159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.25980514) q[2];
sx q[2];
rz(-1.2299579) q[2];
sx q[2];
rz(1.7477431) q[2];
rz(-1.7143837) q[3];
sx q[3];
rz(-2.2673159) q[3];
sx q[3];
rz(-2.8425596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3706936) q[0];
sx q[0];
rz(-2.735266) q[0];
sx q[0];
rz(-2.4942177) q[0];
rz(-0.24230832) q[1];
sx q[1];
rz(-1.7055885) q[1];
sx q[1];
rz(1.6829596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1636617) q[0];
sx q[0];
rz(-0.35172281) q[0];
sx q[0];
rz(0.39885421) q[0];
rz(-0.91084852) q[2];
sx q[2];
rz(-1.0648515) q[2];
sx q[2];
rz(0.3163704) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7902271) q[1];
sx q[1];
rz(-2.0560762) q[1];
sx q[1];
rz(-0.041260551) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69497739) q[3];
sx q[3];
rz(-0.91013346) q[3];
sx q[3];
rz(-1.7373178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4316537) q[2];
sx q[2];
rz(-0.84448758) q[2];
sx q[2];
rz(-1.9169774) q[2];
rz(-2.8905408) q[3];
sx q[3];
rz(-2.4243088) q[3];
sx q[3];
rz(0.93799463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(-1.5994023) q[0];
sx q[0];
rz(-1.1325862) q[0];
sx q[0];
rz(0.19790025) q[0];
rz(-1.9056162) q[1];
sx q[1];
rz(-2.7695152) q[1];
sx q[1];
rz(3.1370251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76659996) q[0];
sx q[0];
rz(-0.99762929) q[0];
sx q[0];
rz(-1.8660924) q[0];
rz(-1.4312237) q[2];
sx q[2];
rz(-2.6032631) q[2];
sx q[2];
rz(-1.1280504) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.63268836) q[1];
sx q[1];
rz(-1.4211854) q[1];
sx q[1];
rz(-2.8652111) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3602348) q[3];
sx q[3];
rz(-1.9145085) q[3];
sx q[3];
rz(2.3794019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4602451) q[2];
sx q[2];
rz(-0.494151) q[2];
sx q[2];
rz(-2.8013308) q[2];
rz(1.6081238) q[3];
sx q[3];
rz(-1.7483277) q[3];
sx q[3];
rz(-1.5197915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22572868) q[0];
sx q[0];
rz(-0.69750834) q[0];
sx q[0];
rz(1.9541784) q[0];
rz(1.4215218) q[1];
sx q[1];
rz(-0.41821304) q[1];
sx q[1];
rz(3.0112867) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7447259) q[0];
sx q[0];
rz(-1.5123741) q[0];
sx q[0];
rz(1.6944042) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88779403) q[2];
sx q[2];
rz(-2.0468132) q[2];
sx q[2];
rz(-2.7107875) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7769717) q[1];
sx q[1];
rz(-1.7386308) q[1];
sx q[1];
rz(1.5701152) q[1];
rz(-pi) q[2];
rz(0.0099139307) q[3];
sx q[3];
rz(-1.3080721) q[3];
sx q[3];
rz(-0.81068057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80798951) q[2];
sx q[2];
rz(-1.5895867) q[2];
sx q[2];
rz(-2.1017334) q[2];
rz(-0.17397675) q[3];
sx q[3];
rz(-2.0128553) q[3];
sx q[3];
rz(-0.66966331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.533605) q[0];
sx q[0];
rz(-1.0796115) q[0];
sx q[0];
rz(0.84130353) q[0];
rz(-1.3249506) q[1];
sx q[1];
rz(-0.94580301) q[1];
sx q[1];
rz(-1.673505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8155839) q[0];
sx q[0];
rz(-1.2213314) q[0];
sx q[0];
rz(-0.38527003) q[0];
rz(-3.1300342) q[2];
sx q[2];
rz(-2.0894139) q[2];
sx q[2];
rz(2.8457038) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0206179) q[1];
sx q[1];
rz(-1.0134032) q[1];
sx q[1];
rz(-1.0110028) q[1];
rz(-pi) q[2];
rz(0.87196799) q[3];
sx q[3];
rz(-0.45365327) q[3];
sx q[3];
rz(-3.0986971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86772052) q[2];
sx q[2];
rz(-2.2404859) q[2];
sx q[2];
rz(-0.64485288) q[2];
rz(1.0817179) q[3];
sx q[3];
rz(-2.7223301) q[3];
sx q[3];
rz(-2.4206415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69407392) q[0];
sx q[0];
rz(-0.17386359) q[0];
sx q[0];
rz(-2.7287667) q[0];
rz(1.5326356) q[1];
sx q[1];
rz(-2.2282579) q[1];
sx q[1];
rz(-2.434381) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9291973) q[0];
sx q[0];
rz(-1.2629422) q[0];
sx q[0];
rz(2.3677127) q[0];
rz(-pi) q[1];
rz(-2.3995072) q[2];
sx q[2];
rz(-1.4550503) q[2];
sx q[2];
rz(-1.5884233) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78145786) q[1];
sx q[1];
rz(-0.64720045) q[1];
sx q[1];
rz(-0.57768808) q[1];
rz(-1.8052638) q[3];
sx q[3];
rz(-0.55113652) q[3];
sx q[3];
rz(-0.41195991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0088719) q[2];
sx q[2];
rz(-2.8421695) q[2];
sx q[2];
rz(-0.51187619) q[2];
rz(2.4608608) q[3];
sx q[3];
rz(-1.778089) q[3];
sx q[3];
rz(0.16630047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7452288) q[0];
sx q[0];
rz(-1.5487211) q[0];
sx q[0];
rz(-0.45528278) q[0];
rz(-1.0633172) q[1];
sx q[1];
rz(-0.89439193) q[1];
sx q[1];
rz(3.0696226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088475835) q[0];
sx q[0];
rz(-1.6040816) q[0];
sx q[0];
rz(1.5584591) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44092245) q[2];
sx q[2];
rz(-1.5917695) q[2];
sx q[2];
rz(1.4628589) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8645975) q[1];
sx q[1];
rz(-2.2114843) q[1];
sx q[1];
rz(-2.6990141) q[1];
x q[2];
rz(1.9017327) q[3];
sx q[3];
rz(-1.8171105) q[3];
sx q[3];
rz(-2.7992751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3528184) q[2];
sx q[2];
rz(-0.41663751) q[2];
sx q[2];
rz(-0.62208661) q[2];
rz(-2.3173053) q[3];
sx q[3];
rz(-1.0930748) q[3];
sx q[3];
rz(0.85339671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033009919) q[0];
sx q[0];
rz(-1.1065296) q[0];
sx q[0];
rz(-0.95091096) q[0];
rz(1.4190326) q[1];
sx q[1];
rz(-0.483069) q[1];
sx q[1];
rz(0.61666617) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86980235) q[0];
sx q[0];
rz(-1.5464725) q[0];
sx q[0];
rz(-0.64016827) q[0];
rz(-pi) q[1];
rz(0.60994591) q[2];
sx q[2];
rz(-1.8862572) q[2];
sx q[2];
rz(-1.1300398) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92041278) q[1];
sx q[1];
rz(-1.3592459) q[1];
sx q[1];
rz(-1.5736363) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70722001) q[3];
sx q[3];
rz(-2.0248746) q[3];
sx q[3];
rz(-1.2421158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0509433) q[2];
sx q[2];
rz(-1.9885149) q[2];
sx q[2];
rz(-2.0743745) q[2];
rz(-1.0956592) q[3];
sx q[3];
rz(-1.2376384) q[3];
sx q[3];
rz(2.1777976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6296366) q[0];
sx q[0];
rz(-2.1558599) q[0];
sx q[0];
rz(2.070367) q[0];
rz(1.6596972) q[1];
sx q[1];
rz(-2.0493458) q[1];
sx q[1];
rz(0.58072166) q[1];
rz(-2.8790375) q[2];
sx q[2];
rz(-1.5947744) q[2];
sx q[2];
rz(-2.6626432) q[2];
rz(0.31286851) q[3];
sx q[3];
rz(-1.4405559) q[3];
sx q[3];
rz(1.4961989) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
