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
rz(0.8241325) q[0];
sx q[0];
rz(-1.9793408) q[0];
sx q[0];
rz(0.21221575) q[0];
rz(-2.4556887) q[1];
sx q[1];
rz(-0.75690126) q[1];
sx q[1];
rz(2.8855355) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4859897) q[0];
sx q[0];
rz(-2.0168883) q[0];
sx q[0];
rz(-2.2091876) q[0];
rz(-pi) q[1];
x q[1];
rz(1.009887) q[2];
sx q[2];
rz(-1.6445064) q[2];
sx q[2];
rz(-0.91154237) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5642173) q[1];
sx q[1];
rz(-1.7322473) q[1];
sx q[1];
rz(-0.20906469) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0388548) q[3];
sx q[3];
rz(-1.3016455) q[3];
sx q[3];
rz(1.0894964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.889692) q[2];
sx q[2];
rz(-0.32521853) q[2];
sx q[2];
rz(-2.0528117) q[2];
rz(-0.47993663) q[3];
sx q[3];
rz(-1.1705385) q[3];
sx q[3];
rz(-2.0282733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8460409) q[0];
sx q[0];
rz(-0.22838455) q[0];
sx q[0];
rz(-0.26902714) q[0];
rz(-2.524952) q[1];
sx q[1];
rz(-0.65526217) q[1];
sx q[1];
rz(2.2356967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3377853) q[0];
sx q[0];
rz(-1.644994) q[0];
sx q[0];
rz(0.78671771) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3712641) q[2];
sx q[2];
rz(-0.69122756) q[2];
sx q[2];
rz(0.28797418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0695659) q[1];
sx q[1];
rz(-2.8588738) q[1];
sx q[1];
rz(0.43510423) q[1];
rz(-2.5652418) q[3];
sx q[3];
rz(-0.64223993) q[3];
sx q[3];
rz(-3.0321966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9783832) q[2];
sx q[2];
rz(-0.78028148) q[2];
sx q[2];
rz(3.1261463) q[2];
rz(0.77203006) q[3];
sx q[3];
rz(-2.6431712) q[3];
sx q[3];
rz(-1.8933403) q[3];
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
rz(pi/2) q[0];
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
rz(-0.60927272) q[0];
sx q[0];
rz(-2.5072704) q[0];
sx q[0];
rz(-0.79737216) q[0];
rz(-2.4728921) q[1];
sx q[1];
rz(-2.1233605) q[1];
sx q[1];
rz(0.55577898) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41125339) q[0];
sx q[0];
rz(-0.78124917) q[0];
sx q[0];
rz(2.8528105) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3923116) q[2];
sx q[2];
rz(-2.4694519) q[2];
sx q[2];
rz(-1.7824506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2712471) q[1];
sx q[1];
rz(-1.5588798) q[1];
sx q[1];
rz(2.4153277) q[1];
x q[2];
rz(0.45826073) q[3];
sx q[3];
rz(-0.76253676) q[3];
sx q[3];
rz(0.5173232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3681616) q[2];
sx q[2];
rz(-1.1661466) q[2];
sx q[2];
rz(-1.2097166) q[2];
rz(-3.0220253) q[3];
sx q[3];
rz(-0.60496324) q[3];
sx q[3];
rz(2.1623478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46215737) q[0];
sx q[0];
rz(-1.4428416) q[0];
sx q[0];
rz(2.6642642) q[0];
rz(0.10074549) q[1];
sx q[1];
rz(-2.0858177) q[1];
sx q[1];
rz(0.13564067) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0654146) q[0];
sx q[0];
rz(-1.9913902) q[0];
sx q[0];
rz(-0.89837822) q[0];
rz(-0.56265752) q[2];
sx q[2];
rz(-0.51501319) q[2];
sx q[2];
rz(0.80660179) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1592334) q[1];
sx q[1];
rz(-1.9894588) q[1];
sx q[1];
rz(-1.406698) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78285406) q[3];
sx q[3];
rz(-1.7373891) q[3];
sx q[3];
rz(1.596791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.2123083) q[2];
sx q[2];
rz(-1.4886651) q[2];
sx q[2];
rz(-2.6279214) q[2];
rz(2.9786927) q[3];
sx q[3];
rz(-0.25560156) q[3];
sx q[3];
rz(-2.0388849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5279919) q[0];
sx q[0];
rz(-3.1311212) q[0];
sx q[0];
rz(-0.42798671) q[0];
rz(1.9635268) q[1];
sx q[1];
rz(-1.4617498) q[1];
sx q[1];
rz(1.7286495) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1278541) q[0];
sx q[0];
rz(-1.1098521) q[0];
sx q[0];
rz(2.7147033) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5287077) q[2];
sx q[2];
rz(-1.5344991) q[2];
sx q[2];
rz(-0.55755471) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6742432) q[1];
sx q[1];
rz(-1.269264) q[1];
sx q[1];
rz(-1.5128972) q[1];
rz(2.574202) q[3];
sx q[3];
rz(-1.1674644) q[3];
sx q[3];
rz(-1.698026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.596375) q[2];
sx q[2];
rz(-1.7866106) q[2];
sx q[2];
rz(2.6522816) q[2];
rz(2.4423068) q[3];
sx q[3];
rz(-1.0438865) q[3];
sx q[3];
rz(-1.4780686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(2.3938703) q[0];
sx q[0];
rz(-1.007653) q[0];
sx q[0];
rz(1.3442159) q[0];
rz(-2.4463553) q[1];
sx q[1];
rz(-1.8724915) q[1];
sx q[1];
rz(-2.7183547) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3205954) q[0];
sx q[0];
rz(-1.508983) q[0];
sx q[0];
rz(0.16979102) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0401523) q[2];
sx q[2];
rz(-0.82221088) q[2];
sx q[2];
rz(-1.492983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9452104) q[1];
sx q[1];
rz(-1.7194413) q[1];
sx q[1];
rz(-2.3988612) q[1];
x q[2];
rz(2.845351) q[3];
sx q[3];
rz(-2.9518173) q[3];
sx q[3];
rz(1.1452183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3678579) q[2];
sx q[2];
rz(-2.9025142) q[2];
sx q[2];
rz(-1.9284922) q[2];
rz(0.086094543) q[3];
sx q[3];
rz(-1.2505069) q[3];
sx q[3];
rz(-1.7566173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4940779) q[0];
sx q[0];
rz(-2.8120742) q[0];
sx q[0];
rz(-0.20429097) q[0];
rz(2.792141) q[1];
sx q[1];
rz(-1.5705669) q[1];
sx q[1];
rz(-1.516516) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0030999936) q[0];
sx q[0];
rz(-0.47687832) q[0];
sx q[0];
rz(-1.3365227) q[0];
rz(-pi) q[1];
rz(1.1687241) q[2];
sx q[2];
rz(-1.8622145) q[2];
sx q[2];
rz(0.60044392) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2353846) q[1];
sx q[1];
rz(-1.6170119) q[1];
sx q[1];
rz(1.6711722) q[1];
rz(-2.1488993) q[3];
sx q[3];
rz(-1.5823964) q[3];
sx q[3];
rz(-1.3101526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.99882952) q[2];
sx q[2];
rz(-0.62935722) q[2];
sx q[2];
rz(-0.91267419) q[2];
rz(-3.1257889) q[3];
sx q[3];
rz(-1.088257) q[3];
sx q[3];
rz(2.6921932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34761804) q[0];
sx q[0];
rz(-2.0517218) q[0];
sx q[0];
rz(2.3225978) q[0];
rz(-0.35797572) q[1];
sx q[1];
rz(-1.1281818) q[1];
sx q[1];
rz(-2.5111759) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47641817) q[0];
sx q[0];
rz(-0.074284241) q[0];
sx q[0];
rz(-1.3528385) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6825601) q[2];
sx q[2];
rz(-0.9046281) q[2];
sx q[2];
rz(0.55202548) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9832657) q[1];
sx q[1];
rz(-0.7559146) q[1];
sx q[1];
rz(-0.96228881) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.181545) q[3];
sx q[3];
rz(-2.4195101) q[3];
sx q[3];
rz(2.275265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36998746) q[2];
sx q[2];
rz(-0.5094173) q[2];
sx q[2];
rz(0.39311692) q[2];
rz(0.76817948) q[3];
sx q[3];
rz(-0.79777515) q[3];
sx q[3];
rz(-1.2099077) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2792252) q[0];
sx q[0];
rz(-0.51280642) q[0];
sx q[0];
rz(2.7592036) q[0];
rz(-2.6928316) q[1];
sx q[1];
rz(-0.02350137) q[1];
sx q[1];
rz(-0.99865595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778393) q[0];
sx q[0];
rz(-2.4344517) q[0];
sx q[0];
rz(2.7718443) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3330179) q[2];
sx q[2];
rz(-0.3467803) q[2];
sx q[2];
rz(-1.4245167) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.71357766) q[1];
sx q[1];
rz(-0.35555118) q[1];
sx q[1];
rz(0.6014892) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5019997) q[3];
sx q[3];
rz(-2.2584887) q[3];
sx q[3];
rz(-2.6042134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0465595) q[2];
sx q[2];
rz(-2.3456489) q[2];
sx q[2];
rz(0.98540068) q[2];
rz(-2.4424148) q[3];
sx q[3];
rz(-2.2907084) q[3];
sx q[3];
rz(0.053475577) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92626524) q[0];
sx q[0];
rz(-0.43020058) q[0];
sx q[0];
rz(2.3660124) q[0];
rz(0.27014488) q[1];
sx q[1];
rz(-1.6855449) q[1];
sx q[1];
rz(2.5712579) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9035868) q[0];
sx q[0];
rz(-2.4232583) q[0];
sx q[0];
rz(0.75634457) q[0];
rz(-pi) q[1];
rz(0.10639543) q[2];
sx q[2];
rz(-2.0465851) q[2];
sx q[2];
rz(1.758213) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8034035) q[1];
sx q[1];
rz(-1.8042776) q[1];
sx q[1];
rz(1.6031053) q[1];
rz(-pi) q[2];
rz(1.4943284) q[3];
sx q[3];
rz(-0.88275331) q[3];
sx q[3];
rz(2.9043759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.21738416) q[2];
sx q[2];
rz(-1.9839828) q[2];
sx q[2];
rz(0.50979924) q[2];
rz(0.23586759) q[3];
sx q[3];
rz(-2.9678952) q[3];
sx q[3];
rz(-2.164446) q[3];
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
rz(-1.8248642) q[0];
sx q[0];
rz(-1.4126128) q[0];
sx q[0];
rz(1.8870507) q[0];
rz(-0.53312373) q[1];
sx q[1];
rz(-1.4630547) q[1];
sx q[1];
rz(-2.0933082) q[1];
rz(-1.539408) q[2];
sx q[2];
rz(-2.7161408) q[2];
sx q[2];
rz(-1.3795992) q[2];
rz(0.20797603) q[3];
sx q[3];
rz(-2.6317876) q[3];
sx q[3];
rz(-0.77710487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
