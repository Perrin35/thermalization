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
rz(-0.8166135) q[0];
sx q[0];
rz(-0.54083523) q[0];
sx q[0];
rz(1.1582561) q[0];
rz(6.3732014) q[1];
sx q[1];
rz(6.7572588) q[1];
sx q[1];
rz(13.401539) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1285889) q[0];
sx q[0];
rz(-0.79377257) q[0];
sx q[0];
rz(2.6966266) q[0];
rz(-pi) q[1];
rz(1.8878172) q[2];
sx q[2];
rz(-2.1875513) q[2];
sx q[2];
rz(-2.1932909) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0660634) q[1];
sx q[1];
rz(-1.3261686) q[1];
sx q[1];
rz(-2.7615158) q[1];
x q[2];
rz(-2.4477169) q[3];
sx q[3];
rz(-2.3531647) q[3];
sx q[3];
rz(2.3027248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0917255) q[2];
sx q[2];
rz(-0.89605248) q[2];
sx q[2];
rz(0.33207616) q[2];
rz(0.24886985) q[3];
sx q[3];
rz(-1.9460461) q[3];
sx q[3];
rz(2.2539049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2700972) q[0];
sx q[0];
rz(-0.29093727) q[0];
sx q[0];
rz(-2.6089597) q[0];
rz(3.0337785) q[1];
sx q[1];
rz(-1.1153406) q[1];
sx q[1];
rz(0.32049387) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5607213) q[0];
sx q[0];
rz(-2.221829) q[0];
sx q[0];
rz(0.080115155) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17283792) q[2];
sx q[2];
rz(-1.342475) q[2];
sx q[2];
rz(-0.31712118) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.962226) q[1];
sx q[1];
rz(-1.4186064) q[1];
sx q[1];
rz(2.9801912) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32933195) q[3];
sx q[3];
rz(-1.7747702) q[3];
sx q[3];
rz(-0.70131174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8104441) q[2];
sx q[2];
rz(-0.30511567) q[2];
sx q[2];
rz(1.6395052) q[2];
rz(0.33809996) q[3];
sx q[3];
rz(-0.88518849) q[3];
sx q[3];
rz(2.6147208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6239887) q[0];
sx q[0];
rz(-1.232134) q[0];
sx q[0];
rz(0.35743085) q[0];
rz(-0.39237157) q[1];
sx q[1];
rz(-2.3452499) q[1];
sx q[1];
rz(1.2145112) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93335184) q[0];
sx q[0];
rz(-2.9989971) q[0];
sx q[0];
rz(-1.6928133) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8896595) q[2];
sx q[2];
rz(-2.5802543) q[2];
sx q[2];
rz(2.6643945) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4068661) q[1];
sx q[1];
rz(-2.7541783) q[1];
sx q[1];
rz(-1.563036) q[1];
x q[2];
rz(-1.6271934) q[3];
sx q[3];
rz(-0.58844756) q[3];
sx q[3];
rz(1.1882888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3452722) q[2];
sx q[2];
rz(-0.97769633) q[2];
sx q[2];
rz(-0.39682445) q[2];
rz(-1.8848298) q[3];
sx q[3];
rz(-1.7233012) q[3];
sx q[3];
rz(-0.61029148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9339555) q[0];
sx q[0];
rz(-1.7455245) q[0];
sx q[0];
rz(-1.9130094) q[0];
rz(1.2127016) q[1];
sx q[1];
rz(-0.96114254) q[1];
sx q[1];
rz(-0.62087762) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1771496) q[0];
sx q[0];
rz(-1.4040134) q[0];
sx q[0];
rz(2.0509999) q[0];
rz(0.75350301) q[2];
sx q[2];
rz(-2.916159) q[2];
sx q[2];
rz(1.4788675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.32737) q[1];
sx q[1];
rz(-1.8291706) q[1];
sx q[1];
rz(1.8381005) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7589448) q[3];
sx q[3];
rz(-1.4605902) q[3];
sx q[3];
rz(-0.73303849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88317251) q[2];
sx q[2];
rz(-2.0032538) q[2];
sx q[2];
rz(2.9849198) q[2];
rz(0.91642085) q[3];
sx q[3];
rz(-2.1082924) q[3];
sx q[3];
rz(-2.5887183) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857392) q[0];
sx q[0];
rz(-1.1612949) q[0];
sx q[0];
rz(2.7440985) q[0];
rz(0.022857895) q[1];
sx q[1];
rz(-0.50067478) q[1];
sx q[1];
rz(-0.674725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2215421) q[0];
sx q[0];
rz(-1.3617317) q[0];
sx q[0];
rz(1.8110214) q[0];
x q[1];
rz(3.0176388) q[2];
sx q[2];
rz(-1.6646241) q[2];
sx q[2];
rz(1.7403062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51395638) q[1];
sx q[1];
rz(-1.7412724) q[1];
sx q[1];
rz(1.7254616) q[1];
rz(-pi) q[2];
rz(2.8679315) q[3];
sx q[3];
rz(-1.869264) q[3];
sx q[3];
rz(-0.65200114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5270093) q[2];
sx q[2];
rz(-0.60636568) q[2];
sx q[2];
rz(2.442339) q[2];
rz(-1.7953385) q[3];
sx q[3];
rz(-2.5739539) q[3];
sx q[3];
rz(-0.7114555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72494495) q[0];
sx q[0];
rz(-0.41794932) q[0];
sx q[0];
rz(0.39837343) q[0];
rz(2.1669972) q[1];
sx q[1];
rz(-1.83788) q[1];
sx q[1];
rz(1.7030565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8038626) q[0];
sx q[0];
rz(-0.32050214) q[0];
sx q[0];
rz(2.8021373) q[0];
rz(-pi) q[1];
rz(-1.826615) q[2];
sx q[2];
rz(-1.5607569) q[2];
sx q[2];
rz(-0.61983392) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12425646) q[1];
sx q[1];
rz(-1.3385695) q[1];
sx q[1];
rz(0.88084765) q[1];
x q[2];
rz(-2.6468032) q[3];
sx q[3];
rz(-1.7980709) q[3];
sx q[3];
rz(-2.761043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67821104) q[2];
sx q[2];
rz(-1.3621829) q[2];
sx q[2];
rz(-1.8360651) q[2];
rz(1.8259004) q[3];
sx q[3];
rz(-1.3633599) q[3];
sx q[3];
rz(-0.8684043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797794) q[0];
sx q[0];
rz(-2.7257305) q[0];
sx q[0];
rz(-1.4965936) q[0];
rz(0.98622259) q[1];
sx q[1];
rz(-1.5602427) q[1];
sx q[1];
rz(2.644002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3517036) q[0];
sx q[0];
rz(-1.6663807) q[0];
sx q[0];
rz(1.8259177) q[0];
rz(-0.53908252) q[2];
sx q[2];
rz(-2.1240799) q[2];
sx q[2];
rz(2.6607571) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.033173) q[1];
sx q[1];
rz(-0.42627508) q[1];
sx q[1];
rz(1.1854965) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6646181) q[3];
sx q[3];
rz(-1.2872496) q[3];
sx q[3];
rz(-0.16167262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8387973) q[2];
sx q[2];
rz(-1.5747728) q[2];
sx q[2];
rz(2.8509169) q[2];
rz(0.21026462) q[3];
sx q[3];
rz(-0.99431521) q[3];
sx q[3];
rz(-1.0342342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57598376) q[0];
sx q[0];
rz(-1.9714332) q[0];
sx q[0];
rz(1.1822816) q[0];
rz(2.9871509) q[1];
sx q[1];
rz(-1.4192105) q[1];
sx q[1];
rz(-0.90528893) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2574999) q[0];
sx q[0];
rz(-1.5587036) q[0];
sx q[0];
rz(1.2980677) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6851875) q[2];
sx q[2];
rz(-1.3564566) q[2];
sx q[2];
rz(-2.6021007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.33069995) q[1];
sx q[1];
rz(-2.3726683) q[1];
sx q[1];
rz(-2.9309209) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2298814) q[3];
sx q[3];
rz(-1.3671759) q[3];
sx q[3];
rz(0.019817185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0901383) q[2];
sx q[2];
rz(-1.5608414) q[2];
sx q[2];
rz(2.1913989) q[2];
rz(0.072362445) q[3];
sx q[3];
rz(-1.7642998) q[3];
sx q[3];
rz(0.70934057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7224834) q[0];
sx q[0];
rz(-2.1868732) q[0];
sx q[0];
rz(-1.0954274) q[0];
rz(2.0504045) q[1];
sx q[1];
rz(-0.90091101) q[1];
sx q[1];
rz(-2.1464164) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4596915) q[0];
sx q[0];
rz(-0.90666214) q[0];
sx q[0];
rz(-1.4891529) q[0];
rz(-pi) q[1];
rz(1.9068933) q[2];
sx q[2];
rz(-1.8310412) q[2];
sx q[2];
rz(-2.0356095) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.02444754) q[1];
sx q[1];
rz(-3.0452485) q[1];
sx q[1];
rz(2.6481555) q[1];
rz(-2.9561958) q[3];
sx q[3];
rz(-1.7668006) q[3];
sx q[3];
rz(-2.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5114078) q[2];
sx q[2];
rz(-0.46773043) q[2];
sx q[2];
rz(-0.67997813) q[2];
rz(-1.6857111) q[3];
sx q[3];
rz(-0.89434353) q[3];
sx q[3];
rz(1.1792012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42530123) q[0];
sx q[0];
rz(-2.20708) q[0];
sx q[0];
rz(2.5262078) q[0];
rz(1.3353434) q[1];
sx q[1];
rz(-1.5348624) q[1];
sx q[1];
rz(1.7937484) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16432654) q[0];
sx q[0];
rz(-1.5249671) q[0];
sx q[0];
rz(0.35542506) q[0];
x q[1];
rz(0.87839076) q[2];
sx q[2];
rz(-0.88365245) q[2];
sx q[2];
rz(2.2855482) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3879229) q[1];
sx q[1];
rz(-1.7311454) q[1];
sx q[1];
rz(2.760375) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75611434) q[3];
sx q[3];
rz(-1.4817837) q[3];
sx q[3];
rz(1.4579888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24796692) q[2];
sx q[2];
rz(-1.8411571) q[2];
sx q[2];
rz(2.628053) q[2];
rz(-0.14704554) q[3];
sx q[3];
rz(-0.5439609) q[3];
sx q[3];
rz(1.8769544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87179398) q[0];
sx q[0];
rz(-0.88903058) q[0];
sx q[0];
rz(-0.24118184) q[0];
rz(1.3006032) q[1];
sx q[1];
rz(-2.0722957) q[1];
sx q[1];
rz(3.121079) q[1];
rz(-1.2215963) q[2];
sx q[2];
rz(-0.39633718) q[2];
sx q[2];
rz(0.63971165) q[2];
rz(2.5138598) q[3];
sx q[3];
rz(-1.9398324) q[3];
sx q[3];
rz(2.0331665) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
