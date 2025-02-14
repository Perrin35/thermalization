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
rz(-0.38464889) q[0];
sx q[0];
rz(0.83847133) q[0];
sx q[0];
rz(5.6738927) q[0];
rz(0.72120178) q[1];
sx q[1];
rz(-0.92858044) q[1];
sx q[1];
rz(-2.0452926) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23966889) q[0];
sx q[0];
rz(-0.87920183) q[0];
sx q[0];
rz(-0.14622971) q[0];
rz(-0.9303665) q[2];
sx q[2];
rz(-1.1784679) q[2];
sx q[2];
rz(1.0898255) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.23725736) q[1];
sx q[1];
rz(-1.0544027) q[1];
sx q[1];
rz(-1.8947435) q[1];
rz(-pi) q[2];
rz(1.4601379) q[3];
sx q[3];
rz(-1.4165464) q[3];
sx q[3];
rz(1.339104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91813749) q[2];
sx q[2];
rz(-1.3104985) q[2];
sx q[2];
rz(-1.1892595) q[2];
rz(0.39092815) q[3];
sx q[3];
rz(-1.9783741) q[3];
sx q[3];
rz(-2.0093567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696058) q[0];
sx q[0];
rz(-1.7867418) q[0];
sx q[0];
rz(-1.6433486) q[0];
rz(2.4636726) q[1];
sx q[1];
rz(-1.8410212) q[1];
sx q[1];
rz(-2.4300785) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2133117) q[0];
sx q[0];
rz(-0.60325256) q[0];
sx q[0];
rz(-2.2343982) q[0];
x q[1];
rz(1.8078126) q[2];
sx q[2];
rz(-2.0227602) q[2];
sx q[2];
rz(-2.1353561) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7603858) q[1];
sx q[1];
rz(-0.35522705) q[1];
sx q[1];
rz(0.23204259) q[1];
x q[2];
rz(-1.5270803) q[3];
sx q[3];
rz(-2.5175142) q[3];
sx q[3];
rz(1.3169552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70637643) q[2];
sx q[2];
rz(-1.7855568) q[2];
sx q[2];
rz(1.0630652) q[2];
rz(-1.6478018) q[3];
sx q[3];
rz(-2.6726275) q[3];
sx q[3];
rz(-0.74023214) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4110334) q[0];
sx q[0];
rz(-0.39091245) q[0];
sx q[0];
rz(-2.539769) q[0];
rz(-2.9959294) q[1];
sx q[1];
rz(-1.0258976) q[1];
sx q[1];
rz(2.513733) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8635628) q[0];
sx q[0];
rz(-3.0042227) q[0];
sx q[0];
rz(-0.53047116) q[0];
rz(-pi) q[1];
x q[1];
rz(1.026903) q[2];
sx q[2];
rz(-0.62335194) q[2];
sx q[2];
rz(-2.198213) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.16937373) q[1];
sx q[1];
rz(-1.7986606) q[1];
sx q[1];
rz(-1.9337173) q[1];
x q[2];
rz(-0.95730036) q[3];
sx q[3];
rz(-1.9578551) q[3];
sx q[3];
rz(-1.8107896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7566028) q[2];
sx q[2];
rz(-2.0774697) q[2];
sx q[2];
rz(-2.9260054) q[2];
rz(2.0032517) q[3];
sx q[3];
rz(-1.3201069) q[3];
sx q[3];
rz(0.57083541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0067921) q[0];
sx q[0];
rz(-1.4147867) q[0];
sx q[0];
rz(3.0227645) q[0];
rz(-1.4745332) q[1];
sx q[1];
rz(-2.2524565) q[1];
sx q[1];
rz(0.80783358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0197163) q[0];
sx q[0];
rz(-0.84152475) q[0];
sx q[0];
rz(-2.3725933) q[0];
rz(0.045305552) q[2];
sx q[2];
rz(-1.8938365) q[2];
sx q[2];
rz(1.3458061) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3701757) q[1];
sx q[1];
rz(-0.96730212) q[1];
sx q[1];
rz(-2.9379528) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.019185493) q[3];
sx q[3];
rz(-2.4632255) q[3];
sx q[3];
rz(2.8205591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64041758) q[2];
sx q[2];
rz(-1.7265604) q[2];
sx q[2];
rz(-0.51587063) q[2];
rz(0.094203146) q[3];
sx q[3];
rz(-0.7754063) q[3];
sx q[3];
rz(1.5883821) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3696988) q[0];
sx q[0];
rz(-1.1006681) q[0];
sx q[0];
rz(1.9544253) q[0];
rz(2.5502603) q[1];
sx q[1];
rz(-1.2966803) q[1];
sx q[1];
rz(0.82685131) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7477404) q[0];
sx q[0];
rz(-1.6918139) q[0];
sx q[0];
rz(0.38328247) q[0];
x q[1];
rz(-0.081080699) q[2];
sx q[2];
rz(-1.0141918) q[2];
sx q[2];
rz(-0.48559819) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8137774) q[1];
sx q[1];
rz(-0.89164387) q[1];
sx q[1];
rz(-0.16736729) q[1];
rz(-pi) q[2];
rz(-0.86381819) q[3];
sx q[3];
rz(-0.70138273) q[3];
sx q[3];
rz(0.4566628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7190711) q[2];
sx q[2];
rz(-0.25455385) q[2];
sx q[2];
rz(-2.1264326) q[2];
rz(-0.90967956) q[3];
sx q[3];
rz(-1.2404975) q[3];
sx q[3];
rz(-0.66143405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5231617) q[0];
sx q[0];
rz(-1.2105415) q[0];
sx q[0];
rz(-0.30995187) q[0];
rz(3.1228206) q[1];
sx q[1];
rz(-1.680178) q[1];
sx q[1];
rz(3.054256) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7531484) q[0];
sx q[0];
rz(-2.5678621) q[0];
sx q[0];
rz(0.059772003) q[0];
rz(-pi) q[1];
rz(-0.64962663) q[2];
sx q[2];
rz(-2.1611593) q[2];
sx q[2];
rz(-0.43250674) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.046459196) q[1];
sx q[1];
rz(-1.9659871) q[1];
sx q[1];
rz(-2.1161152) q[1];
x q[2];
rz(1.3052854) q[3];
sx q[3];
rz(-0.54406057) q[3];
sx q[3];
rz(1.8484704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35149082) q[2];
sx q[2];
rz(-2.6678706) q[2];
sx q[2];
rz(2.6215485) q[2];
rz(3.0067048) q[3];
sx q[3];
rz(-1.8205732) q[3];
sx q[3];
rz(-1.7322056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670369) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(-0.38597646) q[0];
rz(2.0416253) q[1];
sx q[1];
rz(-2.4620582) q[1];
sx q[1];
rz(2.3540672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.814491) q[0];
sx q[0];
rz(-1.2145894) q[0];
sx q[0];
rz(-0.30762958) q[0];
rz(-0.50561302) q[2];
sx q[2];
rz(-0.68196251) q[2];
sx q[2];
rz(1.4346892) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0213607) q[1];
sx q[1];
rz(-1.225705) q[1];
sx q[1];
rz(-0.80806513) q[1];
rz(-1.3357032) q[3];
sx q[3];
rz(-1.4054148) q[3];
sx q[3];
rz(-0.38136417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1197352) q[2];
sx q[2];
rz(-1.3081552) q[2];
sx q[2];
rz(-1.5137399) q[2];
rz(1.9205903) q[3];
sx q[3];
rz(-1.9767714) q[3];
sx q[3];
rz(-0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30962238) q[0];
sx q[0];
rz(-1.7124875) q[0];
sx q[0];
rz(0.23304644) q[0];
rz(-1.6538992) q[1];
sx q[1];
rz(-1.033604) q[1];
sx q[1];
rz(2.1344562) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24503532) q[0];
sx q[0];
rz(-0.63525891) q[0];
sx q[0];
rz(0.59478514) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2738259) q[2];
sx q[2];
rz(-2.179702) q[2];
sx q[2];
rz(2.0021653) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1888652) q[1];
sx q[1];
rz(-0.42546526) q[1];
sx q[1];
rz(-1.2098321) q[1];
x q[2];
rz(3.0606224) q[3];
sx q[3];
rz(-0.24484466) q[3];
sx q[3];
rz(-0.096658215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2795589) q[2];
sx q[2];
rz(-1.1125914) q[2];
sx q[2];
rz(-2.2588008) q[2];
rz(2.8629996) q[3];
sx q[3];
rz(-1.9735347) q[3];
sx q[3];
rz(2.6826503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9423264) q[0];
sx q[0];
rz(-2.6658604) q[0];
sx q[0];
rz(1.5717773) q[0];
rz(2.3528631) q[1];
sx q[1];
rz(-2.1756344) q[1];
sx q[1];
rz(-2.4615361) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2325033) q[0];
sx q[0];
rz(-1.5111369) q[0];
sx q[0];
rz(-2.3807314) q[0];
rz(-pi) q[1];
rz(2.3236507) q[2];
sx q[2];
rz(-2.9688058) q[2];
sx q[2];
rz(0.30889749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9763261) q[1];
sx q[1];
rz(-1.9496456) q[1];
sx q[1];
rz(-2.5514609) q[1];
rz(0.96931501) q[3];
sx q[3];
rz(-2.318104) q[3];
sx q[3];
rz(-1.6723417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1532229) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(1.1663158) q[2];
rz(1.204528) q[3];
sx q[3];
rz(-1.322999) q[3];
sx q[3];
rz(2.5148463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.6176497) q[0];
sx q[0];
rz(-0.88215041) q[0];
sx q[0];
rz(-2.7045265) q[0];
rz(2.3433459) q[1];
sx q[1];
rz(-2.6415446) q[1];
sx q[1];
rz(2.9214568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1967314) q[0];
sx q[0];
rz(-1.5240977) q[0];
sx q[0];
rz(0.80094211) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.885845) q[2];
sx q[2];
rz(-1.3761259) q[2];
sx q[2];
rz(-2.6597629) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5497351) q[1];
sx q[1];
rz(-2.3384574) q[1];
sx q[1];
rz(1.1863043) q[1];
x q[2];
rz(1.5256568) q[3];
sx q[3];
rz(-0.76971005) q[3];
sx q[3];
rz(-1.5906281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87307125) q[2];
sx q[2];
rz(-1.2348509) q[2];
sx q[2];
rz(-0.29042563) q[2];
rz(-0.63646603) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(-1.9071473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81454043) q[0];
sx q[0];
rz(-2.4335813) q[0];
sx q[0];
rz(-2.0306564) q[0];
rz(-0.064432714) q[1];
sx q[1];
rz(-1.671052) q[1];
sx q[1];
rz(-0.64238092) q[1];
rz(2.1656373) q[2];
sx q[2];
rz(-1.5968762) q[2];
sx q[2];
rz(2.1322875) q[2];
rz(-2.7148132) q[3];
sx q[3];
rz(-0.58146324) q[3];
sx q[3];
rz(2.0406722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
