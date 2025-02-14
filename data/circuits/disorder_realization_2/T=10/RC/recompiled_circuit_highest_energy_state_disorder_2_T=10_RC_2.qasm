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
rz(2.958137) q[0];
sx q[0];
rz(-2.3975211) q[0];
sx q[0];
rz(-2.0896572) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(2.5084578) q[1];
sx q[1];
rz(8.8517744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6010701) q[0];
sx q[0];
rz(-0.94694505) q[0];
sx q[0];
rz(0.076657587) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0016167) q[2];
sx q[2];
rz(-0.9316906) q[2];
sx q[2];
rz(0.72102816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0172559) q[1];
sx q[1];
rz(-1.6410429) q[1];
sx q[1];
rz(2.6439366) q[1];
rz(-pi) q[2];
rz(-2.4450847) q[3];
sx q[3];
rz(-1.6497496) q[3];
sx q[3];
rz(2.1895308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8594592) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(-1.0563043) q[2];
rz(-0.022196444) q[3];
sx q[3];
rz(-0.76265097) q[3];
sx q[3];
rz(-1.4200042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9098772) q[0];
sx q[0];
rz(-0.48119369) q[0];
sx q[0];
rz(-0.43014446) q[0];
rz(0.12769708) q[1];
sx q[1];
rz(-1.9857429) q[1];
sx q[1];
rz(1.7040303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50629726) q[0];
sx q[0];
rz(-2.7397635) q[0];
sx q[0];
rz(-2.4319629) q[0];
rz(0.84354894) q[2];
sx q[2];
rz(-3.0650716) q[2];
sx q[2];
rz(0.25329548) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.592652) q[1];
sx q[1];
rz(-2.2363064) q[1];
sx q[1];
rz(-1.1489465) q[1];
x q[2];
rz(-0.53967472) q[3];
sx q[3];
rz(-1.7575348) q[3];
sx q[3];
rz(2.5591171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0251856) q[2];
sx q[2];
rz(-1.6802639) q[2];
sx q[2];
rz(-0.44542584) q[2];
rz(0.40924254) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(-0.87944952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5889848) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(-2.4267922) q[0];
rz(-2.0843166) q[1];
sx q[1];
rz(-0.42066586) q[1];
sx q[1];
rz(-0.32726273) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48178534) q[0];
sx q[0];
rz(-2.5707173) q[0];
sx q[0];
rz(-1.2587738) q[0];
x q[1];
rz(2.1959841) q[2];
sx q[2];
rz(-0.66788061) q[2];
sx q[2];
rz(-0.57648522) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.93720528) q[1];
sx q[1];
rz(-1.4193077) q[1];
sx q[1];
rz(-0.020843055) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25423519) q[3];
sx q[3];
rz(-1.2336858) q[3];
sx q[3];
rz(-1.1439307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0032234) q[2];
sx q[2];
rz(-2.1681163) q[2];
sx q[2];
rz(1.6376015) q[2];
rz(1.2184294) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(-2.159582) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686491) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(-2.9856227) q[0];
rz(-2.7929557) q[1];
sx q[1];
rz(-0.60395423) q[1];
sx q[1];
rz(1.2329996) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5208682) q[0];
sx q[0];
rz(-1.3107131) q[0];
sx q[0];
rz(-0.12518945) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21141036) q[2];
sx q[2];
rz(-0.93358913) q[2];
sx q[2];
rz(1.5508955) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8914596) q[1];
sx q[1];
rz(-0.95055184) q[1];
sx q[1];
rz(2.7930082) q[1];
rz(-pi) q[2];
rz(1.5855012) q[3];
sx q[3];
rz(-1.8222408) q[3];
sx q[3];
rz(-0.78212839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6984581) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(-1.6301463) q[2];
rz(2.002142) q[3];
sx q[3];
rz(-1.8099433) q[3];
sx q[3];
rz(0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646506) q[0];
sx q[0];
rz(-0.41780892) q[0];
sx q[0];
rz(-1.1530217) q[0];
rz(-2.5979089) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(0.26184729) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9898997) q[0];
sx q[0];
rz(-1.4873056) q[0];
sx q[0];
rz(-0.033647353) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33436454) q[2];
sx q[2];
rz(-1.0131256) q[2];
sx q[2];
rz(2.0206491) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0386069) q[1];
sx q[1];
rz(-0.9894045) q[1];
sx q[1];
rz(1.8379148) q[1];
rz(0.20528593) q[3];
sx q[3];
rz(-2.3488099) q[3];
sx q[3];
rz(-1.0012116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61949817) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(1.4052793) q[2];
rz(0.74553472) q[3];
sx q[3];
rz(-1.2657974) q[3];
sx q[3];
rz(-0.94329992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0010506823) q[0];
sx q[0];
rz(-0.11740919) q[0];
sx q[0];
rz(2.7897799) q[0];
rz(-0.44031269) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(2.9702759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85799828) q[0];
sx q[0];
rz(-2.4263315) q[0];
sx q[0];
rz(-0.37566988) q[0];
rz(-pi) q[1];
rz(-2.8181658) q[2];
sx q[2];
rz(-2.3617871) q[2];
sx q[2];
rz(-0.16743539) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80444634) q[1];
sx q[1];
rz(-2.4427572) q[1];
sx q[1];
rz(-0.032445907) q[1];
rz(-pi) q[2];
rz(0.8796575) q[3];
sx q[3];
rz(-2.1591957) q[3];
sx q[3];
rz(-3.1373051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.064284023) q[2];
sx q[2];
rz(-1.372154) q[2];
sx q[2];
rz(-2.1997931) q[2];
rz(1.1549548) q[3];
sx q[3];
rz(-2.933511) q[3];
sx q[3];
rz(1.340516) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6963541) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(3.136694) q[0];
rz(2.694963) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(0.68797025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82179994) q[0];
sx q[0];
rz(-1.3988136) q[0];
sx q[0];
rz(-2.4235307) q[0];
rz(-pi) q[1];
rz(1.2254459) q[2];
sx q[2];
rz(-0.98042578) q[2];
sx q[2];
rz(-1.8777868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35536218) q[1];
sx q[1];
rz(-1.7837824) q[1];
sx q[1];
rz(1.1287685) q[1];
rz(-pi) q[2];
rz(-0.68617188) q[3];
sx q[3];
rz(-2.1451575) q[3];
sx q[3];
rz(-2.2597093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61116162) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(-1.9200578) q[3];
sx q[3];
rz(-1.113021) q[3];
sx q[3];
rz(-0.59337029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4596443) q[0];
sx q[0];
rz(-1.3197897) q[0];
sx q[0];
rz(-0.090106877) q[0];
rz(-1.4247165) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(0.43509126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97532192) q[0];
sx q[0];
rz(-1.7164125) q[0];
sx q[0];
rz(0.50384371) q[0];
x q[1];
rz(-2.4736011) q[2];
sx q[2];
rz(-1.4282303) q[2];
sx q[2];
rz(3.0270456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7104946) q[1];
sx q[1];
rz(-2.5113547) q[1];
sx q[1];
rz(-1.789854) q[1];
rz(-pi) q[2];
rz(2.7172312) q[3];
sx q[3];
rz(-1.5699982) q[3];
sx q[3];
rz(-0.80959807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4789751) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(-2.5153861) q[2];
rz(2.9446757) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.7530493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58186746) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(0.054280601) q[0];
rz(-0.42298969) q[1];
sx q[1];
rz(-1.7661679) q[1];
sx q[1];
rz(1.8064226) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2535219) q[0];
sx q[0];
rz(-2.1162863) q[0];
sx q[0];
rz(-0.92407273) q[0];
rz(1.0112052) q[2];
sx q[2];
rz(-1.2647243) q[2];
sx q[2];
rz(-1.0572421) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74354913) q[1];
sx q[1];
rz(-0.7760007) q[1];
sx q[1];
rz(2.3659124) q[1];
rz(-pi) q[2];
rz(1.7087206) q[3];
sx q[3];
rz(-1.63878) q[3];
sx q[3];
rz(1.2757511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6173031) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(-0.35150251) q[2];
rz(-1.3737804) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(-0.26941776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7518625) q[0];
sx q[0];
rz(-0.26079145) q[0];
sx q[0];
rz(-2.3824298) q[0];
rz(2.0579386) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(1.0677451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6747492) q[0];
sx q[0];
rz(-0.7848133) q[0];
sx q[0];
rz(-0.87133566) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2944417) q[2];
sx q[2];
rz(-1.3047555) q[2];
sx q[2];
rz(0.58318116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9107161) q[1];
sx q[1];
rz(-1.9741804) q[1];
sx q[1];
rz(0.97539263) q[1];
x q[2];
rz(1.178252) q[3];
sx q[3];
rz(-2.4053229) q[3];
sx q[3];
rz(-2.6047849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7119673) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(-2.9313226) q[2];
rz(0.78768864) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(2.6113966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44337153) q[0];
sx q[0];
rz(-0.5589232) q[0];
sx q[0];
rz(1.9236175) q[0];
rz(1.5086077) q[1];
sx q[1];
rz(-1.2651545) q[1];
sx q[1];
rz(-1.9427585) q[1];
rz(-0.40661033) q[2];
sx q[2];
rz(-2.2175773) q[2];
sx q[2];
rz(-2.8178136) q[2];
rz(2.0785594) q[3];
sx q[3];
rz(-2.8294143) q[3];
sx q[3];
rz(-2.1635273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
