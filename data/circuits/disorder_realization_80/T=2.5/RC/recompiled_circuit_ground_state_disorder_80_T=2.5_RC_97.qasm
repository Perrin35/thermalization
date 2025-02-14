OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.03117938) q[0];
sx q[0];
rz(2.2194982) q[0];
sx q[0];
rz(10.194869) q[0];
rz(0.71647477) q[1];
sx q[1];
rz(-2.1991576) q[1];
sx q[1];
rz(0.64185774) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19260423) q[0];
sx q[0];
rz(-1.7718713) q[0];
sx q[0];
rz(-1.4549535) q[0];
rz(-pi) q[1];
rz(-1.5719169) q[2];
sx q[2];
rz(-1.5693451) q[2];
sx q[2];
rz(0.077479428) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9145467) q[1];
sx q[1];
rz(-1.7265336) q[1];
sx q[1];
rz(-2.7776633) q[1];
rz(0.58310469) q[3];
sx q[3];
rz(-0.92188406) q[3];
sx q[3];
rz(0.16498868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98051071) q[2];
sx q[2];
rz(-2.0445721) q[2];
sx q[2];
rz(-2.3316627) q[2];
rz(3.1071281) q[3];
sx q[3];
rz(-0.66027111) q[3];
sx q[3];
rz(0.1035498) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.506839) q[0];
sx q[0];
rz(-2.6320808) q[0];
sx q[0];
rz(-2.4822045) q[0];
rz(-1.7022645) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(2.6606681) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56136614) q[0];
sx q[0];
rz(-0.49005383) q[0];
sx q[0];
rz(-1.1183764) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85428736) q[2];
sx q[2];
rz(-2.169974) q[2];
sx q[2];
rz(-0.5391268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2175423) q[1];
sx q[1];
rz(-2.9311507) q[1];
sx q[1];
rz(1.4242875) q[1];
rz(-pi) q[2];
rz(2.0388076) q[3];
sx q[3];
rz(-1.0429405) q[3];
sx q[3];
rz(-0.46609391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3769569) q[2];
sx q[2];
rz(-0.83676338) q[2];
sx q[2];
rz(-0.62977201) q[2];
rz(1.9624286) q[3];
sx q[3];
rz(-2.4383014) q[3];
sx q[3];
rz(-2.0298957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4513007) q[0];
sx q[0];
rz(-3.1307104) q[0];
sx q[0];
rz(-0.96845281) q[0];
rz(0.14006607) q[1];
sx q[1];
rz(-1.3532956) q[1];
sx q[1];
rz(-2.5586939) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2067277) q[0];
sx q[0];
rz(-1.4040274) q[0];
sx q[0];
rz(1.5251446) q[0];
rz(-pi) q[1];
rz(2.3034322) q[2];
sx q[2];
rz(-1.3930905) q[2];
sx q[2];
rz(-0.017824307) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1000468) q[1];
sx q[1];
rz(-2.1963781) q[1];
sx q[1];
rz(2.2872674) q[1];
x q[2];
rz(-2.8606877) q[3];
sx q[3];
rz(-1.2634522) q[3];
sx q[3];
rz(2.8859649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53050238) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(-0.21406847) q[2];
rz(2.8392082) q[3];
sx q[3];
rz(-2.3979135) q[3];
sx q[3];
rz(2.7157057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9532303) q[0];
sx q[0];
rz(-1.0091877) q[0];
sx q[0];
rz(2.0971712) q[0];
rz(-2.2638679) q[1];
sx q[1];
rz(-1.5889771) q[1];
sx q[1];
rz(0.16876076) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8498189) q[0];
sx q[0];
rz(-2.5892604) q[0];
sx q[0];
rz(-2.0122347) q[0];
rz(-pi) q[1];
rz(1.243882) q[2];
sx q[2];
rz(-1.3762646) q[2];
sx q[2];
rz(0.036613001) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16774878) q[1];
sx q[1];
rz(-2.1653321) q[1];
sx q[1];
rz(1.3092625) q[1];
x q[2];
rz(-2.2764858) q[3];
sx q[3];
rz(-1.5388515) q[3];
sx q[3];
rz(0.77199304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3622482) q[2];
sx q[2];
rz(-1.2835953) q[2];
sx q[2];
rz(2.0812422) q[2];
rz(-2.849071) q[3];
sx q[3];
rz(-0.7258324) q[3];
sx q[3];
rz(3.0574851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.5523858) q[0];
sx q[0];
rz(-2.4942106) q[0];
sx q[0];
rz(0.86828434) q[0];
rz(0.6262511) q[1];
sx q[1];
rz(-1.7528844) q[1];
sx q[1];
rz(-1.7832322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043644301) q[0];
sx q[0];
rz(-1.4820423) q[0];
sx q[0];
rz(2.8525145) q[0];
rz(0.50164671) q[2];
sx q[2];
rz(-0.42914647) q[2];
sx q[2];
rz(2.8535064) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89857453) q[1];
sx q[1];
rz(-2.822091) q[1];
sx q[1];
rz(0.77844365) q[1];
rz(0.63774469) q[3];
sx q[3];
rz(-1.2079835) q[3];
sx q[3];
rz(0.90076288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2061578) q[2];
sx q[2];
rz(-0.81659603) q[2];
sx q[2];
rz(2.6182776) q[2];
rz(2.7715136) q[3];
sx q[3];
rz(-2.3612634) q[3];
sx q[3];
rz(-1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0321781) q[0];
sx q[0];
rz(-0.21571708) q[0];
sx q[0];
rz(-2.7874462) q[0];
rz(-2.1996563) q[1];
sx q[1];
rz(-1.4553921) q[1];
sx q[1];
rz(1.8249493) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9958651) q[0];
sx q[0];
rz(-2.6168129) q[0];
sx q[0];
rz(-2.9324233) q[0];
x q[1];
rz(-3.1194341) q[2];
sx q[2];
rz(-0.54931927) q[2];
sx q[2];
rz(1.4469128) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.020654708) q[1];
sx q[1];
rz(-1.5096501) q[1];
sx q[1];
rz(2.967359) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8527669) q[3];
sx q[3];
rz(-2.1957955) q[3];
sx q[3];
rz(0.19240141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3256623) q[2];
sx q[2];
rz(-0.38620913) q[2];
sx q[2];
rz(-0.33829921) q[2];
rz(2.6650688) q[3];
sx q[3];
rz(-0.75348133) q[3];
sx q[3];
rz(-2.8009955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7020096) q[0];
sx q[0];
rz(-1.5583353) q[0];
sx q[0];
rz(2.6038792) q[0];
rz(1.733755) q[1];
sx q[1];
rz(-2.6815963) q[1];
sx q[1];
rz(2.5163311) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940373) q[0];
sx q[0];
rz(-0.68935822) q[0];
sx q[0];
rz(2.5563142) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4971259) q[2];
sx q[2];
rz(-1.0331312) q[2];
sx q[2];
rz(-1.0755838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3846674) q[1];
sx q[1];
rz(-2.1672492) q[1];
sx q[1];
rz(-1.3104964) q[1];
rz(2.7688249) q[3];
sx q[3];
rz(-2.3298752) q[3];
sx q[3];
rz(1.2889287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9739146) q[2];
sx q[2];
rz(-2.6748071) q[2];
sx q[2];
rz(-1.7612877) q[2];
rz(-2.7050833) q[3];
sx q[3];
rz(-2.1197539) q[3];
sx q[3];
rz(0.71353394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1768782) q[0];
sx q[0];
rz(-0.62129337) q[0];
sx q[0];
rz(0.51625133) q[0];
rz(0.77975726) q[1];
sx q[1];
rz(-2.1596491) q[1];
sx q[1];
rz(1.0468743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8314227) q[0];
sx q[0];
rz(-2.8084233) q[0];
sx q[0];
rz(-2.6875671) q[0];
rz(-pi) q[1];
rz(0.5811695) q[2];
sx q[2];
rz(-2.0072674) q[2];
sx q[2];
rz(1.798686) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.13345737) q[1];
sx q[1];
rz(-1.5127458) q[1];
sx q[1];
rz(-1.4769555) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47886301) q[3];
sx q[3];
rz(-1.2148464) q[3];
sx q[3];
rz(1.009481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9376935) q[2];
sx q[2];
rz(-2.9475309) q[2];
sx q[2];
rz(0.95721179) q[2];
rz(0.29414487) q[3];
sx q[3];
rz(-1.1295986) q[3];
sx q[3];
rz(-2.8637776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99943632) q[0];
sx q[0];
rz(-0.23040982) q[0];
sx q[0];
rz(0.18705046) q[0];
rz(0.43141463) q[1];
sx q[1];
rz(-0.43771935) q[1];
sx q[1];
rz(-1.4923219) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6607912) q[0];
sx q[0];
rz(-1.4993164) q[0];
sx q[0];
rz(-3.0707703) q[0];
rz(1.9052986) q[2];
sx q[2];
rz(-1.750947) q[2];
sx q[2];
rz(-2.1366773) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8111251) q[1];
sx q[1];
rz(-2.6973675) q[1];
sx q[1];
rz(-1.1796239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1077376) q[3];
sx q[3];
rz(-0.18165043) q[3];
sx q[3];
rz(-2.8665115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6734068) q[2];
sx q[2];
rz(-0.91172051) q[2];
sx q[2];
rz(-2.1569596) q[2];
rz(0.51472384) q[3];
sx q[3];
rz(-0.52557164) q[3];
sx q[3];
rz(-1.908186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.89727) q[0];
sx q[0];
rz(-1.6615302) q[0];
sx q[0];
rz(0.75862128) q[0];
rz(1.9482535) q[1];
sx q[1];
rz(-1.1921644) q[1];
sx q[1];
rz(1.4512482) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8817026) q[0];
sx q[0];
rz(-0.30891793) q[0];
sx q[0];
rz(2.3007459) q[0];
rz(-pi) q[1];
rz(-0.042165857) q[2];
sx q[2];
rz(-1.3683967) q[2];
sx q[2];
rz(0.037994904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2047233) q[1];
sx q[1];
rz(-1.544049) q[1];
sx q[1];
rz(2.1863947) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2281353) q[3];
sx q[3];
rz(-2.6413915) q[3];
sx q[3];
rz(0.72266912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.54185581) q[2];
sx q[2];
rz(-2.2208322) q[2];
sx q[2];
rz(0.80121458) q[2];
rz(2.2090705) q[3];
sx q[3];
rz(-1.2152117) q[3];
sx q[3];
rz(-2.7264989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5781317) q[0];
sx q[0];
rz(-1.7628071) q[0];
sx q[0];
rz(2.5264869) q[0];
rz(0.32147944) q[1];
sx q[1];
rz(-2.165806) q[1];
sx q[1];
rz(1.4478366) q[1];
rz(-0.64340016) q[2];
sx q[2];
rz(-2.9360124) q[2];
sx q[2];
rz(0.58694466) q[2];
rz(0.23033167) q[3];
sx q[3];
rz(-1.2970222) q[3];
sx q[3];
rz(2.6877689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
