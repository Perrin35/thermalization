OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8653712) q[0];
sx q[0];
rz(-2.2844391) q[0];
sx q[0];
rz(3.0091118) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(-0.58499709) q[1];
sx q[1];
rz(2.4490228) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24774691) q[0];
sx q[0];
rz(-2.6129299) q[0];
sx q[0];
rz(1.371944) q[0];
x q[1];
rz(-2.5508467) q[2];
sx q[2];
rz(-2.4060537) q[2];
sx q[2];
rz(0.22434805) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8323648) q[1];
sx q[1];
rz(-1.4432866) q[1];
sx q[1];
rz(-3.068919) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0971783) q[3];
sx q[3];
rz(-1.3485104) q[3];
sx q[3];
rz(-1.8412631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0341558) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.5365323) q[2];
rz(-1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-0.03604123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543106) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(1.8923627) q[0];
rz(2.5800887) q[1];
sx q[1];
rz(-2.3655472) q[1];
sx q[1];
rz(0.5805648) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3378355) q[0];
sx q[0];
rz(-0.32807402) q[0];
sx q[0];
rz(0.34123811) q[0];
rz(1.7891907) q[2];
sx q[2];
rz(-0.94422715) q[2];
sx q[2];
rz(-0.44109694) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3483352) q[1];
sx q[1];
rz(-2.1264592) q[1];
sx q[1];
rz(2.4130915) q[1];
rz(0.35238102) q[3];
sx q[3];
rz(-1.1412732) q[3];
sx q[3];
rz(-0.22892117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7291752) q[2];
sx q[2];
rz(-0.20755945) q[2];
sx q[2];
rz(-0.87835971) q[2];
rz(-2.7495524) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(-0.6033321) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47675258) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-2.3679249) q[0];
rz(3.1402918) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(0.032827854) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5884018) q[0];
sx q[0];
rz(-1.8817888) q[0];
sx q[0];
rz(1.3283967) q[0];
rz(-pi) q[1];
rz(-0.29067729) q[2];
sx q[2];
rz(-2.103984) q[2];
sx q[2];
rz(0.79613396) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58363885) q[1];
sx q[1];
rz(-1.4065521) q[1];
sx q[1];
rz(0.79973952) q[1];
rz(-pi) q[2];
rz(1.0333943) q[3];
sx q[3];
rz(-0.61363797) q[3];
sx q[3];
rz(-2.9463241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34439987) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(0.27734217) q[2];
rz(-0.39595655) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(-2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70401496) q[0];
sx q[0];
rz(-2.8058348) q[0];
sx q[0];
rz(0.26279703) q[0];
rz(2.908005) q[1];
sx q[1];
rz(-2.3065152) q[1];
sx q[1];
rz(2.3707726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9642826) q[0];
sx q[0];
rz(-1.5543803) q[0];
sx q[0];
rz(-1.55127) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0609444) q[2];
sx q[2];
rz(-1.66301) q[2];
sx q[2];
rz(0.25445081) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0499038) q[1];
sx q[1];
rz(-2.4871475) q[1];
sx q[1];
rz(2.5212538) q[1];
rz(2.2724857) q[3];
sx q[3];
rz(-2.8512555) q[3];
sx q[3];
rz(1.7770191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(2.4285994) q[2];
rz(-1.0130079) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7825496) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(-1.7011401) q[0];
rz(0.094093181) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(0.17000155) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1124135) q[0];
sx q[0];
rz(-2.3110483) q[0];
sx q[0];
rz(-1.889617) q[0];
rz(-pi) q[1];
rz(-2.0427809) q[2];
sx q[2];
rz(-1.3655647) q[2];
sx q[2];
rz(-0.99265487) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1483037) q[1];
sx q[1];
rz(-1.9699886) q[1];
sx q[1];
rz(-2.1972375) q[1];
rz(0.65808987) q[3];
sx q[3];
rz(-0.98539017) q[3];
sx q[3];
rz(0.64774367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1468982) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(-1.7374932) q[2];
rz(-1.5480301) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(-1.7061957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4858522) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(0.45853841) q[0];
rz(-0.25587747) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(-2.4564254) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2049094) q[0];
sx q[0];
rz(-0.65403599) q[0];
sx q[0];
rz(1.4564287) q[0];
rz(-0.14641996) q[2];
sx q[2];
rz(-2.358987) q[2];
sx q[2];
rz(0.27461068) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8243858) q[1];
sx q[1];
rz(-2.3134391) q[1];
sx q[1];
rz(-0.87888996) q[1];
rz(-pi) q[2];
rz(-0.018304304) q[3];
sx q[3];
rz(-0.98494512) q[3];
sx q[3];
rz(-2.626112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7489862) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(1.7123327) q[2];
rz(1.0990934) q[3];
sx q[3];
rz(-0.50656879) q[3];
sx q[3];
rz(-2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012506164) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(0.72934735) q[0];
rz(0.29306456) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(1.1475295) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.731819) q[0];
sx q[0];
rz(-1.496135) q[0];
sx q[0];
rz(-2.7274107) q[0];
rz(-pi) q[1];
rz(1.0853026) q[2];
sx q[2];
rz(-2.448423) q[2];
sx q[2];
rz(-1.8144516) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9784669) q[1];
sx q[1];
rz(-1.6118057) q[1];
sx q[1];
rz(0.48467111) q[1];
x q[2];
rz(-2.494874) q[3];
sx q[3];
rz(-2.5463856) q[3];
sx q[3];
rz(-1.7531542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3532233) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(0.74907556) q[2];
rz(0.64368147) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(0.914004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1483243) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(0.83129445) q[0];
rz(-1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(-2.7430699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7300028) q[0];
sx q[0];
rz(-1.411479) q[0];
sx q[0];
rz(2.0429862) q[0];
rz(-pi) q[1];
rz(2.5404846) q[2];
sx q[2];
rz(-1.483344) q[2];
sx q[2];
rz(0.89278883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41259137) q[1];
sx q[1];
rz(-1.593643) q[1];
sx q[1];
rz(-2.1680135) q[1];
rz(-2.1274444) q[3];
sx q[3];
rz(-2.017052) q[3];
sx q[3];
rz(2.2198912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1901671) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(1.8438967) q[2];
rz(1.1249582) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012638906) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(-1.3893611) q[0];
rz(1.6268436) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(-1.0983889) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59745715) q[0];
sx q[0];
rz(-1.8034593) q[0];
sx q[0];
rz(-1.0730037) q[0];
rz(-pi) q[1];
rz(2.056698) q[2];
sx q[2];
rz(-2.3294318) q[2];
sx q[2];
rz(-2.95129) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2213649) q[1];
sx q[1];
rz(-2.5623294) q[1];
sx q[1];
rz(3.1406162) q[1];
rz(-pi) q[2];
rz(-1.3994201) q[3];
sx q[3];
rz(-1.3203353) q[3];
sx q[3];
rz(2.7385538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2417458) q[2];
sx q[2];
rz(-1.8490303) q[2];
sx q[2];
rz(-0.51952726) q[2];
rz(1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(-1.3174723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7463503) q[0];
sx q[0];
rz(-2.0210176) q[0];
sx q[0];
rz(-0.46646068) q[0];
rz(0.17164224) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-0.62896532) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2254612) q[0];
sx q[0];
rz(-2.7636508) q[0];
sx q[0];
rz(1.4378689) q[0];
x q[1];
rz(1.8877108) q[2];
sx q[2];
rz(-2.4744518) q[2];
sx q[2];
rz(-2.4184879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4911463) q[1];
sx q[1];
rz(-1.568734) q[1];
sx q[1];
rz(2.6994929) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4217989) q[3];
sx q[3];
rz(-1.3070953) q[3];
sx q[3];
rz(-1.7566453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(-2.0824599) q[2];
rz(0.6774261) q[3];
sx q[3];
rz(-0.99223653) q[3];
sx q[3];
rz(-2.2476851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8582936) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(0.25390608) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(-1.614744) q[2];
sx q[2];
rz(-0.84779253) q[2];
sx q[2];
rz(0.070889125) q[2];
rz(-0.98388381) q[3];
sx q[3];
rz(-1.6279396) q[3];
sx q[3];
rz(2.770594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];