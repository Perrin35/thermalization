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
rz(-0.55627745) q[0];
sx q[0];
rz(-0.039160691) q[0];
sx q[0];
rz(-0.66443366) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(-1.4540949) q[1];
sx q[1];
rz(1.7323642) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7519386) q[0];
sx q[0];
rz(-2.1028127) q[0];
sx q[0];
rz(-1.3774894) q[0];
rz(1.9033236) q[2];
sx q[2];
rz(-1.1477787) q[2];
sx q[2];
rz(0.11655434) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8027026) q[1];
sx q[1];
rz(-2.4842508) q[1];
sx q[1];
rz(0.53513066) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9326747) q[3];
sx q[3];
rz(-1.3986949) q[3];
sx q[3];
rz(2.4042396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0802143) q[2];
sx q[2];
rz(-0.55157101) q[2];
sx q[2];
rz(2.3870094) q[2];
rz(1.4590229) q[3];
sx q[3];
rz(-0.85586923) q[3];
sx q[3];
rz(-1.4065019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9262806) q[0];
sx q[0];
rz(-2.5975241) q[0];
sx q[0];
rz(-1.1625483) q[0];
rz(-2.4303719) q[1];
sx q[1];
rz(-2.4237207) q[1];
sx q[1];
rz(-2.9971163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.719559) q[0];
sx q[0];
rz(-1.2782017) q[0];
sx q[0];
rz(1.6574615) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63074449) q[2];
sx q[2];
rz(-2.1877383) q[2];
sx q[2];
rz(1.601905) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4481758) q[1];
sx q[1];
rz(-0.87478335) q[1];
sx q[1];
rz(1.0820746) q[1];
rz(2.8322093) q[3];
sx q[3];
rz(-0.9321292) q[3];
sx q[3];
rz(2.2974599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7084536) q[2];
sx q[2];
rz(-1.3044367) q[2];
sx q[2];
rz(-0.30678314) q[2];
rz(-0.77543801) q[3];
sx q[3];
rz(-2.756835) q[3];
sx q[3];
rz(-2.9966089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48258346) q[0];
sx q[0];
rz(-0.44335303) q[0];
sx q[0];
rz(-0.79202598) q[0];
rz(-1.5576942) q[1];
sx q[1];
rz(-1.3521103) q[1];
sx q[1];
rz(-2.1477594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64354831) q[0];
sx q[0];
rz(-2.1442765) q[0];
sx q[0];
rz(1.9924966) q[0];
x q[1];
rz(-2.5517919) q[2];
sx q[2];
rz(-1.6720534) q[2];
sx q[2];
rz(0.73039215) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5842918) q[1];
sx q[1];
rz(-1.2879432) q[1];
sx q[1];
rz(-1.1206579) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.951521) q[3];
sx q[3];
rz(-1.1255923) q[3];
sx q[3];
rz(-1.0560738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1899679) q[2];
sx q[2];
rz(-2.3051395) q[2];
sx q[2];
rz(-0.17568976) q[2];
rz(2.9133993) q[3];
sx q[3];
rz(-2.5722645) q[3];
sx q[3];
rz(0.5051676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2100385) q[0];
sx q[0];
rz(-0.84489548) q[0];
sx q[0];
rz(1.5616052) q[0];
rz(2.2700894) q[1];
sx q[1];
rz(-1.611064) q[1];
sx q[1];
rz(2.9759488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3627421) q[0];
sx q[0];
rz(-1.5691994) q[0];
sx q[0];
rz(-1.5719747) q[0];
rz(-pi) q[1];
rz(2.0164967) q[2];
sx q[2];
rz(-1.5786849) q[2];
sx q[2];
rz(-2.5119022) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1042479) q[1];
sx q[1];
rz(-2.3583721) q[1];
sx q[1];
rz(1.4595637) q[1];
rz(0.40940855) q[3];
sx q[3];
rz(-2.5207075) q[3];
sx q[3];
rz(1.3001022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44276253) q[2];
sx q[2];
rz(-2.5593968) q[2];
sx q[2];
rz(-0.77787918) q[2];
rz(2.2569518) q[3];
sx q[3];
rz(-1.1529461) q[3];
sx q[3];
rz(0.9790498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1409461) q[0];
sx q[0];
rz(-2.1053173) q[0];
sx q[0];
rz(-0.86877862) q[0];
rz(-0.90256214) q[1];
sx q[1];
rz(-1.2259918) q[1];
sx q[1];
rz(2.5696519) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71170744) q[0];
sx q[0];
rz(-2.0458303) q[0];
sx q[0];
rz(0.82127969) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4221572) q[2];
sx q[2];
rz(-2.509382) q[2];
sx q[2];
rz(-2.3001463) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8780788) q[1];
sx q[1];
rz(-2.0058419) q[1];
sx q[1];
rz(0.33532354) q[1];
rz(-pi) q[2];
rz(1.9492416) q[3];
sx q[3];
rz(-1.8009543) q[3];
sx q[3];
rz(-0.6164906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.14704412) q[2];
sx q[2];
rz(-2.2037555) q[2];
sx q[2];
rz(-0.63329548) q[2];
rz(0.58081943) q[3];
sx q[3];
rz(-1.7827026) q[3];
sx q[3];
rz(0.66494989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(3.0112515) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(-1.8036386) q[0];
rz(-1.7300216) q[1];
sx q[1];
rz(-2.7346225) q[1];
sx q[1];
rz(0.79708797) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6501084) q[0];
sx q[0];
rz(-0.31768018) q[0];
sx q[0];
rz(-2.0965212) q[0];
rz(-1.448579) q[2];
sx q[2];
rz(-2.5725992) q[2];
sx q[2];
rz(0.48556604) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.63070272) q[1];
sx q[1];
rz(-1.4107009) q[1];
sx q[1];
rz(-1.2994205) q[1];
rz(-pi) q[2];
rz(-2.09054) q[3];
sx q[3];
rz(-1.4134348) q[3];
sx q[3];
rz(1.5998154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.36593124) q[2];
sx q[2];
rz(-0.66880995) q[2];
sx q[2];
rz(0.68525165) q[2];
rz(-0.25964409) q[3];
sx q[3];
rz(-0.97340596) q[3];
sx q[3];
rz(-1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.892136) q[0];
sx q[0];
rz(-0.15661713) q[0];
sx q[0];
rz(2.6020965) q[0];
rz(0.19142137) q[1];
sx q[1];
rz(-1.4861264) q[1];
sx q[1];
rz(-0.44245455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434261) q[0];
sx q[0];
rz(-1.5616722) q[0];
sx q[0];
rz(0.00034279963) q[0];
rz(-pi) q[1];
rz(-2.1947144) q[2];
sx q[2];
rz(-1.8995452) q[2];
sx q[2];
rz(0.88916949) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.328622) q[1];
sx q[1];
rz(-0.96862205) q[1];
sx q[1];
rz(-0.41090907) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4367853) q[3];
sx q[3];
rz(-0.92105812) q[3];
sx q[3];
rz(0.55834246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29360867) q[2];
sx q[2];
rz(-0.30985761) q[2];
sx q[2];
rz(2.4388745) q[2];
rz(-2.8637049) q[3];
sx q[3];
rz(-1.0669471) q[3];
sx q[3];
rz(-1.8028629) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24669312) q[0];
sx q[0];
rz(-2.2791857) q[0];
sx q[0];
rz(0.53584164) q[0];
rz(2.4193173) q[1];
sx q[1];
rz(-1.4403789) q[1];
sx q[1];
rz(2.3586418) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7039258) q[0];
sx q[0];
rz(-1.9294318) q[0];
sx q[0];
rz(0.093067972) q[0];
rz(3.0083904) q[2];
sx q[2];
rz(-2.7410772) q[2];
sx q[2];
rz(2.7513964) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7069426) q[1];
sx q[1];
rz(-2.1044113) q[1];
sx q[1];
rz(0.60551079) q[1];
x q[2];
rz(-2.1184741) q[3];
sx q[3];
rz(-1.9503106) q[3];
sx q[3];
rz(-1.7051769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7213514) q[2];
sx q[2];
rz(-2.6768117) q[2];
sx q[2];
rz(2.6381524) q[2];
rz(0.33468801) q[3];
sx q[3];
rz(-1.3379593) q[3];
sx q[3];
rz(2.9065342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47401416) q[0];
sx q[0];
rz(-2.773371) q[0];
sx q[0];
rz(1.0598805) q[0];
rz(2.8267951) q[1];
sx q[1];
rz(-1.3455201) q[1];
sx q[1];
rz(1.5816636) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81875694) q[0];
sx q[0];
rz(-2.4222932) q[0];
sx q[0];
rz(-0.34157217) q[0];
rz(-pi) q[1];
rz(3.1078075) q[2];
sx q[2];
rz(-1.8352766) q[2];
sx q[2];
rz(-2.1839301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.25645721) q[1];
sx q[1];
rz(-0.68084913) q[1];
sx q[1];
rz(1.225994) q[1];
rz(-pi) q[2];
rz(1.520789) q[3];
sx q[3];
rz(-2.7259318) q[3];
sx q[3];
rz(2.5396944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9554837) q[2];
sx q[2];
rz(-1.1649705) q[2];
sx q[2];
rz(1.0355518) q[2];
rz(-1.1459076) q[3];
sx q[3];
rz(-2.0290387) q[3];
sx q[3];
rz(0.15315332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2027407) q[0];
sx q[0];
rz(-0.11883141) q[0];
sx q[0];
rz(-2.3760997) q[0];
rz(-2.1047523) q[1];
sx q[1];
rz(-1.8427269) q[1];
sx q[1];
rz(-3.0293005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0033007) q[0];
sx q[0];
rz(-1.1213741) q[0];
sx q[0];
rz(-1.2076735) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77941497) q[2];
sx q[2];
rz(-2.0287435) q[2];
sx q[2];
rz(0.13546695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8143554) q[1];
sx q[1];
rz(-0.8443409) q[1];
sx q[1];
rz(2.8721456) q[1];
rz(-2.3577317) q[3];
sx q[3];
rz(-1.6039768) q[3];
sx q[3];
rz(-2.9110094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74467337) q[2];
sx q[2];
rz(-1.8645218) q[2];
sx q[2];
rz(3.0199158) q[2];
rz(0.35950288) q[3];
sx q[3];
rz(-0.72327852) q[3];
sx q[3];
rz(3.1239037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.551238) q[0];
sx q[0];
rz(-1.6726765) q[0];
sx q[0];
rz(-1.8400675) q[0];
rz(1.3941258) q[1];
sx q[1];
rz(-1.196741) q[1];
sx q[1];
rz(-1.7653042) q[1];
rz(2.0040705) q[2];
sx q[2];
rz(-2.1480297) q[2];
sx q[2];
rz(0.99110023) q[2];
rz(0.13124851) q[3];
sx q[3];
rz(-2.003142) q[3];
sx q[3];
rz(3.005645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
