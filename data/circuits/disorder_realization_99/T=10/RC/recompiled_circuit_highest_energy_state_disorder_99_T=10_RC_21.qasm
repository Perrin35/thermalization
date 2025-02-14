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
rz(2.5853152) q[0];
sx q[0];
rz(-3.102432) q[0];
sx q[0];
rz(-2.477159) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(-1.4540949) q[1];
sx q[1];
rz(1.7323642) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8614757) q[0];
sx q[0];
rz(-1.7371157) q[0];
sx q[0];
rz(-2.6013199) q[0];
x q[1];
rz(-2.514226) q[2];
sx q[2];
rz(-0.53178501) q[2];
sx q[2];
rz(0.58284679) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4709027) q[1];
sx q[1];
rz(-1.2539314) q[1];
sx q[1];
rz(2.5554727) q[1];
x q[2];
rz(1.9326747) q[3];
sx q[3];
rz(-1.3986949) q[3];
sx q[3];
rz(-0.73735305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0613784) q[2];
sx q[2];
rz(-0.55157101) q[2];
sx q[2];
rz(0.75458327) q[2];
rz(1.4590229) q[3];
sx q[3];
rz(-2.2857234) q[3];
sx q[3];
rz(1.4065019) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.215312) q[0];
sx q[0];
rz(-0.54406852) q[0];
sx q[0];
rz(1.9790443) q[0];
rz(-0.7112208) q[1];
sx q[1];
rz(-0.7178719) q[1];
sx q[1];
rz(0.1444764) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1294589) q[0];
sx q[0];
rz(-0.30480614) q[0];
sx q[0];
rz(-0.27979677) q[0];
rz(-pi) q[1];
rz(0.8501022) q[2];
sx q[2];
rz(-1.0689702) q[2];
sx q[2];
rz(-0.36862954) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6934168) q[1];
sx q[1];
rz(-0.87478335) q[1];
sx q[1];
rz(-1.0820746) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8322093) q[3];
sx q[3];
rz(-2.2094634) q[3];
sx q[3];
rz(-0.84413278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7084536) q[2];
sx q[2];
rz(-1.837156) q[2];
sx q[2];
rz(2.8348095) q[2];
rz(2.3661546) q[3];
sx q[3];
rz(-0.38475761) q[3];
sx q[3];
rz(-0.14498372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6590092) q[0];
sx q[0];
rz(-0.44335303) q[0];
sx q[0];
rz(2.3495667) q[0];
rz(1.5576942) q[1];
sx q[1];
rz(-1.3521103) q[1];
sx q[1];
rz(2.1477594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64354831) q[0];
sx q[0];
rz(-2.1442765) q[0];
sx q[0];
rz(-1.9924966) q[0];
rz(1.6924528) q[2];
sx q[2];
rz(-0.98441974) q[2];
sx q[2];
rz(-0.77285484) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.12056226) q[1];
sx q[1];
rz(-1.1397727) q[1];
sx q[1];
rz(-0.31224183) q[1];
rz(-pi) q[2];
rz(-1.1937856) q[3];
sx q[3];
rz(-2.6600231) q[3];
sx q[3];
rz(-1.4762312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1899679) q[2];
sx q[2];
rz(-0.8364532) q[2];
sx q[2];
rz(2.9659029) q[2];
rz(2.9133993) q[3];
sx q[3];
rz(-2.5722645) q[3];
sx q[3];
rz(0.5051676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.93155414) q[0];
sx q[0];
rz(-0.84489548) q[0];
sx q[0];
rz(1.5616052) q[0];
rz(2.2700894) q[1];
sx q[1];
rz(-1.5305287) q[1];
sx q[1];
rz(-2.9759488) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4145395) q[0];
sx q[0];
rz(-3.139608) q[0];
sx q[0];
rz(0.63568799) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0087426337) q[2];
sx q[2];
rz(-1.1251108) q[2];
sx q[2];
rz(-0.937337) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45446515) q[1];
sx q[1];
rz(-1.6491967) q[1];
sx q[1];
rz(-2.3509174) q[1];
rz(-pi) q[2];
rz(0.40940855) q[3];
sx q[3];
rz(-2.5207075) q[3];
sx q[3];
rz(1.3001022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6988301) q[2];
sx q[2];
rz(-0.58219588) q[2];
sx q[2];
rz(-2.3637135) q[2];
rz(-0.88464087) q[3];
sx q[3];
rz(-1.9886465) q[3];
sx q[3];
rz(2.1625429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1409461) q[0];
sx q[0];
rz(-1.0362754) q[0];
sx q[0];
rz(-0.86877862) q[0];
rz(-0.90256214) q[1];
sx q[1];
rz(-1.2259918) q[1];
sx q[1];
rz(-0.57194078) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2615271) q[0];
sx q[0];
rz(-0.91996096) q[0];
sx q[0];
rz(-0.61246792) q[0];
rz(1.1210905) q[2];
sx q[2];
rz(-2.0313946) q[2];
sx q[2];
rz(3.1267303) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1871698) q[1];
sx q[1];
rz(-2.5989418) q[1];
sx q[1];
rz(-2.1869248) q[1];
rz(-1.192351) q[3];
sx q[3];
rz(-1.3406383) q[3];
sx q[3];
rz(-2.5251021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14704412) q[2];
sx q[2];
rz(-2.2037555) q[2];
sx q[2];
rz(-0.63329548) q[2];
rz(2.5607732) q[3];
sx q[3];
rz(-1.3588901) q[3];
sx q[3];
rz(0.66494989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13034114) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(-1.3379541) q[0];
rz(1.411571) q[1];
sx q[1];
rz(-2.7346225) q[1];
sx q[1];
rz(-2.3445047) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0397698) q[0];
sx q[0];
rz(-1.8443786) q[0];
sx q[0];
rz(0.16353971) q[0];
x q[1];
rz(3.0637805) q[2];
sx q[2];
rz(-1.0065662) q[2];
sx q[2];
rz(2.5112453) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89576777) q[1];
sx q[1];
rz(-1.3029769) q[1];
sx q[1];
rz(2.9755249) q[1];
rz(-pi) q[2];
rz(-1.8800297) q[3];
sx q[3];
rz(-0.54094523) q[3];
sx q[3];
rz(-0.23829392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7756614) q[2];
sx q[2];
rz(-2.4727827) q[2];
sx q[2];
rz(-2.456341) q[2];
rz(-0.25964409) q[3];
sx q[3];
rz(-2.1681867) q[3];
sx q[3];
rz(1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.892136) q[0];
sx q[0];
rz(-0.15661713) q[0];
sx q[0];
rz(-2.6020965) q[0];
rz(2.9501713) q[1];
sx q[1];
rz(-1.6554662) q[1];
sx q[1];
rz(2.6991381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060613077) q[0];
sx q[0];
rz(-0.0091305841) q[0];
sx q[0];
rz(-1.5332444) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94687825) q[2];
sx q[2];
rz(-1.2420474) q[2];
sx q[2];
rz(0.88916949) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.48414184) q[1];
sx q[1];
rz(-1.2353578) q[1];
sx q[1];
rz(0.92745933) q[1];
rz(-pi) q[2];
rz(2.4367853) q[3];
sx q[3];
rz(-0.92105812) q[3];
sx q[3];
rz(2.5832502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.847984) q[2];
sx q[2];
rz(-2.831735) q[2];
sx q[2];
rz(0.7027182) q[2];
rz(2.8637049) q[3];
sx q[3];
rz(-2.0746456) q[3];
sx q[3];
rz(-1.8028629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8948995) q[0];
sx q[0];
rz(-2.2791857) q[0];
sx q[0];
rz(2.605751) q[0];
rz(-0.72227532) q[1];
sx q[1];
rz(-1.4403789) q[1];
sx q[1];
rz(-0.78295082) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0412126) q[0];
sx q[0];
rz(-1.6579275) q[0];
sx q[0];
rz(1.9308596) q[0];
rz(-0.39733072) q[2];
sx q[2];
rz(-1.6226007) q[2];
sx q[2];
rz(-2.0837633) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47470505) q[1];
sx q[1];
rz(-2.0829446) q[1];
sx q[1];
rz(-2.1938503) q[1];
x q[2];
rz(-0.43704982) q[3];
sx q[3];
rz(-2.0756222) q[3];
sx q[3];
rz(-2.7849891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7213514) q[2];
sx q[2];
rz(-0.46478096) q[2];
sx q[2];
rz(2.6381524) q[2];
rz(2.8069046) q[3];
sx q[3];
rz(-1.8036333) q[3];
sx q[3];
rz(-0.23505841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6675785) q[0];
sx q[0];
rz(-0.36822167) q[0];
sx q[0];
rz(2.0817122) q[0];
rz(0.31479752) q[1];
sx q[1];
rz(-1.7960725) q[1];
sx q[1];
rz(1.5816636) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49071872) q[0];
sx q[0];
rz(-1.3482674) q[0];
sx q[0];
rz(-0.68993802) q[0];
rz(-pi) q[1];
x q[1];
rz(1.306172) q[2];
sx q[2];
rz(-1.6034063) q[2];
sx q[2];
rz(-2.5196241) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25645721) q[1];
sx q[1];
rz(-0.68084913) q[1];
sx q[1];
rz(1.225994) q[1];
rz(-pi) q[2];
rz(1.9859953) q[3];
sx q[3];
rz(-1.5506107) q[3];
sx q[3];
rz(-2.2184499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.18610893) q[2];
sx q[2];
rz(-1.9766221) q[2];
sx q[2];
rz(-2.1060409) q[2];
rz(-1.1459076) q[3];
sx q[3];
rz(-1.112554) q[3];
sx q[3];
rz(-0.15315332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.938852) q[0];
sx q[0];
rz(-3.0227612) q[0];
sx q[0];
rz(2.3760997) q[0];
rz(2.1047523) q[1];
sx q[1];
rz(-1.2988657) q[1];
sx q[1];
rz(-3.0293005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0033007) q[0];
sx q[0];
rz(-2.0202185) q[0];
sx q[0];
rz(1.9339191) q[0];
rz(-pi) q[1];
rz(-2.3621777) q[2];
sx q[2];
rz(-1.1128491) q[2];
sx q[2];
rz(-3.0061257) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0794509) q[1];
sx q[1];
rz(-1.7711258) q[1];
sx q[1];
rz(2.315534) q[1];
x q[2];
rz(3.0946135) q[3];
sx q[3];
rz(-2.3571797) q[3];
sx q[3];
rz(1.3734773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74467337) q[2];
sx q[2];
rz(-1.2770709) q[2];
sx q[2];
rz(0.12167682) q[2];
rz(2.7820898) q[3];
sx q[3];
rz(-2.4183141) q[3];
sx q[3];
rz(-0.017688964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59035463) q[0];
sx q[0];
rz(-1.4689162) q[0];
sx q[0];
rz(1.3015251) q[0];
rz(1.7474668) q[1];
sx q[1];
rz(-1.9448517) q[1];
sx q[1];
rz(1.3762884) q[1];
rz(1.1375221) q[2];
sx q[2];
rz(-0.99356298) q[2];
sx q[2];
rz(-2.1504924) q[2];
rz(-1.1351552) q[3];
sx q[3];
rz(-1.451685) q[3];
sx q[3];
rz(1.4901037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
