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
rz(2.477159) q[0];
rz(2.9593664) q[1];
sx q[1];
rz(-1.6874977) q[1];
sx q[1];
rz(1.4092285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8614757) q[0];
sx q[0];
rz(-1.7371157) q[0];
sx q[0];
rz(0.54027277) q[0];
rz(-pi) q[1];
rz(2.6970942) q[2];
sx q[2];
rz(-1.8730361) q[2];
sx q[2];
rz(-1.5465178) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8027026) q[1];
sx q[1];
rz(-0.65734184) q[1];
sx q[1];
rz(2.606462) q[1];
x q[2];
rz(-2.9578311) q[3];
sx q[3];
rz(-1.9270883) q[3];
sx q[3];
rz(-2.3728865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0613784) q[2];
sx q[2];
rz(-2.5900216) q[2];
sx q[2];
rz(2.3870094) q[2];
rz(1.4590229) q[3];
sx q[3];
rz(-2.2857234) q[3];
sx q[3];
rz(-1.7350908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.215312) q[0];
sx q[0];
rz(-0.54406852) q[0];
sx q[0];
rz(-1.1625483) q[0];
rz(-0.7112208) q[1];
sx q[1];
rz(-2.4237207) q[1];
sx q[1];
rz(-0.1444764) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677749) q[0];
sx q[0];
rz(-1.4878232) q[0];
sx q[0];
rz(-2.847958) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63074449) q[2];
sx q[2];
rz(-2.1877383) q[2];
sx q[2];
rz(-1.5396876) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.45119845) q[1];
sx q[1];
rz(-1.2022126) q[1];
sx q[1];
rz(-0.75779961) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30938332) q[3];
sx q[3];
rz(-2.2094634) q[3];
sx q[3];
rz(2.2974599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43313906) q[2];
sx q[2];
rz(-1.837156) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48258346) q[0];
sx q[0];
rz(-0.44335303) q[0];
sx q[0];
rz(0.79202598) q[0];
rz(1.5576942) q[1];
sx q[1];
rz(-1.7894824) q[1];
sx q[1];
rz(0.99383324) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3344468) q[0];
sx q[0];
rz(-0.69753555) q[0];
sx q[0];
rz(-0.56484449) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5517919) q[2];
sx q[2];
rz(-1.4695393) q[2];
sx q[2];
rz(-0.73039215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12056226) q[1];
sx q[1];
rz(-2.0018199) q[1];
sx q[1];
rz(-0.31224183) q[1];
rz(-2.023104) q[3];
sx q[3];
rz(-1.7421466) q[3];
sx q[3];
rz(2.7095344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9516248) q[2];
sx q[2];
rz(-0.8364532) q[2];
sx q[2];
rz(2.9659029) q[2];
rz(0.22819337) q[3];
sx q[3];
rz(-0.56932813) q[3];
sx q[3];
rz(0.5051676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2100385) q[0];
sx q[0];
rz(-0.84489548) q[0];
sx q[0];
rz(1.5616052) q[0];
rz(-0.87150323) q[1];
sx q[1];
rz(-1.611064) q[1];
sx q[1];
rz(-0.16564381) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77885056) q[0];
sx q[0];
rz(-1.5723933) q[0];
sx q[0];
rz(1.5719747) q[0];
rz(-pi) q[1];
rz(1.125096) q[2];
sx q[2];
rz(-1.5629077) q[2];
sx q[2];
rz(-2.5119022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6871275) q[1];
sx q[1];
rz(-1.6491967) q[1];
sx q[1];
rz(0.7906753) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40940855) q[3];
sx q[3];
rz(-2.5207075) q[3];
sx q[3];
rz(1.8414904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44276253) q[2];
sx q[2];
rz(-0.58219588) q[2];
sx q[2];
rz(-2.3637135) q[2];
rz(-2.2569518) q[3];
sx q[3];
rz(-1.1529461) q[3];
sx q[3];
rz(2.1625429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00064656249) q[0];
sx q[0];
rz(-2.1053173) q[0];
sx q[0];
rz(-0.86877862) q[0];
rz(2.2390305) q[1];
sx q[1];
rz(-1.2259918) q[1];
sx q[1];
rz(2.5696519) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4298852) q[0];
sx q[0];
rz(-1.0957624) q[0];
sx q[0];
rz(2.320313) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4221572) q[2];
sx q[2];
rz(-0.63221064) q[2];
sx q[2];
rz(-2.3001463) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8780788) q[1];
sx q[1];
rz(-2.0058419) q[1];
sx q[1];
rz(0.33532354) q[1];
rz(-1.9492416) q[3];
sx q[3];
rz(-1.3406383) q[3];
sx q[3];
rz(-0.6164906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14704412) q[2];
sx q[2];
rz(-2.2037555) q[2];
sx q[2];
rz(-0.63329548) q[2];
rz(-2.5607732) q[3];
sx q[3];
rz(-1.7827026) q[3];
sx q[3];
rz(0.66494989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13034114) q[0];
sx q[0];
rz(-0.3891775) q[0];
sx q[0];
rz(-1.8036386) q[0];
rz(-1.411571) q[1];
sx q[1];
rz(-2.7346225) q[1];
sx q[1];
rz(2.3445047) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0397698) q[0];
sx q[0];
rz(-1.2972141) q[0];
sx q[0];
rz(-0.16353971) q[0];
rz(-pi) q[1];
rz(2.1363968) q[2];
sx q[2];
rz(-1.6365286) q[2];
sx q[2];
rz(2.1594723) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5108899) q[1];
sx q[1];
rz(-1.4107009) q[1];
sx q[1];
rz(1.2994205) q[1];
rz(-pi) q[2];
rz(-0.18081801) q[3];
sx q[3];
rz(-1.0581087) q[3];
sx q[3];
rz(-0.11845438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36593124) q[2];
sx q[2];
rz(-2.4727827) q[2];
sx q[2];
rz(0.68525165) q[2];
rz(-2.8819486) q[3];
sx q[3];
rz(-0.97340596) q[3];
sx q[3];
rz(1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.892136) q[0];
sx q[0];
rz(-0.15661713) q[0];
sx q[0];
rz(2.6020965) q[0];
rz(-2.9501713) q[1];
sx q[1];
rz(-1.4861264) q[1];
sx q[1];
rz(-0.44245455) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434261) q[0];
sx q[0];
rz(-1.5616722) q[0];
sx q[0];
rz(-3.1412499) q[0];
x q[1];
rz(0.94687825) q[2];
sx q[2];
rz(-1.2420474) q[2];
sx q[2];
rz(-0.88916949) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8129706) q[1];
sx q[1];
rz(-0.96862205) q[1];
sx q[1];
rz(-0.41090907) q[1];
rz(-pi) q[2];
rz(-0.70480736) q[3];
sx q[3];
rz(-2.2205345) q[3];
sx q[3];
rz(-2.5832502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29360867) q[2];
sx q[2];
rz(-2.831735) q[2];
sx q[2];
rz(0.7027182) q[2];
rz(0.27788776) q[3];
sx q[3];
rz(-1.0669471) q[3];
sx q[3];
rz(-1.8028629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24669312) q[0];
sx q[0];
rz(-0.86240697) q[0];
sx q[0];
rz(0.53584164) q[0];
rz(2.4193173) q[1];
sx q[1];
rz(-1.7012137) q[1];
sx q[1];
rz(-2.3586418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0412126) q[0];
sx q[0];
rz(-1.4836652) q[0];
sx q[0];
rz(1.210733) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7442619) q[2];
sx q[2];
rz(-1.518992) q[2];
sx q[2];
rz(-1.0578294) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6668876) q[1];
sx q[1];
rz(-1.0586481) q[1];
sx q[1];
rz(0.94774232) q[1];
rz(-1.0231185) q[3];
sx q[3];
rz(-1.9503106) q[3];
sx q[3];
rz(-1.4364157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7213514) q[2];
sx q[2];
rz(-0.46478096) q[2];
sx q[2];
rz(2.6381524) q[2];
rz(-0.33468801) q[3];
sx q[3];
rz(-1.8036333) q[3];
sx q[3];
rz(-0.23505841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.47401416) q[0];
sx q[0];
rz(-2.773371) q[0];
sx q[0];
rz(-1.0598805) q[0];
rz(-0.31479752) q[1];
sx q[1];
rz(-1.3455201) q[1];
sx q[1];
rz(-1.559929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2602291) q[0];
sx q[0];
rz(-0.90103982) q[0];
sx q[0];
rz(1.8561646) q[0];
rz(-pi) q[1];
x q[1];
rz(1.694881) q[2];
sx q[2];
rz(-0.26657924) q[2];
sx q[2];
rz(-0.82908344) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0421903) q[1];
sx q[1];
rz(-1.785197) q[1];
sx q[1];
rz(-0.91941388) q[1];
rz(3.1195333) q[3];
sx q[3];
rz(-1.1556871) q[3];
sx q[3];
rz(-0.65655113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.18610893) q[2];
sx q[2];
rz(-1.9766221) q[2];
sx q[2];
rz(1.0355518) q[2];
rz(1.995685) q[3];
sx q[3];
rz(-1.112554) q[3];
sx q[3];
rz(-0.15315332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.938852) q[0];
sx q[0];
rz(-0.11883141) q[0];
sx q[0];
rz(0.76549292) q[0];
rz(-1.0368404) q[1];
sx q[1];
rz(-1.8427269) q[1];
sx q[1];
rz(3.0293005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85690391) q[0];
sx q[0];
rz(-0.56988003) q[0];
sx q[0];
rz(0.63473746) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1767581) q[2];
sx q[2];
rz(-0.8886742) q[2];
sx q[2];
rz(1.8471931) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0621418) q[1];
sx q[1];
rz(-1.3704668) q[1];
sx q[1];
rz(-0.82605861) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0946135) q[3];
sx q[3];
rz(-0.78441294) q[3];
sx q[3];
rz(1.3734773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74467337) q[2];
sx q[2];
rz(-1.8645218) q[2];
sx q[2];
rz(0.12167682) q[2];
rz(-0.35950288) q[3];
sx q[3];
rz(-2.4183141) q[3];
sx q[3];
rz(3.1239037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59035463) q[0];
sx q[0];
rz(-1.6726765) q[0];
sx q[0];
rz(-1.8400675) q[0];
rz(-1.3941258) q[1];
sx q[1];
rz(-1.9448517) q[1];
sx q[1];
rz(1.3762884) q[1];
rz(-0.6223901) q[2];
sx q[2];
rz(-1.9303106) q[2];
sx q[2];
rz(2.8091693) q[2];
rz(-1.2944503) q[3];
sx q[3];
rz(-0.45062267) q[3];
sx q[3];
rz(-2.9723321) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
