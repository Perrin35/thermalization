OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(-0.72609225) q[0];
sx q[0];
rz(-0.2015764) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(2.6013241) q[1];
sx q[1];
rz(7.2202914) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0721604) q[0];
sx q[0];
rz(-1.1496815) q[0];
sx q[0];
rz(-0.62299563) q[0];
rz(-pi) q[1];
rz(-0.61203981) q[2];
sx q[2];
rz(-1.9844374) q[2];
sx q[2];
rz(-2.1940102) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.01552445) q[1];
sx q[1];
rz(-1.4926732) q[1];
sx q[1];
rz(1.7541581) q[1];
rz(-1.8398383) q[3];
sx q[3];
rz(-1.6645414) q[3];
sx q[3];
rz(-0.71414381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9822838) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(3.0554331) q[2];
rz(2.384095) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(2.0479726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57698292) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(1.3264867) q[0];
rz(1.8857229) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(2.870141) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87646987) q[0];
sx q[0];
rz(-1.0713097) q[0];
sx q[0];
rz(-2.6861517) q[0];
rz(-0.51241264) q[2];
sx q[2];
rz(-2.5645442) q[2];
sx q[2];
rz(-1.7713828) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93310624) q[1];
sx q[1];
rz(-2.3098364) q[1];
sx q[1];
rz(2.0887124) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80941697) q[3];
sx q[3];
rz(-1.3076233) q[3];
sx q[3];
rz(-1.1337048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0945956) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(-0.57717741) q[2];
rz(0.92352891) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975824) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(2.999021) q[0];
rz(1.7890731) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(2.9325063) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7798628) q[0];
sx q[0];
rz(-1.0233101) q[0];
sx q[0];
rz(-0.89625432) q[0];
rz(-pi) q[1];
rz(-2.694391) q[2];
sx q[2];
rz(-2.3502091) q[2];
sx q[2];
rz(0.42391047) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1413893) q[1];
sx q[1];
rz(-0.63265002) q[1];
sx q[1];
rz(-0.81694095) q[1];
rz(-pi) q[2];
rz(2.8385931) q[3];
sx q[3];
rz(-1.3649366) q[3];
sx q[3];
rz(2.1752165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7769988) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(2.7374632) q[2];
rz(1.2858307) q[3];
sx q[3];
rz(-1.1288246) q[3];
sx q[3];
rz(-2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362815) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(-2.1787815) q[0];
rz(-2.6722233) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(-0.00096360047) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92361802) q[0];
sx q[0];
rz(-2.1305363) q[0];
sx q[0];
rz(-0.864242) q[0];
x q[1];
rz(1.5713111) q[2];
sx q[2];
rz(-1.4467903) q[2];
sx q[2];
rz(-2.4027783) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5597792) q[1];
sx q[1];
rz(-1.0043136) q[1];
sx q[1];
rz(-2.986746) q[1];
rz(-1.6603052) q[3];
sx q[3];
rz(-2.138391) q[3];
sx q[3];
rz(-1.2168509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0756388) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(3.0299419) q[2];
rz(0.81104898) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(3.1276935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951293) q[0];
sx q[0];
rz(-0.90068978) q[0];
sx q[0];
rz(-2.7850889) q[0];
rz(0.50645343) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(-2.8809663) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168419) q[0];
sx q[0];
rz(-1.5289613) q[0];
sx q[0];
rz(0.14194685) q[0];
rz(2.3512958) q[2];
sx q[2];
rz(-2.7919263) q[2];
sx q[2];
rz(-0.5860354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.063559859) q[1];
sx q[1];
rz(-1.5884001) q[1];
sx q[1];
rz(-2.5672008) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14560933) q[3];
sx q[3];
rz(-1.7123316) q[3];
sx q[3];
rz(1.4841929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9115209) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(2.9122706) q[2];
rz(-2.5991332) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(-1.0361766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(-1.649958) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.7274436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8261674) q[0];
sx q[0];
rz(-0.27914473) q[0];
sx q[0];
rz(0.17255737) q[0];
rz(-pi) q[1];
rz(2.6944689) q[2];
sx q[2];
rz(-0.56338718) q[2];
sx q[2];
rz(-0.901957) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9767178) q[1];
sx q[1];
rz(-2.4396982) q[1];
sx q[1];
rz(-2.328842) q[1];
x q[2];
rz(1.9023444) q[3];
sx q[3];
rz(-0.22833951) q[3];
sx q[3];
rz(2.6339298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(-3.1138528) q[2];
rz(-2.6489143) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(-1.3940575) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7572927) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(3.0199155) q[0];
rz(-1.9901468) q[1];
sx q[1];
rz(-0.45184389) q[1];
sx q[1];
rz(-0.40245232) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6399022) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(-0.34696607) q[0];
rz(2.9629838) q[2];
sx q[2];
rz(-3.0948234) q[2];
sx q[2];
rz(0.59431078) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.788159) q[1];
sx q[1];
rz(-2.806059) q[1];
sx q[1];
rz(2.5787756) q[1];
rz(-pi) q[2];
rz(-1.6472858) q[3];
sx q[3];
rz(-0.34919958) q[3];
sx q[3];
rz(2.4793712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1272614) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(-2.2793615) q[2];
rz(-0.47752738) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(-2.2235218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35571337) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(-2.4086337) q[0];
rz(3.0015302) q[1];
sx q[1];
rz(-0.99761325) q[1];
sx q[1];
rz(-1.0345116) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1326133) q[0];
sx q[0];
rz(-1.7404557) q[0];
sx q[0];
rz(-1.3404113) q[0];
rz(-pi) q[1];
rz(0.10266281) q[2];
sx q[2];
rz(-1.6409988) q[2];
sx q[2];
rz(1.3753124) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.58971436) q[1];
sx q[1];
rz(-2.1363746) q[1];
sx q[1];
rz(-3.0977071) q[1];
rz(-pi) q[2];
rz(-2.7233852) q[3];
sx q[3];
rz(-1.6697262) q[3];
sx q[3];
rz(1.4560771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2404279) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(-0.77587664) q[2];
rz(-2.4173229) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(-1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6999321) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(3.124776) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-2.3628078) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9433702) q[0];
sx q[0];
rz(-1.3226042) q[0];
sx q[0];
rz(-2.7944837) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8674556) q[2];
sx q[2];
rz(-2.0275653) q[2];
sx q[2];
rz(2.0704839) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3719337) q[1];
sx q[1];
rz(-1.9366591) q[1];
sx q[1];
rz(-0.2497754) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0770244) q[3];
sx q[3];
rz(-1.350292) q[3];
sx q[3];
rz(1.8622423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-2.4251535) q[2];
rz(-1.4572432) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(-1.8523857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(-2.9901436) q[0];
rz(-0.17627136) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(-0.72296468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3176214) q[0];
sx q[0];
rz(-0.184632) q[0];
sx q[0];
rz(2.1872107) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9724571) q[2];
sx q[2];
rz(-1.198487) q[2];
sx q[2];
rz(1.362209) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42230095) q[1];
sx q[1];
rz(-2.0740168) q[1];
sx q[1];
rz(3.0784688) q[1];
x q[2];
rz(1.1687752) q[3];
sx q[3];
rz(-2.3178181) q[3];
sx q[3];
rz(-1.0812425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67939776) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(-0.52552137) q[2];
rz(-0.28371352) q[3];
sx q[3];
rz(-2.0339537) q[3];
sx q[3];
rz(0.39653683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239607) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(2.0422968) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(0.47808403) q[2];
sx q[2];
rz(-2.1838084) q[2];
sx q[2];
rz(2.9391391) q[2];
rz(-2.5476534) q[3];
sx q[3];
rz(-2.754302) q[3];
sx q[3];
rz(2.1828628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
