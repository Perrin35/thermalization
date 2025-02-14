OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.60431689) q[0];
sx q[0];
rz(3.3831626) q[0];
sx q[0];
rz(9.0917505) q[0];
rz(-1.1355407) q[1];
sx q[1];
rz(3.9685213) q[1];
sx q[1];
rz(10.068738) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0892093) q[0];
sx q[0];
rz(-0.96714562) q[0];
sx q[0];
rz(0.33120819) q[0];
rz(-1.7276554) q[2];
sx q[2];
rz(-1.2741977) q[2];
sx q[2];
rz(0.083904249) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2702858) q[1];
sx q[1];
rz(-1.8075004) q[1];
sx q[1];
rz(-0.38856296) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9831564) q[3];
sx q[3];
rz(-1.8915911) q[3];
sx q[3];
rz(0.12925805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48646271) q[2];
sx q[2];
rz(-0.4643521) q[2];
sx q[2];
rz(0.74074024) q[2];
rz(-1.5898534) q[3];
sx q[3];
rz(-2.430075) q[3];
sx q[3];
rz(-0.5927425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084259) q[0];
sx q[0];
rz(-1.1335224) q[0];
sx q[0];
rz(0.11696996) q[0];
rz(2.6372657) q[1];
sx q[1];
rz(-1.802899) q[1];
sx q[1];
rz(0.28409827) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5829651) q[0];
sx q[0];
rz(-2.5694048) q[0];
sx q[0];
rz(1.7250604) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5710377) q[2];
sx q[2];
rz(-0.88316702) q[2];
sx q[2];
rz(0.90435435) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85728474) q[1];
sx q[1];
rz(-1.4487639) q[1];
sx q[1];
rz(-0.034842592) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3044182) q[3];
sx q[3];
rz(-2.7751013) q[3];
sx q[3];
rz(-0.68658406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.32610193) q[2];
sx q[2];
rz(-0.88560605) q[2];
sx q[2];
rz(2.3808114) q[2];
rz(0.86959362) q[3];
sx q[3];
rz(-1.2660916) q[3];
sx q[3];
rz(-0.45309711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78075439) q[0];
sx q[0];
rz(-2.0212845) q[0];
sx q[0];
rz(0.25217062) q[0];
rz(-1.699327) q[1];
sx q[1];
rz(-2.1863054) q[1];
sx q[1];
rz(-0.35626492) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8667135) q[0];
sx q[0];
rz(-1.5646114) q[0];
sx q[0];
rz(-1.5112108) q[0];
rz(-2.5926931) q[2];
sx q[2];
rz(-1.2147012) q[2];
sx q[2];
rz(-0.20129542) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60881847) q[1];
sx q[1];
rz(-1.7131117) q[1];
sx q[1];
rz(-2.2102093) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3144794) q[3];
sx q[3];
rz(-1.5689092) q[3];
sx q[3];
rz(2.6607115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0226655) q[2];
sx q[2];
rz(-1.3902731) q[2];
sx q[2];
rz(-1.5125795) q[2];
rz(0.0096983612) q[3];
sx q[3];
rz(-1.9111948) q[3];
sx q[3];
rz(-0.62140083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464722) q[0];
sx q[0];
rz(-0.54238129) q[0];
sx q[0];
rz(-2.8908308) q[0];
rz(-1.3465025) q[1];
sx q[1];
rz(-0.52918068) q[1];
sx q[1];
rz(1.6983324) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8478717) q[0];
sx q[0];
rz(-1.3547055) q[0];
sx q[0];
rz(-2.3225075) q[0];
rz(-pi) q[1];
rz(1.0403149) q[2];
sx q[2];
rz(-1.217739) q[2];
sx q[2];
rz(-0.7952035) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5233744) q[1];
sx q[1];
rz(-2.187145) q[1];
sx q[1];
rz(-0.90202721) q[1];
rz(2.3425927) q[3];
sx q[3];
rz(-2.7741787) q[3];
sx q[3];
rz(-1.1297117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.5556339) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(0.24331681) q[2];
rz(2.5569052) q[3];
sx q[3];
rz(-0.57716113) q[3];
sx q[3];
rz(3.1134591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.993416) q[0];
sx q[0];
rz(-2.7758444) q[0];
sx q[0];
rz(-2.962033) q[0];
rz(-1.9163632) q[1];
sx q[1];
rz(-1.4045249) q[1];
sx q[1];
rz(-0.45825759) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62406926) q[0];
sx q[0];
rz(-1.6485414) q[0];
sx q[0];
rz(0.29645424) q[0];
x q[1];
rz(-1.304428) q[2];
sx q[2];
rz(-2.2812216) q[2];
sx q[2];
rz(-1.2183587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1128487) q[1];
sx q[1];
rz(-2.9874871) q[1];
sx q[1];
rz(-1.0636368) q[1];
rz(-pi) q[2];
rz(-2.7814976) q[3];
sx q[3];
rz(-0.90237877) q[3];
sx q[3];
rz(-0.51262142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3852343) q[2];
sx q[2];
rz(-0.42432722) q[2];
sx q[2];
rz(2.2198086) q[2];
rz(1.8900169) q[3];
sx q[3];
rz(-1.7332964) q[3];
sx q[3];
rz(3.0799227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09403041) q[0];
sx q[0];
rz(-2.229409) q[0];
sx q[0];
rz(-0.14866522) q[0];
rz(-0.31691638) q[1];
sx q[1];
rz(-2.6628559) q[1];
sx q[1];
rz(-2.5247578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7148833) q[0];
sx q[0];
rz(-1.5093818) q[0];
sx q[0];
rz(1.4922569) q[0];
x q[1];
rz(-0.89646879) q[2];
sx q[2];
rz(-2.5988262) q[2];
sx q[2];
rz(1.8405869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1779536) q[1];
sx q[1];
rz(-0.63040367) q[1];
sx q[1];
rz(-2.0042581) q[1];
rz(-pi) q[2];
rz(2.8657593) q[3];
sx q[3];
rz(-1.0636998) q[3];
sx q[3];
rz(2.3292993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2769015) q[2];
sx q[2];
rz(-1.7128877) q[2];
sx q[2];
rz(-0.0083262715) q[2];
rz(0.62638038) q[3];
sx q[3];
rz(-2.4255987) q[3];
sx q[3];
rz(-3.1410419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83547) q[0];
sx q[0];
rz(-1.1518814) q[0];
sx q[0];
rz(2.6595111) q[0];
rz(0.029190633) q[1];
sx q[1];
rz(-1.8455285) q[1];
sx q[1];
rz(2.3775502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0917863) q[0];
sx q[0];
rz(-0.37702628) q[0];
sx q[0];
rz(-1.754712) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5381728) q[2];
sx q[2];
rz(-0.66002405) q[2];
sx q[2];
rz(0.78885022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.205207) q[1];
sx q[1];
rz(-2.1644785) q[1];
sx q[1];
rz(2.9365345) q[1];
rz(-pi) q[2];
rz(1.589631) q[3];
sx q[3];
rz(-0.40098396) q[3];
sx q[3];
rz(-2.8249521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46334106) q[2];
sx q[2];
rz(-1.8222787) q[2];
sx q[2];
rz(0.48416644) q[2];
rz(2.4593027) q[3];
sx q[3];
rz(-0.65410084) q[3];
sx q[3];
rz(0.13474034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084455647) q[0];
sx q[0];
rz(-0.77780044) q[0];
sx q[0];
rz(-1.8498259) q[0];
rz(-2.6765587) q[1];
sx q[1];
rz(-0.51883042) q[1];
sx q[1];
rz(-2.8344287) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7028593) q[0];
sx q[0];
rz(-2.4552058) q[0];
sx q[0];
rz(0.86875654) q[0];
rz(-0.50339209) q[2];
sx q[2];
rz(-2.1648228) q[2];
sx q[2];
rz(2.4979748) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7502082) q[1];
sx q[1];
rz(-0.5942052) q[1];
sx q[1];
rz(-1.8529525) q[1];
x q[2];
rz(-1.5504595) q[3];
sx q[3];
rz(-2.1957948) q[3];
sx q[3];
rz(2.60256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9610567) q[2];
sx q[2];
rz(-2.7013216) q[2];
sx q[2];
rz(-0.72009909) q[2];
rz(-0.85047203) q[3];
sx q[3];
rz(-1.1668147) q[3];
sx q[3];
rz(1.4115964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20486031) q[0];
sx q[0];
rz(-2.9731049) q[0];
sx q[0];
rz(-2.4334461) q[0];
rz(-1.5180961) q[1];
sx q[1];
rz(-2.203439) q[1];
sx q[1];
rz(0.35619563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1461648) q[0];
sx q[0];
rz(-2.0588027) q[0];
sx q[0];
rz(1.0900709) q[0];
x q[1];
rz(0.61788606) q[2];
sx q[2];
rz(-1.5466585) q[2];
sx q[2];
rz(1.2302421) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6345919) q[1];
sx q[1];
rz(-2.2193877) q[1];
sx q[1];
rz(-1.533094) q[1];
x q[2];
rz(2.0354) q[3];
sx q[3];
rz(-1.937791) q[3];
sx q[3];
rz(1.3657686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0922962) q[2];
sx q[2];
rz(-2.3342817) q[2];
sx q[2];
rz(2.2558007) q[2];
rz(0.93585912) q[3];
sx q[3];
rz(-1.1805308) q[3];
sx q[3];
rz(-2.9203171) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7037999) q[0];
sx q[0];
rz(-3.0942823) q[0];
sx q[0];
rz(1.4495151) q[0];
rz(-1.0225147) q[1];
sx q[1];
rz(-0.48373628) q[1];
sx q[1];
rz(1.449301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3614907) q[0];
sx q[0];
rz(-1.6204483) q[0];
sx q[0];
rz(-1.1583369) q[0];
x q[1];
rz(1.9841393) q[2];
sx q[2];
rz(-1.1755694) q[2];
sx q[2];
rz(-0.64973598) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8359747) q[1];
sx q[1];
rz(-0.59594369) q[1];
sx q[1];
rz(1.9553595) q[1];
rz(-2.9935097) q[3];
sx q[3];
rz(-1.3080773) q[3];
sx q[3];
rz(2.1488291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.25541043) q[2];
sx q[2];
rz(-1.9516727) q[2];
sx q[2];
rz(-0.6748684) q[2];
rz(2.1700962) q[3];
sx q[3];
rz(-0.70236218) q[3];
sx q[3];
rz(-1.6500047) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3424727) q[0];
sx q[0];
rz(-1.5958888) q[0];
sx q[0];
rz(2.2842443) q[0];
rz(-2.9091861) q[1];
sx q[1];
rz(-2.1288165) q[1];
sx q[1];
rz(-1.7850599) q[1];
rz(2.7276298) q[2];
sx q[2];
rz(-2.064075) q[2];
sx q[2];
rz(1.0309564) q[2];
rz(3.0633756) q[3];
sx q[3];
rz(-0.84368869) q[3];
sx q[3];
rz(0.61144184) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
