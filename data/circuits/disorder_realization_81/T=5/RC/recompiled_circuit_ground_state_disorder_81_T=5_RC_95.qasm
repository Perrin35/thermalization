OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0185952) q[0];
sx q[0];
rz(-1.3774435) q[0];
sx q[0];
rz(-0.42940816) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(-1.4338926) q[1];
sx q[1];
rz(-1.5996999) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7943952) q[0];
sx q[0];
rz(-1.7849677) q[0];
sx q[0];
rz(1.8632351) q[0];
x q[1];
rz(0.1437896) q[2];
sx q[2];
rz(-1.5359288) q[2];
sx q[2];
rz(-0.26511017) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1356537) q[1];
sx q[1];
rz(-2.7126813) q[1];
sx q[1];
rz(-3.0687544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2774994) q[3];
sx q[3];
rz(-1.2240922) q[3];
sx q[3];
rz(-1.5145921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8371007) q[2];
sx q[2];
rz(-1.203275) q[2];
sx q[2];
rz(-2.4486747) q[2];
rz(2.9414226) q[3];
sx q[3];
rz(-2.9514511) q[3];
sx q[3];
rz(2.0076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8189341) q[0];
sx q[0];
rz(-0.19911961) q[0];
sx q[0];
rz(1.06485) q[0];
rz(0.62243593) q[1];
sx q[1];
rz(-1.7935926) q[1];
sx q[1];
rz(0.49076864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9483295) q[0];
sx q[0];
rz(-1.5800467) q[0];
sx q[0];
rz(-3.1178283) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6726637) q[2];
sx q[2];
rz(-0.46655077) q[2];
sx q[2];
rz(0.55794898) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.307617) q[1];
sx q[1];
rz(-0.81980356) q[1];
sx q[1];
rz(1.4403254) q[1];
rz(-pi) q[2];
x q[2];
rz(1.845188) q[3];
sx q[3];
rz(-0.17582045) q[3];
sx q[3];
rz(-0.71064132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48851594) q[2];
sx q[2];
rz(-1.4986897) q[2];
sx q[2];
rz(-0.24031362) q[2];
rz(-0.49992418) q[3];
sx q[3];
rz(-2.5559055) q[3];
sx q[3];
rz(1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2568473) q[0];
sx q[0];
rz(-2.186543) q[0];
sx q[0];
rz(0.98168674) q[0];
rz(0.49199545) q[1];
sx q[1];
rz(-2.056608) q[1];
sx q[1];
rz(1.6384151) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5971736) q[0];
sx q[0];
rz(-1.0962631) q[0];
sx q[0];
rz(-2.3719792) q[0];
x q[1];
rz(-2.7415034) q[2];
sx q[2];
rz(-1.5201668) q[2];
sx q[2];
rz(-2.7226677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3635849) q[1];
sx q[1];
rz(-0.52138803) q[1];
sx q[1];
rz(-0.85659493) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7812626) q[3];
sx q[3];
rz(-0.099848824) q[3];
sx q[3];
rz(2.9304402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4846399) q[2];
sx q[2];
rz(-2.6159365) q[2];
sx q[2];
rz(3.0204115) q[2];
rz(-1.0080522) q[3];
sx q[3];
rz(-1.1153778) q[3];
sx q[3];
rz(-2.0602267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0722395) q[0];
sx q[0];
rz(-3.1407052) q[0];
sx q[0];
rz(-2.6278611) q[0];
rz(2.1836102) q[1];
sx q[1];
rz(-1.999141) q[1];
sx q[1];
rz(-1.6987919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0009036) q[0];
sx q[0];
rz(-1.1620518) q[0];
sx q[0];
rz(-1.7501786) q[0];
rz(-pi) q[1];
rz(0.53057495) q[2];
sx q[2];
rz(-1.3791891) q[2];
sx q[2];
rz(2.8838188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.629238) q[1];
sx q[1];
rz(-1.4886453) q[1];
sx q[1];
rz(-0.95463971) q[1];
rz(-pi) q[2];
rz(1.4738655) q[3];
sx q[3];
rz(-0.91141322) q[3];
sx q[3];
rz(1.3644548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55502597) q[2];
sx q[2];
rz(-2.1777966) q[2];
sx q[2];
rz(-1.2653992) q[2];
rz(1.966656) q[3];
sx q[3];
rz(-2.411071) q[3];
sx q[3];
rz(-1.9780212) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2555399) q[0];
sx q[0];
rz(-0.20296725) q[0];
sx q[0];
rz(1.8170005) q[0];
rz(-2.1728204) q[1];
sx q[1];
rz(-1.6561534) q[1];
sx q[1];
rz(-2.103215) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7234112) q[0];
sx q[0];
rz(-1.1901642) q[0];
sx q[0];
rz(1.7443954) q[0];
rz(-pi) q[1];
rz(0.71717306) q[2];
sx q[2];
rz(-0.96276564) q[2];
sx q[2];
rz(-3.0423683) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15971359) q[1];
sx q[1];
rz(-2.7687589) q[1];
sx q[1];
rz(0.11319686) q[1];
x q[2];
rz(-2.0821003) q[3];
sx q[3];
rz(-2.1645774) q[3];
sx q[3];
rz(-0.64835129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5021299) q[2];
sx q[2];
rz(-0.30439964) q[2];
sx q[2];
rz(-2.2501865) q[2];
rz(-1.8799479) q[3];
sx q[3];
rz(-1.6311878) q[3];
sx q[3];
rz(0.6663028) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0797743) q[0];
sx q[0];
rz(-2.6462055) q[0];
sx q[0];
rz(2.1451982) q[0];
rz(0.90006104) q[1];
sx q[1];
rz(-2.3521017) q[1];
sx q[1];
rz(-0.17734227) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19781659) q[0];
sx q[0];
rz(-1.4192441) q[0];
sx q[0];
rz(3.1413548) q[0];
rz(-1.5992237) q[2];
sx q[2];
rz(-0.6387944) q[2];
sx q[2];
rz(-0.5109238) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23456043) q[1];
sx q[1];
rz(-1.8563358) q[1];
sx q[1];
rz(-1.0726867) q[1];
rz(-pi) q[2];
rz(1.3979891) q[3];
sx q[3];
rz(-2.6263642) q[3];
sx q[3];
rz(2.990641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29048723) q[2];
sx q[2];
rz(-1.1815716) q[2];
sx q[2];
rz(-2.8372852) q[2];
rz(-2.815222) q[3];
sx q[3];
rz(-0.54309741) q[3];
sx q[3];
rz(0.91199005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07311634) q[0];
sx q[0];
rz(-0.40962064) q[0];
sx q[0];
rz(-0.9325183) q[0];
rz(3.0724691) q[1];
sx q[1];
rz(-0.2176452) q[1];
sx q[1];
rz(2.1014012) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5297663) q[0];
sx q[0];
rz(-1.5965723) q[0];
sx q[0];
rz(-2.3759205) q[0];
rz(1.3142775) q[2];
sx q[2];
rz(-1.4019012) q[2];
sx q[2];
rz(-1.0948563) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2582764) q[1];
sx q[1];
rz(-0.81352106) q[1];
sx q[1];
rz(1.798428) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4916441) q[3];
sx q[3];
rz(-2.3293265) q[3];
sx q[3];
rz(-2.5358729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13504623) q[2];
sx q[2];
rz(-0.73837215) q[2];
sx q[2];
rz(-2.4086003) q[2];
rz(2.0586355) q[3];
sx q[3];
rz(-1.0561918) q[3];
sx q[3];
rz(-2.1347031) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55423823) q[0];
sx q[0];
rz(-0.08389689) q[0];
sx q[0];
rz(-2.6485637) q[0];
rz(2.9250277) q[1];
sx q[1];
rz(-1.1410057) q[1];
sx q[1];
rz(-2.1176178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7465377) q[0];
sx q[0];
rz(-2.1667272) q[0];
sx q[0];
rz(1.5851969) q[0];
rz(-pi) q[1];
rz(-0.31957133) q[2];
sx q[2];
rz(-2.7929578) q[2];
sx q[2];
rz(-0.92952585) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1962016) q[1];
sx q[1];
rz(-1.7111254) q[1];
sx q[1];
rz(2.7811471) q[1];
rz(-pi) q[2];
rz(2.2129503) q[3];
sx q[3];
rz(-1.3285884) q[3];
sx q[3];
rz(0.95706576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0563125) q[2];
sx q[2];
rz(-1.5093426) q[2];
sx q[2];
rz(-2.4658266) q[2];
rz(-2.7285649) q[3];
sx q[3];
rz(-1.4512117) q[3];
sx q[3];
rz(-2.5954424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5125047) q[0];
sx q[0];
rz(-0.65852037) q[0];
sx q[0];
rz(-0.034544695) q[0];
rz(0.054556219) q[1];
sx q[1];
rz(-1.1628393) q[1];
sx q[1];
rz(-2.3202855) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6550044) q[0];
sx q[0];
rz(-1.6704217) q[0];
sx q[0];
rz(2.5667046) q[0];
x q[1];
rz(1.5615145) q[2];
sx q[2];
rz(-1.6595652) q[2];
sx q[2];
rz(2.2616158) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7682338) q[1];
sx q[1];
rz(-1.7204509) q[1];
sx q[1];
rz(1.4537158) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48252941) q[3];
sx q[3];
rz(-0.97597968) q[3];
sx q[3];
rz(-1.9687259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7811232) q[2];
sx q[2];
rz(-2.1253773) q[2];
sx q[2];
rz(-3.0086369) q[2];
rz(-0.41108701) q[3];
sx q[3];
rz(-1.5186331) q[3];
sx q[3];
rz(1.0857922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0947615) q[0];
sx q[0];
rz(-0.60698858) q[0];
sx q[0];
rz(-1.2588311) q[0];
rz(-1.6350485) q[1];
sx q[1];
rz(-2.4848487) q[1];
sx q[1];
rz(0.81319317) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5482322) q[0];
sx q[0];
rz(-1.5357067) q[0];
sx q[0];
rz(3.0635053) q[0];
rz(-pi) q[1];
rz(-2.2426064) q[2];
sx q[2];
rz(-1.4788663) q[2];
sx q[2];
rz(-1.7882333) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.836536) q[1];
sx q[1];
rz(-0.99545762) q[1];
sx q[1];
rz(2.1161377) q[1];
rz(-pi) q[2];
rz(1.7045295) q[3];
sx q[3];
rz(-1.1220891) q[3];
sx q[3];
rz(2.1577842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4930111) q[2];
sx q[2];
rz(-1.1375256) q[2];
sx q[2];
rz(-1.0515593) q[2];
rz(2.1227396) q[3];
sx q[3];
rz(-2.2994883) q[3];
sx q[3];
rz(-2.4837608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1558341) q[0];
sx q[0];
rz(-0.65407615) q[0];
sx q[0];
rz(0.88055897) q[0];
rz(1.3195994) q[1];
sx q[1];
rz(-1.1450014) q[1];
sx q[1];
rz(-3.0897279) q[1];
rz(0.20959494) q[2];
sx q[2];
rz(-1.4996698) q[2];
sx q[2];
rz(-1.3347129) q[2];
rz(-1.6307835) q[3];
sx q[3];
rz(-1.7675478) q[3];
sx q[3];
rz(0.02789733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
