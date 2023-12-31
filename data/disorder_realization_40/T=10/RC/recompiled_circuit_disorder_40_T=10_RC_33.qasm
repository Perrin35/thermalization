OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(-3.0298046) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(-2.9878374) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.114404) q[0];
sx q[0];
rz(-1.1429042) q[0];
sx q[0];
rz(1.2501636) q[0];
rz(-3.0646677) q[2];
sx q[2];
rz(-1.3300606) q[2];
sx q[2];
rz(-1.392729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1195214) q[1];
sx q[1];
rz(-1.6011366) q[1];
sx q[1];
rz(1.2519217) q[1];
x q[2];
rz(2.8044756) q[3];
sx q[3];
rz(-1.1532591) q[3];
sx q[3];
rz(1.2085714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6979606) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(1.704818) q[2];
rz(-0.73389655) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88965082) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(0.92457986) q[0];
rz(2.1444767) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(1.3234214) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2572718) q[0];
sx q[0];
rz(-1.1182251) q[0];
sx q[0];
rz(1.0898468) q[0];
rz(2.7919865) q[2];
sx q[2];
rz(-1.6125624) q[2];
sx q[2];
rz(-1.2431527) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4393649) q[1];
sx q[1];
rz(-0.80184466) q[1];
sx q[1];
rz(-0.37978362) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38815659) q[3];
sx q[3];
rz(-2.1106488) q[3];
sx q[3];
rz(2.848958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9374342) q[2];
sx q[2];
rz(-1.5519451) q[2];
sx q[2];
rz(-0.75817529) q[2];
rz(0.6289064) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(1.1531856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.73308289) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(-2.4531903) q[0];
rz(0.06772659) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(-0.53007954) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353697) q[0];
sx q[0];
rz(-0.3100937) q[0];
sx q[0];
rz(-1.9276516) q[0];
rz(2.7128503) q[2];
sx q[2];
rz(-2.0566166) q[2];
sx q[2];
rz(0.57810099) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0878151) q[1];
sx q[1];
rz(-0.74790819) q[1];
sx q[1];
rz(2.0762216) q[1];
rz(1.8102874) q[3];
sx q[3];
rz(-1.4756087) q[3];
sx q[3];
rz(0.59323192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3337341) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(-0.90144908) q[2];
rz(0.83550134) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9505342) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(0.78432551) q[0];
rz(-3.0803608) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(-0.13664666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2480859) q[0];
sx q[0];
rz(-2.6419905) q[0];
sx q[0];
rz(1.8462371) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6776671) q[2];
sx q[2];
rz(-1.8564965) q[2];
sx q[2];
rz(2.6505016) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54388753) q[1];
sx q[1];
rz(-1.502431) q[1];
sx q[1];
rz(-0.023936546) q[1];
rz(0.09282077) q[3];
sx q[3];
rz(-0.99675677) q[3];
sx q[3];
rz(-2.1294347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0121997) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(-2.5811035) q[2];
rz(-3.1292606) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(2.0509317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43301582) q[0];
sx q[0];
rz(-0.60537678) q[0];
sx q[0];
rz(-2.3204455) q[0];
rz(-0.87617809) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(1.4076153) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9579904) q[0];
sx q[0];
rz(-1.9200268) q[0];
sx q[0];
rz(-0.41718418) q[0];
x q[1];
rz(-0.9887092) q[2];
sx q[2];
rz(-1.8277797) q[2];
sx q[2];
rz(0.35494057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.69406063) q[1];
sx q[1];
rz(-0.58528712) q[1];
sx q[1];
rz(-1.8156169) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1122583) q[3];
sx q[3];
rz(-2.7280305) q[3];
sx q[3];
rz(-0.13973164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.451482) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(-3.0991128) q[2];
rz(0.63043198) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(-0.49155864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38934389) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(2.65843) q[0];
rz(2.0893611) q[1];
sx q[1];
rz(-1.1455043) q[1];
sx q[1];
rz(2.5767456) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0739581) q[0];
sx q[0];
rz(-1.7059776) q[0];
sx q[0];
rz(1.0413175) q[0];
x q[1];
rz(-0.075750307) q[2];
sx q[2];
rz(-1.1510013) q[2];
sx q[2];
rz(-0.23725739) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1419303) q[1];
sx q[1];
rz(-1.6301486) q[1];
sx q[1];
rz(1.1430986) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23128831) q[3];
sx q[3];
rz(-1.7219208) q[3];
sx q[3];
rz(0.95201991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5715282) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(1.338039) q[2];
rz(1.3048874) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(0.6750955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682482) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(0.54779732) q[0];
rz(0.785218) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(2.8731667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7086605) q[0];
sx q[0];
rz(-1.4888568) q[0];
sx q[0];
rz(-0.24631484) q[0];
rz(-2.3586876) q[2];
sx q[2];
rz(-2.1856538) q[2];
sx q[2];
rz(2.8564786) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.06936) q[1];
sx q[1];
rz(-1.4727122) q[1];
sx q[1];
rz(-2.1588438) q[1];
x q[2];
rz(-1.0309585) q[3];
sx q[3];
rz(-1.8317878) q[3];
sx q[3];
rz(-2.6574082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61775529) q[2];
sx q[2];
rz(-0.80344168) q[2];
sx q[2];
rz(2.3274373) q[2];
rz(-0.37627775) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(3.0686839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58586621) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(-1.0193753) q[0];
rz(-2.2881919) q[1];
sx q[1];
rz(-1.9995721) q[1];
sx q[1];
rz(-2.6928435) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73483) q[0];
sx q[0];
rz(-1.3376298) q[0];
sx q[0];
rz(2.0582817) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2737695) q[2];
sx q[2];
rz(-0.38735577) q[2];
sx q[2];
rz(0.58434904) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6905578) q[1];
sx q[1];
rz(-1.0447377) q[1];
sx q[1];
rz(-1.8655538) q[1];
rz(-2.2711146) q[3];
sx q[3];
rz(-0.92910367) q[3];
sx q[3];
rz(2.2023647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.86924187) q[2];
sx q[2];
rz(-1.7607471) q[2];
sx q[2];
rz(0.31420079) q[2];
rz(-2.3172486) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(-0.79469386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0712414) q[0];
sx q[0];
rz(-3.0817139) q[0];
sx q[0];
rz(1.2605793) q[0];
rz(-2.4977327) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(0.018928122) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2429598) q[0];
sx q[0];
rz(-1.8343933) q[0];
sx q[0];
rz(1.5522869) q[0];
x q[1];
rz(0.55982121) q[2];
sx q[2];
rz(-0.18351843) q[2];
sx q[2];
rz(-0.57932094) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70871204) q[1];
sx q[1];
rz(-1.5203262) q[1];
sx q[1];
rz(-0.18458472) q[1];
rz(-2.4631259) q[3];
sx q[3];
rz(-2.7760091) q[3];
sx q[3];
rz(-1.4640704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54946047) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(2.4712759) q[2];
rz(-0.51236764) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5877514) q[0];
sx q[0];
rz(-1.2441664) q[0];
sx q[0];
rz(-0.49945369) q[0];
rz(-1.5669426) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(2.0589028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4944045) q[0];
sx q[0];
rz(-1.6196961) q[0];
sx q[0];
rz(0.58736579) q[0];
rz(1.535365) q[2];
sx q[2];
rz(-0.75181585) q[2];
sx q[2];
rz(-2.5978136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4808828) q[1];
sx q[1];
rz(-2.1143267) q[1];
sx q[1];
rz(-1.5488312) q[1];
rz(2.4389078) q[3];
sx q[3];
rz(-2.6583238) q[3];
sx q[3];
rz(-2.0600704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3988951) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(2.6297074) q[2];
rz(-0.39294696) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.186541) q[0];
sx q[0];
rz(-1.825009) q[0];
sx q[0];
rz(0.64074989) q[0];
rz(-0.74116771) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(-2.7411186) q[2];
sx q[2];
rz(-0.99292143) q[2];
sx q[2];
rz(-1.6510486) q[2];
rz(-0.18556553) q[3];
sx q[3];
rz(-2.3864828) q[3];
sx q[3];
rz(-1.3226487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
