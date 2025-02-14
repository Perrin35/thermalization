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
rz(-1.8259814) q[0];
sx q[0];
rz(-2.313518) q[0];
sx q[0];
rz(0.84870422) q[0];
rz(1.4915713) q[1];
sx q[1];
rz(2.4529011) q[1];
sx q[1];
rz(12.139763) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4147581) q[0];
sx q[0];
rz(-1.6462741) q[0];
sx q[0];
rz(-3.0407136) q[0];
x q[1];
rz(2.3057904) q[2];
sx q[2];
rz(-2.9907132) q[2];
sx q[2];
rz(2.7250233) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9200807) q[1];
sx q[1];
rz(-2.6903841) q[1];
sx q[1];
rz(0.16772224) q[1];
rz(-0.60517444) q[3];
sx q[3];
rz(-0.89720336) q[3];
sx q[3];
rz(1.0104346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2171057) q[2];
sx q[2];
rz(-1.3873528) q[2];
sx q[2];
rz(3.1021049) q[2];
rz(-1.0314137) q[3];
sx q[3];
rz(-1.02905) q[3];
sx q[3];
rz(1.903681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51204387) q[0];
sx q[0];
rz(-1.5271674) q[0];
sx q[0];
rz(-1.398983) q[0];
rz(-1.6276739) q[1];
sx q[1];
rz(-2.2986423) q[1];
sx q[1];
rz(1.7248076) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2279204) q[0];
sx q[0];
rz(-0.79288542) q[0];
sx q[0];
rz(-0.31358729) q[0];
x q[1];
rz(-0.26160364) q[2];
sx q[2];
rz(-0.79643476) q[2];
sx q[2];
rz(-1.4858044) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4128215) q[1];
sx q[1];
rz(-1.9660453) q[1];
sx q[1];
rz(-1.7186307) q[1];
rz(-pi) q[2];
rz(2.7339277) q[3];
sx q[3];
rz(-1.614794) q[3];
sx q[3];
rz(-0.10213479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17025718) q[2];
sx q[2];
rz(-1.1032666) q[2];
sx q[2];
rz(-2.4309168) q[2];
rz(-3.1273048) q[3];
sx q[3];
rz(-2.894214) q[3];
sx q[3];
rz(-1.8054731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22064848) q[0];
sx q[0];
rz(-2.1823688) q[0];
sx q[0];
rz(-0.4162108) q[0];
rz(-1.8572218) q[1];
sx q[1];
rz(-1.6651848) q[1];
sx q[1];
rz(-2.823337) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6775433) q[0];
sx q[0];
rz(-2.661672) q[0];
sx q[0];
rz(0.93911688) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2936965) q[2];
sx q[2];
rz(-2.090914) q[2];
sx q[2];
rz(-2.0682356) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3102195) q[1];
sx q[1];
rz(-1.4351658) q[1];
sx q[1];
rz(0.65717213) q[1];
rz(-pi) q[2];
rz(1.2622358) q[3];
sx q[3];
rz(-1.8073544) q[3];
sx q[3];
rz(-2.6738338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2661813) q[2];
sx q[2];
rz(-2.3616932) q[2];
sx q[2];
rz(1.2904588) q[2];
rz(-3.087888) q[3];
sx q[3];
rz(-1.3196245) q[3];
sx q[3];
rz(1.3857589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-1.5586435) q[0];
sx q[0];
rz(-2.4531893) q[0];
sx q[0];
rz(2.131856) q[0];
rz(0.91521493) q[1];
sx q[1];
rz(-1.5292294) q[1];
sx q[1];
rz(2.7161652) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0260507) q[0];
sx q[0];
rz(-2.4049691) q[0];
sx q[0];
rz(1.8202341) q[0];
x q[1];
rz(1.6507829) q[2];
sx q[2];
rz(-0.74466193) q[2];
sx q[2];
rz(-2.3415945) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3784619) q[1];
sx q[1];
rz(-1.9768856) q[1];
sx q[1];
rz(-2.2927845) q[1];
rz(-1.5365547) q[3];
sx q[3];
rz(-0.16290191) q[3];
sx q[3];
rz(0.77252856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2414744) q[2];
sx q[2];
rz(-1.590531) q[2];
sx q[2];
rz(-3.0305064) q[2];
rz(-0.19242081) q[3];
sx q[3];
rz(-0.40168732) q[3];
sx q[3];
rz(0.69453159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38839328) q[0];
sx q[0];
rz(-2.41112) q[0];
sx q[0];
rz(2.2764192) q[0];
rz(-2.6490037) q[1];
sx q[1];
rz(-2.045423) q[1];
sx q[1];
rz(1.7561087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3139112) q[0];
sx q[0];
rz(-1.7075141) q[0];
sx q[0];
rz(0.78929957) q[0];
rz(-2.0086914) q[2];
sx q[2];
rz(-0.81134331) q[2];
sx q[2];
rz(1.0824301) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.046987586) q[1];
sx q[1];
rz(-1.054227) q[1];
sx q[1];
rz(2.0265685) q[1];
x q[2];
rz(-2.288743) q[3];
sx q[3];
rz(-2.2519023) q[3];
sx q[3];
rz(-1.8913325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2546786) q[2];
sx q[2];
rz(-0.46923894) q[2];
sx q[2];
rz(0.706642) q[2];
rz(-0.32736579) q[3];
sx q[3];
rz(-1.2923765) q[3];
sx q[3];
rz(-2.4842026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3051598) q[0];
sx q[0];
rz(-1.7921472) q[0];
sx q[0];
rz(0.35655546) q[0];
rz(0.56218475) q[1];
sx q[1];
rz(-1.1327344) q[1];
sx q[1];
rz(-2.1551989) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89740035) q[0];
sx q[0];
rz(-2.2119105) q[0];
sx q[0];
rz(-2.2050307) q[0];
x q[1];
rz(1.2012901) q[2];
sx q[2];
rz(-1.8942648) q[2];
sx q[2];
rz(-0.37959129) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2960259) q[1];
sx q[1];
rz(-2.3289776) q[1];
sx q[1];
rz(-1.6065099) q[1];
x q[2];
rz(-2.9602157) q[3];
sx q[3];
rz(-1.3311609) q[3];
sx q[3];
rz(-2.5549623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76546136) q[2];
sx q[2];
rz(-2.4262846) q[2];
sx q[2];
rz(-2.4410655) q[2];
rz(0.96945196) q[3];
sx q[3];
rz(-2.1078096) q[3];
sx q[3];
rz(-2.2112924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9483865) q[0];
sx q[0];
rz(-0.32873118) q[0];
sx q[0];
rz(-2.9902003) q[0];
rz(1.5104177) q[1];
sx q[1];
rz(-2.0163586) q[1];
sx q[1];
rz(-1.4600533) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9064925) q[0];
sx q[0];
rz(-1.0851881) q[0];
sx q[0];
rz(-0.67295815) q[0];
x q[1];
rz(2.1489819) q[2];
sx q[2];
rz(-0.60727233) q[2];
sx q[2];
rz(-0.36888514) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5826539) q[1];
sx q[1];
rz(-2.451921) q[1];
sx q[1];
rz(0.93780545) q[1];
x q[2];
rz(-0.93029706) q[3];
sx q[3];
rz(-0.2411763) q[3];
sx q[3];
rz(-1.1080139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.90525642) q[2];
sx q[2];
rz(-2.4589804) q[2];
sx q[2];
rz(0.35161099) q[2];
rz(0.44175276) q[3];
sx q[3];
rz(-0.92919246) q[3];
sx q[3];
rz(0.15171224) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0996967) q[0];
sx q[0];
rz(-1.2956887) q[0];
sx q[0];
rz(-1.9805441) q[0];
rz(-2.4940122) q[1];
sx q[1];
rz(-1.3787965) q[1];
sx q[1];
rz(0.82239282) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8063342) q[0];
sx q[0];
rz(-2.0762692) q[0];
sx q[0];
rz(-2.6837206) q[0];
x q[1];
rz(2.4057503) q[2];
sx q[2];
rz(-1.9050361) q[2];
sx q[2];
rz(1.4683958) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31588512) q[1];
sx q[1];
rz(-1.2824821) q[1];
sx q[1];
rz(2.350239) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39607285) q[3];
sx q[3];
rz(-1.1669945) q[3];
sx q[3];
rz(-1.0084821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9245727) q[2];
sx q[2];
rz(-1.1447039) q[2];
sx q[2];
rz(-1.8355969) q[2];
rz(2.4235587) q[3];
sx q[3];
rz(-2.3170203) q[3];
sx q[3];
rz(-3.0927299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43593916) q[0];
sx q[0];
rz(-2.0581364) q[0];
sx q[0];
rz(3.1245681) q[0];
rz(-1.9450933) q[1];
sx q[1];
rz(-2.7462609) q[1];
sx q[1];
rz(-2.8065525) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3075313) q[0];
sx q[0];
rz(-2.1387707) q[0];
sx q[0];
rz(-0.012004367) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4990988) q[2];
sx q[2];
rz(-1.2351002) q[2];
sx q[2];
rz(0.030980008) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10427231) q[1];
sx q[1];
rz(-1.3824711) q[1];
sx q[1];
rz(0.89481797) q[1];
x q[2];
rz(-2.3565759) q[3];
sx q[3];
rz(-2.5015273) q[3];
sx q[3];
rz(-2.9736116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42387858) q[2];
sx q[2];
rz(-2.3953891) q[2];
sx q[2];
rz(-0.27900532) q[2];
rz(-1.772607) q[3];
sx q[3];
rz(-1.5149346) q[3];
sx q[3];
rz(-0.92931187) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5170711) q[0];
sx q[0];
rz(-2.9947424) q[0];
sx q[0];
rz(0.76706925) q[0];
rz(-2.7686139) q[1];
sx q[1];
rz(-1.4992799) q[1];
sx q[1];
rz(0.17072089) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6304008) q[0];
sx q[0];
rz(-2.3053279) q[0];
sx q[0];
rz(1.5319351) q[0];
x q[1];
rz(-1.4095979) q[2];
sx q[2];
rz(-1.1606163) q[2];
sx q[2];
rz(0.31482492) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4258044) q[1];
sx q[1];
rz(-1.8192776) q[1];
sx q[1];
rz(1.2475509) q[1];
rz(0.46782084) q[3];
sx q[3];
rz(-1.7138283) q[3];
sx q[3];
rz(-2.9506872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3501222) q[2];
sx q[2];
rz(-0.56362027) q[2];
sx q[2];
rz(-0.0027837022) q[2];
rz(2.2400098) q[3];
sx q[3];
rz(-1.0222579) q[3];
sx q[3];
rz(-0.80488747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053454178) q[0];
sx q[0];
rz(-2.8202941) q[0];
sx q[0];
rz(1.1711076) q[0];
rz(-2.3114655) q[1];
sx q[1];
rz(-1.8142038) q[1];
sx q[1];
rz(1.289191) q[1];
rz(-0.51132497) q[2];
sx q[2];
rz(-2.1187388) q[2];
sx q[2];
rz(-0.50730898) q[2];
rz(-2.3385901) q[3];
sx q[3];
rz(-1.0754801) q[3];
sx q[3];
rz(0.51343244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
