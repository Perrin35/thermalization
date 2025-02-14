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
rz(-2.6920707) q[0];
sx q[0];
rz(-1.7188526) q[0];
sx q[0];
rz(2.0916405) q[0];
rz(2.6226251) q[1];
sx q[1];
rz(-1.608404) q[1];
sx q[1];
rz(3.0906711) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1919132) q[0];
sx q[0];
rz(-2.5891719) q[0];
sx q[0];
rz(2.4925553) q[0];
rz(2.9257183) q[2];
sx q[2];
rz(-1.2435438) q[2];
sx q[2];
rz(0.33282166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4733153) q[1];
sx q[1];
rz(-0.44679579) q[1];
sx q[1];
rz(0.75832851) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2818358) q[3];
sx q[3];
rz(-0.70631344) q[3];
sx q[3];
rz(3.0970124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8672987) q[2];
sx q[2];
rz(-2.2653502) q[2];
sx q[2];
rz(-1.743861) q[2];
rz(-0.45462576) q[3];
sx q[3];
rz(-0.8780829) q[3];
sx q[3];
rz(1.4286058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6468663) q[0];
sx q[0];
rz(-2.2541101) q[0];
sx q[0];
rz(2.5383762) q[0];
rz(-1.8869205) q[1];
sx q[1];
rz(-2.5277977) q[1];
sx q[1];
rz(-2.5616554) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7042687) q[0];
sx q[0];
rz(-0.56341972) q[0];
sx q[0];
rz(-2.826344) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5712319) q[2];
sx q[2];
rz(-1.959539) q[2];
sx q[2];
rz(-2.027957) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9442288) q[1];
sx q[1];
rz(-1.6104638) q[1];
sx q[1];
rz(0.96645379) q[1];
x q[2];
rz(0.5671366) q[3];
sx q[3];
rz(-1.2492164) q[3];
sx q[3];
rz(-2.0547607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4846399) q[2];
sx q[2];
rz(-0.75289774) q[2];
sx q[2];
rz(-2.7638655) q[2];
rz(1.7083302) q[3];
sx q[3];
rz(-1.6323615) q[3];
sx q[3];
rz(-1.1908971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3483873) q[0];
sx q[0];
rz(-3.0866525) q[0];
sx q[0];
rz(1.6435664) q[0];
rz(1.9801697) q[1];
sx q[1];
rz(-1.8893416) q[1];
sx q[1];
rz(-0.84322554) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0204165) q[0];
sx q[0];
rz(-1.6303697) q[0];
sx q[0];
rz(2.1942433) q[0];
rz(-pi) q[1];
rz(2.9239028) q[2];
sx q[2];
rz(-1.9799616) q[2];
sx q[2];
rz(-0.37693757) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.49065152) q[1];
sx q[1];
rz(-0.89812121) q[1];
sx q[1];
rz(2.9789799) q[1];
rz(2.7494861) q[3];
sx q[3];
rz(-2.4014353) q[3];
sx q[3];
rz(1.8405434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8680385) q[2];
sx q[2];
rz(-0.66323438) q[2];
sx q[2];
rz(-0.794945) q[2];
rz(-2.3718209) q[3];
sx q[3];
rz(-2.218518) q[3];
sx q[3];
rz(-2.9845089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(1.5539826) q[0];
sx q[0];
rz(-1.8164604) q[0];
sx q[0];
rz(0.28833589) q[0];
rz(-2.4544857) q[1];
sx q[1];
rz(-1.4930875) q[1];
sx q[1];
rz(-1.7126602) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16609678) q[0];
sx q[0];
rz(-1.5642484) q[0];
sx q[0];
rz(-3.1295144) q[0];
x q[1];
rz(2.7393514) q[2];
sx q[2];
rz(-1.7000442) q[2];
sx q[2];
rz(2.8856173) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7894788) q[1];
sx q[1];
rz(-1.6360456) q[1];
sx q[1];
rz(-0.25635527) q[1];
rz(-2.5023584) q[3];
sx q[3];
rz(-2.1393288) q[3];
sx q[3];
rz(-0.45445874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.571542) q[2];
sx q[2];
rz(-0.96633458) q[2];
sx q[2];
rz(-0.16560444) q[2];
rz(-0.064519493) q[3];
sx q[3];
rz(-3.0118628) q[3];
sx q[3];
rz(-1.8701514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22873838) q[0];
sx q[0];
rz(-1.1930635) q[0];
sx q[0];
rz(-2.2304529) q[0];
rz(2.1707824) q[1];
sx q[1];
rz(-1.3263005) q[1];
sx q[1];
rz(2.2784746) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9476833) q[0];
sx q[0];
rz(-1.2811617) q[0];
sx q[0];
rz(-1.0632443) q[0];
x q[1];
rz(1.3023071) q[2];
sx q[2];
rz(-1.5932065) q[2];
sx q[2];
rz(-2.7145998) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99030322) q[1];
sx q[1];
rz(-2.5545729) q[1];
sx q[1];
rz(2.7128316) q[1];
x q[2];
rz(1.4310525) q[3];
sx q[3];
rz(-0.31523963) q[3];
sx q[3];
rz(-1.8342575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0866278) q[2];
sx q[2];
rz(-1.2530155) q[2];
sx q[2];
rz(0.48913726) q[2];
rz(-2.1431811) q[3];
sx q[3];
rz(-2.9676134) q[3];
sx q[3];
rz(3.0424931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973706) q[0];
sx q[0];
rz(-0.64592823) q[0];
sx q[0];
rz(-0.55322629) q[0];
rz(-0.79752254) q[1];
sx q[1];
rz(-0.99634606) q[1];
sx q[1];
rz(-1.1513938) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0259375) q[0];
sx q[0];
rz(-2.135072) q[0];
sx q[0];
rz(1.7274117) q[0];
rz(1.3998447) q[2];
sx q[2];
rz(-1.0264077) q[2];
sx q[2];
rz(-2.483727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0480736) q[1];
sx q[1];
rz(-2.8052605) q[1];
sx q[1];
rz(-2.5079355) q[1];
x q[2];
rz(2.0460412) q[3];
sx q[3];
rz(-2.4184347) q[3];
sx q[3];
rz(-2.1840546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79661757) q[2];
sx q[2];
rz(-1.7285708) q[2];
sx q[2];
rz(-0.83149347) q[2];
rz(1.8031395) q[3];
sx q[3];
rz(-1.0748539) q[3];
sx q[3];
rz(-0.89053806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(-1.3727342) q[0];
sx q[0];
rz(-1.2059728) q[0];
sx q[0];
rz(-0.15705577) q[0];
rz(1.318469) q[1];
sx q[1];
rz(-0.942197) q[1];
sx q[1];
rz(2.7838321) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45159949) q[0];
sx q[0];
rz(-2.4781961) q[0];
sx q[0];
rz(-2.3729352) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6951201) q[2];
sx q[2];
rz(-2.1887767) q[2];
sx q[2];
rz(-0.22836049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50510943) q[1];
sx q[1];
rz(-0.92527522) q[1];
sx q[1];
rz(-2.4277975) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5671183) q[3];
sx q[3];
rz(-1.7041429) q[3];
sx q[3];
rz(2.3337618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99178189) q[2];
sx q[2];
rz(-1.8842183) q[2];
sx q[2];
rz(-3.120976) q[2];
rz(1.2425544) q[3];
sx q[3];
rz(-0.81762448) q[3];
sx q[3];
rz(-0.29153618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5307411) q[0];
sx q[0];
rz(-1.956097) q[0];
sx q[0];
rz(-0.32671842) q[0];
rz(-2.2945981) q[1];
sx q[1];
rz(-1.0495443) q[1];
sx q[1];
rz(1.0583896) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055870067) q[0];
sx q[0];
rz(-0.51858178) q[0];
sx q[0];
rz(-1.7118171) q[0];
rz(2.9061163) q[2];
sx q[2];
rz(-2.4889766) q[2];
sx q[2];
rz(2.1760698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10524532) q[1];
sx q[1];
rz(-0.39064841) q[1];
sx q[1];
rz(-1.3643144) q[1];
rz(-pi) q[2];
rz(2.8260569) q[3];
sx q[3];
rz(-2.1140751) q[3];
sx q[3];
rz(1.735294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1191001) q[2];
sx q[2];
rz(-2.7969226) q[2];
sx q[2];
rz(-1.9692839) q[2];
rz(3.1398224) q[3];
sx q[3];
rz(-2.6383196) q[3];
sx q[3];
rz(-2.1625904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2718647) q[0];
sx q[0];
rz(-2.3379022) q[0];
sx q[0];
rz(1.8029689) q[0];
rz(-1.9610693) q[1];
sx q[1];
rz(-2.0484643) q[1];
sx q[1];
rz(0.48042935) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051078253) q[0];
sx q[0];
rz(-1.773293) q[0];
sx q[0];
rz(-2.0351054) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17002885) q[2];
sx q[2];
rz(-1.5289834) q[2];
sx q[2];
rz(-2.5332076) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1611966) q[1];
sx q[1];
rz(-2.5727486) q[1];
sx q[1];
rz(-1.7561164) q[1];
rz(0.35828405) q[3];
sx q[3];
rz(-1.7833774) q[3];
sx q[3];
rz(-2.6019118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7052445) q[2];
sx q[2];
rz(-2.1465325) q[2];
sx q[2];
rz(-2.056541) q[2];
rz(-2.0294225) q[3];
sx q[3];
rz(-2.2472436) q[3];
sx q[3];
rz(-0.81356847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28806624) q[0];
sx q[0];
rz(-0.94678322) q[0];
sx q[0];
rz(-1.0119525) q[0];
rz(0.90156737) q[1];
sx q[1];
rz(-2.8755867) q[1];
sx q[1];
rz(0.56545767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0386977) q[0];
sx q[0];
rz(-0.92864543) q[0];
sx q[0];
rz(-0.26102926) q[0];
rz(-1.9996634) q[2];
sx q[2];
rz(-1.8151917) q[2];
sx q[2];
rz(-1.3809057) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9931075) q[1];
sx q[1];
rz(-2.4186181) q[1];
sx q[1];
rz(-1.430368) q[1];
rz(-3.0363085) q[3];
sx q[3];
rz(-1.255688) q[3];
sx q[3];
rz(-1.1694825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7908343) q[2];
sx q[2];
rz(-1.5956722) q[2];
sx q[2];
rz(-1.8809543) q[2];
rz(-2.6993921) q[3];
sx q[3];
rz(-1.0298157) q[3];
sx q[3];
rz(1.7650167) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5762536) q[0];
sx q[0];
rz(-2.2591142) q[0];
sx q[0];
rz(2.2478065) q[0];
rz(-1.5743938) q[1];
sx q[1];
rz(-1.6849453) q[1];
sx q[1];
rz(-1.4243855) q[1];
rz(-2.7161063) q[2];
sx q[2];
rz(-1.8271227) q[2];
sx q[2];
rz(2.9572077) q[2];
rz(-1.1461729) q[3];
sx q[3];
rz(-0.38214798) q[3];
sx q[3];
rz(1.8761763) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
