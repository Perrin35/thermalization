OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(3.6448195) q[0];
sx q[0];
rz(10.148944) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(-2.35676) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0137579) q[0];
sx q[0];
rz(-2.7779967) q[0];
sx q[0];
rz(-2.512393) q[0];
rz(-pi) q[1];
rz(0.41854026) q[2];
sx q[2];
rz(-1.6633908) q[2];
sx q[2];
rz(1.6469524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19476328) q[1];
sx q[1];
rz(-2.0954872) q[1];
sx q[1];
rz(0.25804934) q[1];
rz(-0.033693245) q[3];
sx q[3];
rz(-1.9851306) q[3];
sx q[3];
rz(1.8309483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5518387) q[2];
sx q[2];
rz(-1.4171615) q[2];
sx q[2];
rz(-0.067967728) q[2];
rz(3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(-1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2215866) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(0.24366972) q[0];
rz(-0.63175732) q[1];
sx q[1];
rz(-1.4032204) q[1];
sx q[1];
rz(1.3557281) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6512017) q[0];
sx q[0];
rz(-2.8848007) q[0];
sx q[0];
rz(-1.4635565) q[0];
rz(-pi) q[1];
rz(-1.182105) q[2];
sx q[2];
rz(-1.9027862) q[2];
sx q[2];
rz(-1.4093083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.3416459) q[1];
sx q[1];
rz(-0.46657545) q[1];
sx q[1];
rz(-1.8599618) q[1];
rz(1.5242819) q[3];
sx q[3];
rz(-2.1624613) q[3];
sx q[3];
rz(0.072629645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(2.8919354) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(2.8095968) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24519414) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(2.2431592) q[0];
rz(-1.8067182) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(1.2737087) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0231409) q[0];
sx q[0];
rz(-1.3190735) q[0];
sx q[0];
rz(-1.203042) q[0];
rz(2.8461371) q[2];
sx q[2];
rz(-1.6647415) q[2];
sx q[2];
rz(-0.23677793) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0395567) q[1];
sx q[1];
rz(-1.2219056) q[1];
sx q[1];
rz(-0.7015014) q[1];
rz(-pi) q[2];
rz(2.5898315) q[3];
sx q[3];
rz(-2.3380087) q[3];
sx q[3];
rz(1.4078275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0818103) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(0.95345062) q[2];
rz(3.1070784) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(-2.9147193) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8614486) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(-2.8934073) q[0];
rz(2.10363) q[1];
sx q[1];
rz(-1.1231517) q[1];
sx q[1];
rz(3.0674556) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4081057) q[0];
sx q[0];
rz(-1.8251849) q[0];
sx q[0];
rz(2.0612201) q[0];
rz(0.89216994) q[2];
sx q[2];
rz(-1.2954419) q[2];
sx q[2];
rz(2.8868669) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76064199) q[1];
sx q[1];
rz(-1.4529072) q[1];
sx q[1];
rz(-1.8878493) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7449964) q[3];
sx q[3];
rz(-1.3341691) q[3];
sx q[3];
rz(3.0803806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(0.038643535) q[2];
rz(-0.97366992) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(0.27004778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(-1.779153) q[0];
rz(0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(1.1626676) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2813331) q[0];
sx q[0];
rz(-1.5409924) q[0];
sx q[0];
rz(-1.6822862) q[0];
rz(-pi) q[1];
rz(-0.16935279) q[2];
sx q[2];
rz(-1.2522109) q[2];
sx q[2];
rz(-0.40118518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77046493) q[1];
sx q[1];
rz(-2.7264997) q[1];
sx q[1];
rz(2.0626555) q[1];
x q[2];
rz(-2.0883457) q[3];
sx q[3];
rz(-1.187511) q[3];
sx q[3];
rz(2.7305207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5148619) q[2];
sx q[2];
rz(-1.0802439) q[2];
sx q[2];
rz(-2.999372) q[2];
rz(-0.90406117) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267923) q[0];
sx q[0];
rz(-0.27652201) q[0];
sx q[0];
rz(1.4676771) q[0];
rz(-2.5698075) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(-2.8335559) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338617) q[0];
sx q[0];
rz(-2.9836285) q[0];
sx q[0];
rz(2.4979742) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1790444) q[2];
sx q[2];
rz(-1.0481917) q[2];
sx q[2];
rz(2.7163497) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6412515) q[1];
sx q[1];
rz(-0.89741035) q[1];
sx q[1];
rz(0.25232368) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8940053) q[3];
sx q[3];
rz(-0.35841225) q[3];
sx q[3];
rz(0.7022411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77928153) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(1.9936838) q[2];
rz(2.4273196) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(2.4961297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66184735) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(-0.1299783) q[0];
rz(0.030844363) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(-2.470509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2790047) q[0];
sx q[0];
rz(-0.064084856) q[0];
sx q[0];
rz(2.5628753) q[0];
x q[1];
rz(1.7550049) q[2];
sx q[2];
rz(-1.6008018) q[2];
sx q[2];
rz(-1.1785933) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8819067) q[1];
sx q[1];
rz(-0.38609186) q[1];
sx q[1];
rz(-1.8465471) q[1];
rz(-3.0942261) q[3];
sx q[3];
rz(-2.2455375) q[3];
sx q[3];
rz(2.588152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0188296) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(2.5637131) q[2];
rz(3.1130062) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9776483) q[0];
sx q[0];
rz(-2.2387235) q[0];
sx q[0];
rz(0.41241616) q[0];
rz(1.4498129) q[1];
sx q[1];
rz(-1.342536) q[1];
sx q[1];
rz(-1.1669881) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3229423) q[0];
sx q[0];
rz(-0.65781051) q[0];
sx q[0];
rz(1.9239182) q[0];
rz(-pi) q[1];
rz(1.251986) q[2];
sx q[2];
rz(-2.406771) q[2];
sx q[2];
rz(-0.97359818) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3762714) q[1];
sx q[1];
rz(-1.6011392) q[1];
sx q[1];
rz(0.18860753) q[1];
rz(1.0409045) q[3];
sx q[3];
rz(-1.7003254) q[3];
sx q[3];
rz(1.1216175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2010487) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-0.18903014) q[2];
rz(0.14686251) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(-1.3930901) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1259595) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(0.92700672) q[0];
rz(-1.758763) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(1.6519201) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72367523) q[0];
sx q[0];
rz(-1.105504) q[0];
sx q[0];
rz(-0.76405163) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.050612014) q[2];
sx q[2];
rz(-1.0992556) q[2];
sx q[2];
rz(-2.5052349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8330194) q[1];
sx q[1];
rz(-0.25990572) q[1];
sx q[1];
rz(-0.72867568) q[1];
x q[2];
rz(-2.3143523) q[3];
sx q[3];
rz(-2.4628371) q[3];
sx q[3];
rz(2.7620706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5513409) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(-1.7112973) q[2];
rz(-0.57724214) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(-1.757471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5883314) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(-2.8531895) q[0];
rz(-0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(-0.14702252) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0518236) q[0];
sx q[0];
rz(-2.4252486) q[0];
sx q[0];
rz(-2.1728974) q[0];
rz(-1.64098) q[2];
sx q[2];
rz(-1.8101705) q[2];
sx q[2];
rz(-1.0603051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8204931) q[1];
sx q[1];
rz(-1.610678) q[1];
sx q[1];
rz(-0.95221968) q[1];
rz(0.57659984) q[3];
sx q[3];
rz(-1.7017662) q[3];
sx q[3];
rz(3.1104345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0599351) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(0.74404136) q[2];
rz(0.75731164) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(-1.3967167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.025678) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(0.81746447) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(-1.5031917) q[2];
sx q[2];
rz(-1.1195782) q[2];
sx q[2];
rz(-1.3200214) q[2];
rz(1.0088624) q[3];
sx q[3];
rz(-1.6812134) q[3];
sx q[3];
rz(-2.5440661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
