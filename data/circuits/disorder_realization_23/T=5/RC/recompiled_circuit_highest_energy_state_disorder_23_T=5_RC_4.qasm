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
rz(0.55573207) q[0];
sx q[0];
rz(-1.8620123) q[0];
sx q[0];
rz(-0.32787856) q[0];
rz(-2.9887587) q[1];
sx q[1];
rz(-2.6522377) q[1];
sx q[1];
rz(2.1305003) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41844925) q[0];
sx q[0];
rz(-1.4108424) q[0];
sx q[0];
rz(2.5470995) q[0];
rz(-pi) q[1];
rz(2.7332441) q[2];
sx q[2];
rz(-1.7350001) q[2];
sx q[2];
rz(-2.4105031) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7465621) q[1];
sx q[1];
rz(-0.85828188) q[1];
sx q[1];
rz(1.0814352) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23787161) q[3];
sx q[3];
rz(-1.9783881) q[3];
sx q[3];
rz(-1.067124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27935394) q[2];
sx q[2];
rz(-2.3591159) q[2];
sx q[2];
rz(-1.5834825) q[2];
rz(-2.8067348) q[3];
sx q[3];
rz(-2.0575276) q[3];
sx q[3];
rz(-0.28086942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25664499) q[0];
sx q[0];
rz(-2.5769951) q[0];
sx q[0];
rz(2.8096492) q[0];
rz(-0.360082) q[1];
sx q[1];
rz(-1.8372476) q[1];
sx q[1];
rz(-2.8538381) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1126039) q[0];
sx q[0];
rz(-1.3224012) q[0];
sx q[0];
rz(-1.364014) q[0];
rz(-pi) q[1];
rz(-1.8910725) q[2];
sx q[2];
rz(-0.8647635) q[2];
sx q[2];
rz(1.4780188) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33720484) q[1];
sx q[1];
rz(-1.8998977) q[1];
sx q[1];
rz(1.9371402) q[1];
rz(-pi) q[2];
rz(2.5310181) q[3];
sx q[3];
rz(-2.861851) q[3];
sx q[3];
rz(2.4795462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.45431367) q[2];
sx q[2];
rz(-1.6087029) q[2];
sx q[2];
rz(-1.7506556) q[2];
rz(-0.85919291) q[3];
sx q[3];
rz(-1.2733368) q[3];
sx q[3];
rz(0.19101645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39621064) q[0];
sx q[0];
rz(-1.5353545) q[0];
sx q[0];
rz(2.9111653) q[0];
rz(2.6241809) q[1];
sx q[1];
rz(-1.0942065) q[1];
sx q[1];
rz(2.8527625) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6073501) q[0];
sx q[0];
rz(-1.6301883) q[0];
sx q[0];
rz(1.1618105) q[0];
x q[1];
rz(-1.4963228) q[2];
sx q[2];
rz(-1.1033325) q[2];
sx q[2];
rz(2.4592768) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3687262) q[1];
sx q[1];
rz(-1.9802367) q[1];
sx q[1];
rz(0.057842908) q[1];
rz(-2.9656319) q[3];
sx q[3];
rz(-1.7909414) q[3];
sx q[3];
rz(-0.99316059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1032224) q[2];
sx q[2];
rz(-2.4028845) q[2];
sx q[2];
rz(-0.040741097) q[2];
rz(-0.1013969) q[3];
sx q[3];
rz(-1.3359759) q[3];
sx q[3];
rz(2.0749157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3876225) q[0];
sx q[0];
rz(-0.0059703537) q[0];
sx q[0];
rz(-2.907584) q[0];
rz(0.19800828) q[1];
sx q[1];
rz(-1.0375236) q[1];
sx q[1];
rz(2.1580946) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6925544) q[0];
sx q[0];
rz(-2.6561167) q[0];
sx q[0];
rz(-1.6621157) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6621236) q[2];
sx q[2];
rz(-1.0589561) q[2];
sx q[2];
rz(0.78081607) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2605839) q[1];
sx q[1];
rz(-0.99956341) q[1];
sx q[1];
rz(-0.84239475) q[1];
x q[2];
rz(-1.4928905) q[3];
sx q[3];
rz(-0.39038218) q[3];
sx q[3];
rz(1.4161033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46821758) q[2];
sx q[2];
rz(-2.8523291) q[2];
sx q[2];
rz(2.3667228) q[2];
rz(-1.3261999) q[3];
sx q[3];
rz(-1.1893136) q[3];
sx q[3];
rz(-2.6302122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73959094) q[0];
sx q[0];
rz(-2.2690161) q[0];
sx q[0];
rz(-3.1091029) q[0];
rz(-1.912311) q[1];
sx q[1];
rz(-2.2132497) q[1];
sx q[1];
rz(0.083018735) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5755324) q[0];
sx q[0];
rz(-2.0216746) q[0];
sx q[0];
rz(-1.321248) q[0];
rz(1.4013702) q[2];
sx q[2];
rz(-1.8352785) q[2];
sx q[2];
rz(-1.5096017) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8889035) q[1];
sx q[1];
rz(-1.1071883) q[1];
sx q[1];
rz(-2.3867334) q[1];
rz(2.2076603) q[3];
sx q[3];
rz(-1.5832925) q[3];
sx q[3];
rz(-3.0325471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7897537) q[2];
sx q[2];
rz(-2.3408076) q[2];
sx q[2];
rz(0.16072533) q[2];
rz(-1.9488581) q[3];
sx q[3];
rz(-0.21218097) q[3];
sx q[3];
rz(0.76782697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.47650325) q[0];
sx q[0];
rz(-1.1192717) q[0];
sx q[0];
rz(0.29092586) q[0];
rz(1.3833969) q[1];
sx q[1];
rz(-1.8427126) q[1];
sx q[1];
rz(-2.6417522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6747777) q[0];
sx q[0];
rz(-1.6143911) q[0];
sx q[0];
rz(1.6131667) q[0];
rz(-pi) q[1];
rz(0.96237125) q[2];
sx q[2];
rz(-0.82686868) q[2];
sx q[2];
rz(-0.27368557) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9466797) q[1];
sx q[1];
rz(-2.1482067) q[1];
sx q[1];
rz(-1.0275847) q[1];
rz(-pi) q[2];
rz(2.1736702) q[3];
sx q[3];
rz(-0.80652666) q[3];
sx q[3];
rz(1.4462808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9342039) q[2];
sx q[2];
rz(-2.5298205) q[2];
sx q[2];
rz(-0.70154166) q[2];
rz(-0.97405854) q[3];
sx q[3];
rz(-2.3249224) q[3];
sx q[3];
rz(2.2402703) q[3];
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
rz(-1.8262254) q[0];
sx q[0];
rz(-0.81804818) q[0];
sx q[0];
rz(0.65959626) q[0];
rz(1.9024128) q[1];
sx q[1];
rz(-1.0737123) q[1];
sx q[1];
rz(2.0674131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.488193) q[0];
sx q[0];
rz(-1.5563413) q[0];
sx q[0];
rz(-0.01343347) q[0];
rz(-0.72788179) q[2];
sx q[2];
rz(-2.5932576) q[2];
sx q[2];
rz(-0.20073433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3201808) q[1];
sx q[1];
rz(-1.786725) q[1];
sx q[1];
rz(1.126299) q[1];
x q[2];
rz(2.3775565) q[3];
sx q[3];
rz(-1.2437399) q[3];
sx q[3];
rz(-0.59248176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85019511) q[2];
sx q[2];
rz(-2.3829298) q[2];
sx q[2];
rz(0.58471739) q[2];
rz(1.0662474) q[3];
sx q[3];
rz(-1.2277579) q[3];
sx q[3];
rz(2.7339981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.9503815) q[0];
sx q[0];
rz(-1.9712912) q[0];
sx q[0];
rz(2.6608652) q[0];
rz(1.2113384) q[1];
sx q[1];
rz(-1.6119266) q[1];
sx q[1];
rz(2.5239677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6473501) q[0];
sx q[0];
rz(-1.2456919) q[0];
sx q[0];
rz(-3.1371497) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4489538) q[2];
sx q[2];
rz(-2.4071781) q[2];
sx q[2];
rz(-1.8633597) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9165914) q[1];
sx q[1];
rz(-1.3559506) q[1];
sx q[1];
rz(1.8720354) q[1];
rz(2.9493773) q[3];
sx q[3];
rz(-1.6055067) q[3];
sx q[3];
rz(-2.672003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1787313) q[2];
sx q[2];
rz(-1.7469254) q[2];
sx q[2];
rz(-0.89087957) q[2];
rz(-1.0293055) q[3];
sx q[3];
rz(-1.7512713) q[3];
sx q[3];
rz(1.5911969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35128281) q[0];
sx q[0];
rz(-2.2064378) q[0];
sx q[0];
rz(-1.1736322) q[0];
rz(1.952518) q[1];
sx q[1];
rz(-0.74780858) q[1];
sx q[1];
rz(2.660451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0521942) q[0];
sx q[0];
rz(-1.3476831) q[0];
sx q[0];
rz(2.5714178) q[0];
rz(-pi) q[1];
rz(0.16002197) q[2];
sx q[2];
rz(-0.56567398) q[2];
sx q[2];
rz(-0.46505798) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6509241) q[1];
sx q[1];
rz(-1.4238402) q[1];
sx q[1];
rz(2.1747194) q[1];
rz(-0.88925006) q[3];
sx q[3];
rz(-2.4947531) q[3];
sx q[3];
rz(-1.3726039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9159307) q[2];
sx q[2];
rz(-1.2292726) q[2];
sx q[2];
rz(1.6602328) q[2];
rz(-3.0952752) q[3];
sx q[3];
rz(-1.1888844) q[3];
sx q[3];
rz(1.2272629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34761053) q[0];
sx q[0];
rz(-1.0337669) q[0];
sx q[0];
rz(3.0718497) q[0];
rz(-2.8522885) q[1];
sx q[1];
rz(-1.8569088) q[1];
sx q[1];
rz(2.0893673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493794) q[0];
sx q[0];
rz(-1.2468296) q[0];
sx q[0];
rz(2.762336) q[0];
x q[1];
rz(-1.7830816) q[2];
sx q[2];
rz(-2.4683778) q[2];
sx q[2];
rz(-0.5897943) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6695646) q[1];
sx q[1];
rz(-1.2109562) q[1];
sx q[1];
rz(-1.9739975) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7131548) q[3];
sx q[3];
rz(-1.1827381) q[3];
sx q[3];
rz(-0.65754277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5436486) q[2];
sx q[2];
rz(-1.1755377) q[2];
sx q[2];
rz(3.0808466) q[2];
rz(0.28901035) q[3];
sx q[3];
rz(-1.3649536) q[3];
sx q[3];
rz(1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.65171) q[0];
sx q[0];
rz(-1.6429506) q[0];
sx q[0];
rz(2.0605675) q[0];
rz(-1.0501077) q[1];
sx q[1];
rz(-2.0790015) q[1];
sx q[1];
rz(1.9008295) q[1];
rz(-2.0403258) q[2];
sx q[2];
rz(-1.4360089) q[2];
sx q[2];
rz(0.77965005) q[2];
rz(-2.859194) q[3];
sx q[3];
rz(-1.9871319) q[3];
sx q[3];
rz(2.7642767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
