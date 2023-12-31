OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4857793) q[0];
sx q[0];
rz(-1.9137148) q[0];
sx q[0];
rz(-1.2113843) q[0];
x q[1];
rz(2.3762796) q[2];
sx q[2];
rz(-1.0368477) q[2];
sx q[2];
rz(-3.090976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1882602) q[1];
sx q[1];
rz(-0.42020513) q[1];
sx q[1];
rz(0.62700595) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2075495) q[3];
sx q[3];
rz(-1.7540635) q[3];
sx q[3];
rz(-2.7024384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2661665) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(-2.0092633) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(0.18584132) q[0];
rz(0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(0.21683189) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8685535) q[0];
sx q[0];
rz(-0.93398636) q[0];
sx q[0];
rz(2.7815232) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0160604) q[2];
sx q[2];
rz(-1.1604571) q[2];
sx q[2];
rz(1.0085269) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9456957) q[1];
sx q[1];
rz(-2.2356114) q[1];
sx q[1];
rz(2.871454) q[1];
rz(-pi) q[2];
rz(2.300755) q[3];
sx q[3];
rz(-2.3008122) q[3];
sx q[3];
rz(-2.9527612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(-1.8537834) q[2];
rz(0.76256049) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(-1.8151981) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-0.93260971) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1754477) q[0];
sx q[0];
rz(-1.6245337) q[0];
sx q[0];
rz(1.2530112) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2432125) q[2];
sx q[2];
rz(-0.2873688) q[2];
sx q[2];
rz(-0.76295602) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31308094) q[1];
sx q[1];
rz(-2.5701437) q[1];
sx q[1];
rz(0.22027318) q[1];
x q[2];
rz(-0.91498418) q[3];
sx q[3];
rz(-2.4768156) q[3];
sx q[3];
rz(-1.8613601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(-2.0489342) q[2];
rz(-2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(2.6413667) q[0];
rz(2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(1.6436228) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.995979) q[0];
sx q[0];
rz(-1.0756452) q[0];
sx q[0];
rz(2.8061295) q[0];
rz(1.7613212) q[2];
sx q[2];
rz(-1.6263783) q[2];
sx q[2];
rz(1.3560825) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8151617) q[1];
sx q[1];
rz(-1.8306499) q[1];
sx q[1];
rz(1.8111147) q[1];
x q[2];
rz(-1.1180531) q[3];
sx q[3];
rz(-1.1467883) q[3];
sx q[3];
rz(2.8297569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(-0.70181075) q[2];
rz(-2.3102405) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(-2.519616) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410626) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(2.3262614) q[0];
rz(-1.5218081) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(-2.0933847) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5883023) q[0];
sx q[0];
rz(-2.7748845) q[0];
sx q[0];
rz(-2.328863) q[0];
x q[1];
rz(-3.1254966) q[2];
sx q[2];
rz(-2.490009) q[2];
sx q[2];
rz(-0.96166699) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9427467) q[1];
sx q[1];
rz(-1.8837351) q[1];
sx q[1];
rz(1.8206157) q[1];
rz(-0.95543315) q[3];
sx q[3];
rz(-1.9121998) q[3];
sx q[3];
rz(-2.0379025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0734171) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(0.90240479) q[0];
rz(-1.0166608) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(0.12983233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5141402) q[0];
sx q[0];
rz(-0.99533004) q[0];
sx q[0];
rz(2.0155725) q[0];
rz(-pi) q[1];
rz(0.99545698) q[2];
sx q[2];
rz(-1.2577004) q[2];
sx q[2];
rz(-2.7196333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9784769) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(-2.5549868) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.767166) q[3];
sx q[3];
rz(-1.0305627) q[3];
sx q[3];
rz(-1.2353209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(0.20425805) q[2];
rz(-1.2060818) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181353) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(-1.4468505) q[0];
rz(-1.2591259) q[1];
sx q[1];
rz(-2.1513758) q[1];
sx q[1];
rz(-2.4553305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9335564) q[0];
sx q[0];
rz(-0.66970034) q[0];
sx q[0];
rz(2.3963388) q[0];
rz(-pi) q[1];
rz(-0.83001901) q[2];
sx q[2];
rz(-2.0247211) q[2];
sx q[2];
rz(1.9759076) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2974907) q[1];
sx q[1];
rz(-2.9417848) q[1];
sx q[1];
rz(1.1872477) q[1];
x q[2];
rz(2.7191914) q[3];
sx q[3];
rz(-1.8654612) q[3];
sx q[3];
rz(0.15448031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4454322) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(-0.0017722842) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(-1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6034265) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(-0.7810477) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(-1.5015645) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43270375) q[0];
sx q[0];
rz(-0.49312691) q[0];
sx q[0];
rz(-1.5296442) q[0];
rz(1.976622) q[2];
sx q[2];
rz(-0.067194447) q[2];
sx q[2];
rz(-0.42025987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1135243) q[1];
sx q[1];
rz(-1.4816195) q[1];
sx q[1];
rz(2.813617) q[1];
rz(-pi) q[2];
rz(-0.99174188) q[3];
sx q[3];
rz(-1.0240882) q[3];
sx q[3];
rz(-2.9255097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7897196) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(1.8224576) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(-0.31931988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(1.2040899) q[0];
rz(-0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(0.35167545) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4597804) q[0];
sx q[0];
rz(-0.60428719) q[0];
sx q[0];
rz(-0.87862815) q[0];
rz(-2.7792764) q[2];
sx q[2];
rz(-3*pi/13) q[2];
sx q[2];
rz(2.97646) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.44605276) q[1];
sx q[1];
rz(-2.3791168) q[1];
sx q[1];
rz(-1.4941077) q[1];
rz(-pi) q[2];
rz(-0.44378186) q[3];
sx q[3];
rz(-0.1212596) q[3];
sx q[3];
rz(-1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7982771) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(1.8593672) q[2];
rz(1.4964237) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(-1.6431817) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(1.0378029) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(-2.1077572) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68034222) q[0];
sx q[0];
rz(-0.67078062) q[0];
sx q[0];
rz(-1.3602507) q[0];
x q[1];
rz(2.6331484) q[2];
sx q[2];
rz(-0.93548453) q[2];
sx q[2];
rz(-2.9490162) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.918805) q[1];
sx q[1];
rz(-2.5606887) q[1];
sx q[1];
rz(1.3851628) q[1];
rz(-pi) q[2];
rz(2.8909573) q[3];
sx q[3];
rz(-0.72892979) q[3];
sx q[3];
rz(0.38754101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(0.6357843) q[2];
rz(0.27030269) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.6939659) q[0];
sx q[0];
rz(-1.3128558) q[0];
sx q[0];
rz(-2.0679612) q[0];
rz(1.4355961) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(0.96799093) q[2];
sx q[2];
rz(-1.5445166) q[2];
sx q[2];
rz(-1.9894285) q[2];
rz(-1.773949) q[3];
sx q[3];
rz(-2.805134) q[3];
sx q[3];
rz(0.72611879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
