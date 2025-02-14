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
rz(-1.5440829) q[0];
sx q[0];
rz(-1.7680327) q[0];
sx q[0];
rz(-1.0184259) q[0];
rz(2.6234558) q[1];
sx q[1];
rz(5.5261373) q[1];
sx q[1];
rz(11.93461) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3903094) q[0];
sx q[0];
rz(-1.7770808) q[0];
sx q[0];
rz(0.12796107) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1167157) q[2];
sx q[2];
rz(-1.8734249) q[2];
sx q[2];
rz(0.51371117) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2871197) q[1];
sx q[1];
rz(-1.8504256) q[1];
sx q[1];
rz(2.2210712) q[1];
rz(-pi) q[2];
rz(1.3034091) q[3];
sx q[3];
rz(-2.8041556) q[3];
sx q[3];
rz(2.9149559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53039256) q[2];
sx q[2];
rz(-0.35553122) q[2];
sx q[2];
rz(1.6543039) q[2];
rz(-0.90855956) q[3];
sx q[3];
rz(-2.9049951) q[3];
sx q[3];
rz(2.1366185) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0205883) q[0];
sx q[0];
rz(-1.7089184) q[0];
sx q[0];
rz(-2.9793136) q[0];
rz(-0.24457112) q[1];
sx q[1];
rz(-1.2030315) q[1];
sx q[1];
rz(0.29168209) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89842105) q[0];
sx q[0];
rz(-2.7968639) q[0];
sx q[0];
rz(0.60144641) q[0];
rz(-pi) q[1];
rz(-0.9082011) q[2];
sx q[2];
rz(-0.35396265) q[2];
sx q[2];
rz(2.5733054) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1202787) q[1];
sx q[1];
rz(-0.81715423) q[1];
sx q[1];
rz(-1.2517125) q[1];
rz(-2.347885) q[3];
sx q[3];
rz(-2.5849403) q[3];
sx q[3];
rz(1.2964378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4633999) q[2];
sx q[2];
rz(-2.2726111) q[2];
sx q[2];
rz(-1.0758859) q[2];
rz(1.894527) q[3];
sx q[3];
rz(-1.5445292) q[3];
sx q[3];
rz(1.6943078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4397044) q[0];
sx q[0];
rz(-1.1264369) q[0];
sx q[0];
rz(2.9578748) q[0];
rz(1.5048997) q[1];
sx q[1];
rz(-0.74438649) q[1];
sx q[1];
rz(0.16477975) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713584) q[0];
sx q[0];
rz(-0.58588282) q[0];
sx q[0];
rz(1.5878994) q[0];
x q[1];
rz(-0.67586502) q[2];
sx q[2];
rz(-1.8177351) q[2];
sx q[2];
rz(0.20154146) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2307869) q[1];
sx q[1];
rz(-1.1890829) q[1];
sx q[1];
rz(0.64818212) q[1];
rz(-pi) q[2];
rz(2.1042473) q[3];
sx q[3];
rz(-1.8672393) q[3];
sx q[3];
rz(0.25755778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.034721) q[2];
sx q[2];
rz(-1.9811337) q[2];
sx q[2];
rz(-1.1064233) q[2];
rz(-0.70118457) q[3];
sx q[3];
rz(-1.3283575) q[3];
sx q[3];
rz(2.6760694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.3727386) q[0];
sx q[0];
rz(-1.1356069) q[0];
sx q[0];
rz(-0.94451529) q[0];
rz(1.5185897) q[1];
sx q[1];
rz(-2.1606162) q[1];
sx q[1];
rz(-1.5400344) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2006783) q[0];
sx q[0];
rz(-0.8884065) q[0];
sx q[0];
rz(-0.81253482) q[0];
rz(1.6883019) q[2];
sx q[2];
rz(-1.4521993) q[2];
sx q[2];
rz(0.30726739) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.62918) q[1];
sx q[1];
rz(-0.80958074) q[1];
sx q[1];
rz(-1.045524) q[1];
rz(0.047961162) q[3];
sx q[3];
rz(-0.75100079) q[3];
sx q[3];
rz(2.6329071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1314142) q[2];
sx q[2];
rz(-0.62128908) q[2];
sx q[2];
rz(-2.9739001) q[2];
rz(-0.0532648) q[3];
sx q[3];
rz(-2.0361418) q[3];
sx q[3];
rz(1.3767327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57529706) q[0];
sx q[0];
rz(-3.0360041) q[0];
sx q[0];
rz(-2.8422624) q[0];
rz(-0.37295595) q[1];
sx q[1];
rz(-1.2495709) q[1];
sx q[1];
rz(-1.4704871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84696124) q[0];
sx q[0];
rz(-1.9982013) q[0];
sx q[0];
rz(-0.60890108) q[0];
rz(-pi) q[1];
rz(1.6909356) q[2];
sx q[2];
rz(-1.1098301) q[2];
sx q[2];
rz(-1.1663988) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8052696) q[1];
sx q[1];
rz(-1.0674369) q[1];
sx q[1];
rz(0.663228) q[1];
x q[2];
rz(-1.3277131) q[3];
sx q[3];
rz(-0.63589261) q[3];
sx q[3];
rz(2.2896374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9386998) q[2];
sx q[2];
rz(-0.6627658) q[2];
sx q[2];
rz(1.1104442) q[2];
rz(2.2796196) q[3];
sx q[3];
rz(-1.5098666) q[3];
sx q[3];
rz(1.2601669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-0.92457572) q[0];
sx q[0];
rz(-2.45455) q[0];
sx q[0];
rz(-2.7365015) q[0];
rz(2.8054667) q[1];
sx q[1];
rz(-1.4210217) q[1];
sx q[1];
rz(2.1902671) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38430957) q[0];
sx q[0];
rz(-1.8258137) q[0];
sx q[0];
rz(0.41476212) q[0];
rz(-pi) q[1];
rz(1.3589922) q[2];
sx q[2];
rz(-0.84583218) q[2];
sx q[2];
rz(-2.447809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7860884) q[1];
sx q[1];
rz(-2.5230683) q[1];
sx q[1];
rz(-0.39822802) q[1];
rz(-2.5441465) q[3];
sx q[3];
rz(-2.8077112) q[3];
sx q[3];
rz(0.88874879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0508017) q[2];
sx q[2];
rz(-1.7369221) q[2];
sx q[2];
rz(-0.23078272) q[2];
rz(2.9595621) q[3];
sx q[3];
rz(-0.56806505) q[3];
sx q[3];
rz(0.52869421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35293216) q[0];
sx q[0];
rz(-3.1020628) q[0];
sx q[0];
rz(-0.93210644) q[0];
rz(0.63356361) q[1];
sx q[1];
rz(-2.0220058) q[1];
sx q[1];
rz(-1.8036141) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18620488) q[0];
sx q[0];
rz(-0.13372955) q[0];
sx q[0];
rz(0.9132847) q[0];
x q[1];
rz(-2.7510178) q[2];
sx q[2];
rz(-2.8468067) q[2];
sx q[2];
rz(-2.8975671) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7576068) q[1];
sx q[1];
rz(-1.9644871) q[1];
sx q[1];
rz(2.2298953) q[1];
rz(-pi) q[2];
rz(1.5389493) q[3];
sx q[3];
rz(-2.0689116) q[3];
sx q[3];
rz(-2.826988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46447095) q[2];
sx q[2];
rz(-1.2268343) q[2];
sx q[2];
rz(1.202549) q[2];
rz(0.36457148) q[3];
sx q[3];
rz(-0.69260827) q[3];
sx q[3];
rz(0.83474368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4778336) q[0];
sx q[0];
rz(-0.57357016) q[0];
sx q[0];
rz(-2.9649576) q[0];
rz(0.44003507) q[1];
sx q[1];
rz(-1.1853848) q[1];
sx q[1];
rz(2.1766591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10213908) q[0];
sx q[0];
rz(-0.23604611) q[0];
sx q[0];
rz(-0.59961598) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72146031) q[2];
sx q[2];
rz(-0.74154343) q[2];
sx q[2];
rz(0.78765819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.392004) q[1];
sx q[1];
rz(-0.67396213) q[1];
sx q[1];
rz(0.96629179) q[1];
x q[2];
rz(1.8550526) q[3];
sx q[3];
rz(-0.26460095) q[3];
sx q[3];
rz(-2.3456169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9289916) q[2];
sx q[2];
rz(-1.6182199) q[2];
sx q[2];
rz(-3.1316481) q[2];
rz(3.0250004) q[3];
sx q[3];
rz(-2.8023585) q[3];
sx q[3];
rz(-2.5998083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5803439) q[0];
sx q[0];
rz(-1.9191701) q[0];
sx q[0];
rz(2.5196581) q[0];
rz(-1.7033345) q[1];
sx q[1];
rz(-2.56918) q[1];
sx q[1];
rz(0.66351801) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2885482) q[0];
sx q[0];
rz(-2.2877734) q[0];
sx q[0];
rz(-1.8621481) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1691888) q[2];
sx q[2];
rz(-0.3065232) q[2];
sx q[2];
rz(2.1307132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92163699) q[1];
sx q[1];
rz(-1.8884648) q[1];
sx q[1];
rz(2.2457473) q[1];
x q[2];
rz(0.31511735) q[3];
sx q[3];
rz(-0.2956008) q[3];
sx q[3];
rz(-0.35741266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5180987) q[2];
sx q[2];
rz(-1.7524717) q[2];
sx q[2];
rz(-0.36250472) q[2];
rz(-2.6643961) q[3];
sx q[3];
rz(-2.0878744) q[3];
sx q[3];
rz(2.0188735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0839194) q[0];
sx q[0];
rz(-2.4102983) q[0];
sx q[0];
rz(0.90173632) q[0];
rz(-2.4665191) q[1];
sx q[1];
rz(-0.98552862) q[1];
sx q[1];
rz(-2.9973082) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0269723) q[0];
sx q[0];
rz(-1.2204224) q[0];
sx q[0];
rz(-3.0536985) q[0];
rz(-pi) q[1];
rz(-2.921943) q[2];
sx q[2];
rz(-1.6703419) q[2];
sx q[2];
rz(-1.8451898) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5965138) q[1];
sx q[1];
rz(-2.1987913) q[1];
sx q[1];
rz(1.2596115) q[1];
rz(0.4422632) q[3];
sx q[3];
rz(-1.398842) q[3];
sx q[3];
rz(0.32915914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6361864) q[2];
sx q[2];
rz(-1.212684) q[2];
sx q[2];
rz(0.20067659) q[2];
rz(1.1096654) q[3];
sx q[3];
rz(-2.6789013) q[3];
sx q[3];
rz(-0.5717352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5526445) q[0];
sx q[0];
rz(-2.1727967) q[0];
sx q[0];
rz(-1.1217242) q[0];
rz(2.6523392) q[1];
sx q[1];
rz(-1.6901292) q[1];
sx q[1];
rz(2.1314175) q[1];
rz(-1.7792173) q[2];
sx q[2];
rz(-0.42283146) q[2];
sx q[2];
rz(-2.6953807) q[2];
rz(2.6993467) q[3];
sx q[3];
rz(-1.5968235) q[3];
sx q[3];
rz(2.1376283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
