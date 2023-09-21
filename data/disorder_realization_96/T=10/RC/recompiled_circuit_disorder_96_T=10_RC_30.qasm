OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(-2.3242216) q[0];
rz(0.983239) q[1];
sx q[1];
rz(-0.53951889) q[1];
sx q[1];
rz(1.9411545) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96785986) q[0];
sx q[0];
rz(-1.1636486) q[0];
sx q[0];
rz(1.3074387) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3517411) q[2];
sx q[2];
rz(-0.54974216) q[2];
sx q[2];
rz(1.4866536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9125036) q[1];
sx q[1];
rz(-0.87516811) q[1];
sx q[1];
rz(1.2234664) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8957542) q[3];
sx q[3];
rz(-0.80245362) q[3];
sx q[3];
rz(1.5047531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1203221) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(1.5585287) q[2];
rz(-0.99672404) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(-0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(-3.0875207) q[0];
rz(1.9460829) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(2.6057459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5553404) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(-2.9931086) q[0];
x q[1];
rz(2.8500154) q[2];
sx q[2];
rz(-0.73086408) q[2];
sx q[2];
rz(-2.4183194) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0332196) q[1];
sx q[1];
rz(-0.48711005) q[1];
sx q[1];
rz(2.5750722) q[1];
rz(-pi) q[2];
rz(-1.4278533) q[3];
sx q[3];
rz(-1.6030451) q[3];
sx q[3];
rz(-1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1198931) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(2.1832441) q[2];
rz(3.0751394) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(2.691793) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41985837) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(-0.15047519) q[0];
rz(-0.45723215) q[1];
sx q[1];
rz(-2.229264) q[1];
sx q[1];
rz(0.025807468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8919825) q[0];
sx q[0];
rz(-1.5907856) q[0];
sx q[0];
rz(1.32094) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0518968) q[2];
sx q[2];
rz(-1.8811474) q[2];
sx q[2];
rz(0.87583625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3734696) q[1];
sx q[1];
rz(-0.79263955) q[1];
sx q[1];
rz(2.5440689) q[1];
rz(-3.0098626) q[3];
sx q[3];
rz(-1.1909435) q[3];
sx q[3];
rz(-0.47117885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0187443) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(2.9600926) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.240775) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(2.9649819) q[0];
rz(0.88090849) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(0.53612971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4767544) q[0];
sx q[0];
rz(-1.9739445) q[0];
sx q[0];
rz(0.82157764) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5136112) q[2];
sx q[2];
rz(-1.3344889) q[2];
sx q[2];
rz(-1.1255217) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22524658) q[1];
sx q[1];
rz(-2.1919474) q[1];
sx q[1];
rz(1.1017373) q[1];
rz(1.2767775) q[3];
sx q[3];
rz(-0.91563581) q[3];
sx q[3];
rz(2.7340207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6716016) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(1.0446576) q[2];
rz(2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(-2.7846591) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9064643) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(-1.0513603) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(3.0984745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2152104) q[0];
sx q[0];
rz(-2.1676817) q[0];
sx q[0];
rz(0.2290639) q[0];
x q[1];
rz(-2.0448858) q[2];
sx q[2];
rz(-0.57069639) q[2];
sx q[2];
rz(-1.8622461) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4823618) q[1];
sx q[1];
rz(-1.8419957) q[1];
sx q[1];
rz(1.8455452) q[1];
rz(-pi) q[2];
rz(-1.7205914) q[3];
sx q[3];
rz(-2.3547223) q[3];
sx q[3];
rz(-2.8926135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0126426) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(0.35219231) q[2];
rz(-0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.5181638) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(0.75025264) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(1.0587143) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1832702) q[0];
sx q[0];
rz(-2.000862) q[0];
sx q[0];
rz(1.3413315) q[0];
rz(1.9361587) q[2];
sx q[2];
rz(-2.0632671) q[2];
sx q[2];
rz(-1.1577215) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8844879) q[1];
sx q[1];
rz(-0.37933644) q[1];
sx q[1];
rz(2.2124955) q[1];
rz(-pi) q[2];
rz(1.3771463) q[3];
sx q[3];
rz(-2.2691257) q[3];
sx q[3];
rz(-0.26102548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6713312) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.2188101) q[2];
rz(-1.1550711) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(-1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0999775) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.51145) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-0.310985) q[1];
sx q[1];
rz(-0.84164936) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0253042) q[0];
sx q[0];
rz(-1.7968654) q[0];
sx q[0];
rz(-2.6095005) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26288962) q[2];
sx q[2];
rz(-2.4979696) q[2];
sx q[2];
rz(-2.7123244) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51534286) q[1];
sx q[1];
rz(-1.6193577) q[1];
sx q[1];
rz(1.6471144) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9229786) q[3];
sx q[3];
rz(-2.300005) q[3];
sx q[3];
rz(-1.5515755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1402309) q[2];
sx q[2];
rz(-1.7732239) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(2.4800381) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(-0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89896232) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(-1.7393973) q[0];
rz(-0.095480355) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(0.41762525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71224371) q[0];
sx q[0];
rz(-1.1375543) q[0];
sx q[0];
rz(0.54746763) q[0];
rz(-pi) q[1];
rz(3.0268961) q[2];
sx q[2];
rz(-1.2071929) q[2];
sx q[2];
rz(0.39322688) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.8199181) q[1];
sx q[1];
rz(-2.3707317) q[1];
sx q[1];
rz(2.2714771) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2291064) q[3];
sx q[3];
rz(-2.6819326) q[3];
sx q[3];
rz(1.0115136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79545704) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(-2.0987089) q[2];
rz(2.4677094) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8326571) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(-0.62966627) q[0];
rz(-0.57485238) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(0.94690698) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6553584) q[0];
sx q[0];
rz(-1.8706733) q[0];
sx q[0];
rz(-0.34098682) q[0];
x q[1];
rz(-0.38893716) q[2];
sx q[2];
rz(-3.0367594) q[2];
sx q[2];
rz(-2.7742085) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.088989181) q[1];
sx q[1];
rz(-0.34480428) q[1];
sx q[1];
rz(2.8903972) q[1];
rz(-1.1287273) q[3];
sx q[3];
rz(-1.3978492) q[3];
sx q[3];
rz(-0.63463075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(1.2488731) q[2];
rz(-2.4272264) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0749851) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(-2.4865436) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(2.4972829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9172168) q[0];
sx q[0];
rz(-1.7443568) q[0];
sx q[0];
rz(-1.6266842) q[0];
x q[1];
rz(2.9115453) q[2];
sx q[2];
rz(-2.2736079) q[2];
sx q[2];
rz(2.7237797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88313738) q[1];
sx q[1];
rz(-1.2346134) q[1];
sx q[1];
rz(1.2698445) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44390042) q[3];
sx q[3];
rz(-2.0022941) q[3];
sx q[3];
rz(2.3865226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3506938) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(0.59990668) q[2];
rz(-0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(-1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29466378) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-0.22944336) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(-2.1951998) q[2];
sx q[2];
rz(-1.4675491) q[2];
sx q[2];
rz(1.8617392) q[2];
rz(1.0944081) q[3];
sx q[3];
rz(-0.73275685) q[3];
sx q[3];
rz(-2.3263596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];