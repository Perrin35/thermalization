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
rz(0.99474466) q[0];
sx q[0];
rz(-2.5543307) q[0];
sx q[0];
rz(0.62556148) q[0];
rz(-2.5211531) q[1];
sx q[1];
rz(-0.99045366) q[1];
sx q[1];
rz(-0.78392309) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7083354) q[0];
sx q[0];
rz(-1.1366708) q[0];
sx q[0];
rz(-1.5668847) q[0];
x q[1];
rz(1.6541846) q[2];
sx q[2];
rz(-1.4321799) q[2];
sx q[2];
rz(2.1237789) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.82134889) q[1];
sx q[1];
rz(-0.87961266) q[1];
sx q[1];
rz(0.76637474) q[1];
rz(-pi) q[2];
rz(0.77772452) q[3];
sx q[3];
rz(-1.6334264) q[3];
sx q[3];
rz(-0.81111139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1938532) q[2];
sx q[2];
rz(-1.0449907) q[2];
sx q[2];
rz(-0.6673153) q[2];
rz(-0.31467485) q[3];
sx q[3];
rz(-2.5421725) q[3];
sx q[3];
rz(0.92710322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84745234) q[0];
sx q[0];
rz(-2.2951234) q[0];
sx q[0];
rz(-2.8417929) q[0];
rz(1.0046129) q[1];
sx q[1];
rz(-2.424898) q[1];
sx q[1];
rz(0.85003781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77281053) q[0];
sx q[0];
rz(-1.6666935) q[0];
sx q[0];
rz(0.96536388) q[0];
rz(-2.7072435) q[2];
sx q[2];
rz(-2.953989) q[2];
sx q[2];
rz(-1.4570731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0850573) q[1];
sx q[1];
rz(-1.3621482) q[1];
sx q[1];
rz(1.2056808) q[1];
x q[2];
rz(-2.7369954) q[3];
sx q[3];
rz(-1.5581774) q[3];
sx q[3];
rz(-0.22443988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5689508) q[2];
sx q[2];
rz(-2.6093542) q[2];
sx q[2];
rz(2.4786095) q[2];
rz(-0.7524544) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(-0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982584) q[0];
sx q[0];
rz(-0.4011811) q[0];
sx q[0];
rz(0.67984003) q[0];
rz(-0.72194779) q[1];
sx q[1];
rz(-2.5847021) q[1];
sx q[1];
rz(-2.360875) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76533045) q[0];
sx q[0];
rz(-1.1701382) q[0];
sx q[0];
rz(-3.1332934) q[0];
rz(1.4846615) q[2];
sx q[2];
rz(-1.5745828) q[2];
sx q[2];
rz(2.6876793) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7739984) q[1];
sx q[1];
rz(-1.9450608) q[1];
sx q[1];
rz(0.67919977) q[1];
rz(-pi) q[2];
rz(2.3060477) q[3];
sx q[3];
rz(-1.7342596) q[3];
sx q[3];
rz(-0.36629656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22689936) q[2];
sx q[2];
rz(-1.2279899) q[2];
sx q[2];
rz(0.66960382) q[2];
rz(-2.2412444) q[3];
sx q[3];
rz(-2.1293631) q[3];
sx q[3];
rz(2.3646234) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8825619) q[0];
sx q[0];
rz(-0.43864033) q[0];
sx q[0];
rz(0.90748179) q[0];
rz(-0.59440815) q[1];
sx q[1];
rz(-2.2120357) q[1];
sx q[1];
rz(2.0118735) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27760273) q[0];
sx q[0];
rz(-2.4765795) q[0];
sx q[0];
rz(-1.8825085) q[0];
rz(1.8215034) q[2];
sx q[2];
rz(-2.5959229) q[2];
sx q[2];
rz(2.9755993) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4817244) q[1];
sx q[1];
rz(-0.65978434) q[1];
sx q[1];
rz(-2.6901812) q[1];
rz(2.3334585) q[3];
sx q[3];
rz(-0.39565797) q[3];
sx q[3];
rz(1.7803223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10924673) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(-2.4692811) q[2];
rz(3.1224871) q[3];
sx q[3];
rz(-2.8002383) q[3];
sx q[3];
rz(0.13489558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371209) q[0];
sx q[0];
rz(-2.3470375) q[0];
sx q[0];
rz(0.16803148) q[0];
rz(-1.2270323) q[1];
sx q[1];
rz(-1.2201759) q[1];
sx q[1];
rz(-2.8438445) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.346031) q[0];
sx q[0];
rz(-2.3712284) q[0];
sx q[0];
rz(2.3576231) q[0];
x q[1];
rz(-0.99397583) q[2];
sx q[2];
rz(-2.0703531) q[2];
sx q[2];
rz(0.039488878) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6444466) q[1];
sx q[1];
rz(-0.20731099) q[1];
sx q[1];
rz(-3.0203692) q[1];
rz(0.00068126531) q[3];
sx q[3];
rz(-1.7295803) q[3];
sx q[3];
rz(-0.53631594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9510522) q[2];
sx q[2];
rz(-2.0420044) q[2];
sx q[2];
rz(2.5229048) q[2];
rz(2.3330073) q[3];
sx q[3];
rz(-2.6638668) q[3];
sx q[3];
rz(-1.8019069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.488778) q[0];
sx q[0];
rz(-1.6365016) q[0];
sx q[0];
rz(-2.6763647) q[0];
rz(-0.48509994) q[1];
sx q[1];
rz(-2.7310889) q[1];
sx q[1];
rz(3.0533275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34644914) q[0];
sx q[0];
rz(-2.4426021) q[0];
sx q[0];
rz(-2.3153852) q[0];
x q[1];
rz(0.012862269) q[2];
sx q[2];
rz(-0.84190997) q[2];
sx q[2];
rz(-0.60573214) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40988556) q[1];
sx q[1];
rz(-0.8250069) q[1];
sx q[1];
rz(-1.3456001) q[1];
rz(1.4824105) q[3];
sx q[3];
rz(-1.2840349) q[3];
sx q[3];
rz(2.4915316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39468592) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(0.49760231) q[2];
rz(0.48007128) q[3];
sx q[3];
rz(-0.92104715) q[3];
sx q[3];
rz(-3.0729496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.8467872) q[0];
sx q[0];
rz(-0.46600309) q[0];
sx q[0];
rz(-1.71126) q[0];
rz(-1.3149186) q[1];
sx q[1];
rz(-0.39350915) q[1];
sx q[1];
rz(1.5550782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67051149) q[0];
sx q[0];
rz(-1.9626612) q[0];
sx q[0];
rz(3.0125822) q[0];
x q[1];
rz(-2.2261593) q[2];
sx q[2];
rz(-0.829773) q[2];
sx q[2];
rz(-2.1754153) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58007121) q[1];
sx q[1];
rz(-2.5513756) q[1];
sx q[1];
rz(1.3886222) q[1];
rz(1.4659381) q[3];
sx q[3];
rz(-0.63408454) q[3];
sx q[3];
rz(0.88632562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84307182) q[2];
sx q[2];
rz(-0.71121794) q[2];
sx q[2];
rz(2.2141875) q[2];
rz(-1.1585506) q[3];
sx q[3];
rz(-0.68803334) q[3];
sx q[3];
rz(-3.0454175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7189099) q[0];
sx q[0];
rz(-2.2490608) q[0];
sx q[0];
rz(-0.476015) q[0];
rz(2.1503275) q[1];
sx q[1];
rz(-2.3920993) q[1];
sx q[1];
rz(-3.057726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1779685) q[0];
sx q[0];
rz(-1.9289079) q[0];
sx q[0];
rz(2.0923418) q[0];
rz(2.6798579) q[2];
sx q[2];
rz(-1.8991915) q[2];
sx q[2];
rz(2.6685614) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7788811) q[1];
sx q[1];
rz(-1.485752) q[1];
sx q[1];
rz(1.9910732) q[1];
rz(3.0552577) q[3];
sx q[3];
rz(-1.207296) q[3];
sx q[3];
rz(-2.8192161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91741651) q[2];
sx q[2];
rz(-2.8755964) q[2];
sx q[2];
rz(0.64845294) q[2];
rz(-2.2367541) q[3];
sx q[3];
rz(-1.6831393) q[3];
sx q[3];
rz(2.8308433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34173486) q[0];
sx q[0];
rz(-2.4315727) q[0];
sx q[0];
rz(2.572686) q[0];
rz(-2.3078602) q[1];
sx q[1];
rz(-1.2833475) q[1];
sx q[1];
rz(-1.0047097) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53225258) q[0];
sx q[0];
rz(-2.1933317) q[0];
sx q[0];
rz(-1.3329766) q[0];
rz(-pi) q[1];
rz(2.9113681) q[2];
sx q[2];
rz(-2.0966999) q[2];
sx q[2];
rz(2.3221225) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8160287) q[1];
sx q[1];
rz(-0.70239151) q[1];
sx q[1];
rz(-1.4588474) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12281872) q[3];
sx q[3];
rz(-2.5505722) q[3];
sx q[3];
rz(1.0367048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4363165) q[2];
sx q[2];
rz(-2.269561) q[2];
sx q[2];
rz(2.7936068) q[2];
rz(-0.82585382) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(2.5364449) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0572877) q[0];
sx q[0];
rz(-1.0658406) q[0];
sx q[0];
rz(-2.9344946) q[0];
rz(2.6026978) q[1];
sx q[1];
rz(-2.5205044) q[1];
sx q[1];
rz(2.5394687) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.759256) q[0];
sx q[0];
rz(-0.70092541) q[0];
sx q[0];
rz(0.78030326) q[0];
rz(-pi) q[1];
rz(-2.5670577) q[2];
sx q[2];
rz(-0.67107397) q[2];
sx q[2];
rz(0.35071638) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94127266) q[1];
sx q[1];
rz(-1.4956022) q[1];
sx q[1];
rz(-1.978051) q[1];
rz(-2.3942457) q[3];
sx q[3];
rz(-1.8149733) q[3];
sx q[3];
rz(-1.3507193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3567317) q[2];
sx q[2];
rz(-2.1897903) q[2];
sx q[2];
rz(-1.7968563) q[2];
rz(-2.6209659) q[3];
sx q[3];
rz(-2.8713363) q[3];
sx q[3];
rz(2.3394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44705924) q[0];
sx q[0];
rz(-0.72493989) q[0];
sx q[0];
rz(-1.3865393) q[0];
rz(-0.78320349) q[1];
sx q[1];
rz(-1.7192817) q[1];
sx q[1];
rz(-1.5789938) q[1];
rz(-0.34306768) q[2];
sx q[2];
rz(-1.82556) q[2];
sx q[2];
rz(-1.0221046) q[2];
rz(-0.047719638) q[3];
sx q[3];
rz(-0.24038355) q[3];
sx q[3];
rz(3.0265831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
