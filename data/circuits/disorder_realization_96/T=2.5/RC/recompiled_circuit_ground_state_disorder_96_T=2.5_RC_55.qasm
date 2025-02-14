OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.628217) q[0];
sx q[0];
rz(-1.4613419) q[0];
sx q[0];
rz(-2.2267377) q[0];
rz(2.272361) q[1];
sx q[1];
rz(2.4183122) q[1];
sx q[1];
rz(8.0719168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3769414) q[0];
sx q[0];
rz(-0.94263715) q[0];
sx q[0];
rz(0.92002042) q[0];
x q[1];
rz(2.3981744) q[2];
sx q[2];
rz(-1.2890491) q[2];
sx q[2];
rz(-0.27458686) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1754419) q[1];
sx q[1];
rz(-1.3805318) q[1];
sx q[1];
rz(-0.89118608) q[1];
rz(1.8736137) q[3];
sx q[3];
rz(-2.510364) q[3];
sx q[3];
rz(-2.6018179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.73504084) q[2];
sx q[2];
rz(-0.81321365) q[2];
sx q[2];
rz(-1.0068033) q[2];
rz(-1.6793647) q[3];
sx q[3];
rz(-1.4454602) q[3];
sx q[3];
rz(1.538895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2891069) q[0];
sx q[0];
rz(-2.2302581) q[0];
sx q[0];
rz(0.12595969) q[0];
rz(2.0055298) q[1];
sx q[1];
rz(-1.845153) q[1];
sx q[1];
rz(2.0819285) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.179716) q[0];
sx q[0];
rz(-1.5700514) q[0];
sx q[0];
rz(3.1396833) q[0];
rz(2.4074209) q[2];
sx q[2];
rz(-2.3245272) q[2];
sx q[2];
rz(1.2881607) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30970022) q[1];
sx q[1];
rz(-0.87617481) q[1];
sx q[1];
rz(-0.58160706) q[1];
rz(3.0094395) q[3];
sx q[3];
rz(-0.36214585) q[3];
sx q[3];
rz(-1.3782448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58553592) q[2];
sx q[2];
rz(-1.6126596) q[2];
sx q[2];
rz(-2.5118828) q[2];
rz(-1.2398047) q[3];
sx q[3];
rz(-2.3972378) q[3];
sx q[3];
rz(1.121421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-3.020551) q[0];
sx q[0];
rz(-2.910055) q[0];
sx q[0];
rz(-1.2303906) q[0];
rz(2.4379099) q[1];
sx q[1];
rz(-0.92862248) q[1];
sx q[1];
rz(-1.5812662) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11008292) q[0];
sx q[0];
rz(-0.39387273) q[0];
sx q[0];
rz(-1.6099694) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4034924) q[2];
sx q[2];
rz(-1.4261275) q[2];
sx q[2];
rz(1.6044782) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0496975) q[1];
sx q[1];
rz(-1.9822243) q[1];
sx q[1];
rz(2.0446214) q[1];
rz(1.1915474) q[3];
sx q[3];
rz(-0.82189188) q[3];
sx q[3];
rz(0.773202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33502093) q[2];
sx q[2];
rz(-2.0609914) q[2];
sx q[2];
rz(3.0873155) q[2];
rz(2.055376) q[3];
sx q[3];
rz(-2.351206) q[3];
sx q[3];
rz(2.6495892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739968) q[0];
sx q[0];
rz(-0.01570276) q[0];
sx q[0];
rz(0.98454654) q[0];
rz(1.4060075) q[1];
sx q[1];
rz(-1.5407298) q[1];
sx q[1];
rz(-2.9071009) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1403951) q[0];
sx q[0];
rz(-0.89522982) q[0];
sx q[0];
rz(2.7864149) q[0];
rz(-pi) q[1];
rz(-0.245204) q[2];
sx q[2];
rz(-1.2916995) q[2];
sx q[2];
rz(-0.14268219) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4602743) q[1];
sx q[1];
rz(-1.7806306) q[1];
sx q[1];
rz(0.55822452) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.029697604) q[3];
sx q[3];
rz(-2.1468024) q[3];
sx q[3];
rz(0.32054361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24388127) q[2];
sx q[2];
rz(-2.477024) q[2];
sx q[2];
rz(2.393874) q[2];
rz(2.054935) q[3];
sx q[3];
rz(-0.99435157) q[3];
sx q[3];
rz(0.35191107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20027593) q[0];
sx q[0];
rz(-0.87012297) q[0];
sx q[0];
rz(1.548832) q[0];
rz(3.0039655) q[1];
sx q[1];
rz(-0.96447861) q[1];
sx q[1];
rz(2.4639938) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1573059) q[0];
sx q[0];
rz(-1.56347) q[0];
sx q[0];
rz(-3.1276032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4196012) q[2];
sx q[2];
rz(-1.3791022) q[2];
sx q[2];
rz(2.5602532) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4487344) q[1];
sx q[1];
rz(-0.13693196) q[1];
sx q[1];
rz(-1.5278234) q[1];
rz(-pi) q[2];
rz(-1.8137535) q[3];
sx q[3];
rz(-0.84872171) q[3];
sx q[3];
rz(2.9606282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64348334) q[2];
sx q[2];
rz(-0.52916873) q[2];
sx q[2];
rz(0.37437487) q[2];
rz(0.74786413) q[3];
sx q[3];
rz(-1.6467843) q[3];
sx q[3];
rz(-1.098746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.13136524) q[0];
sx q[0];
rz(-1.7759198) q[0];
sx q[0];
rz(0.22450547) q[0];
rz(2.7911216) q[1];
sx q[1];
rz(-1.9572565) q[1];
sx q[1];
rz(1.1856461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3315546) q[0];
sx q[0];
rz(-1.7186949) q[0];
sx q[0];
rz(3.0857791) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19059981) q[2];
sx q[2];
rz(-2.6707044) q[2];
sx q[2];
rz(2.6550271) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95690475) q[1];
sx q[1];
rz(-0.84057759) q[1];
sx q[1];
rz(-2.5504677) q[1];
rz(-pi) q[2];
rz(2.3697183) q[3];
sx q[3];
rz(-0.83362416) q[3];
sx q[3];
rz(1.9825359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4993569) q[2];
sx q[2];
rz(-2.2580052) q[2];
sx q[2];
rz(2.9767766) q[2];
rz(2.3757101) q[3];
sx q[3];
rz(-0.82847786) q[3];
sx q[3];
rz(2.5501521) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0484475) q[0];
sx q[0];
rz(-0.026938139) q[0];
sx q[0];
rz(0.41937605) q[0];
rz(1.558149) q[1];
sx q[1];
rz(-0.42834586) q[1];
sx q[1];
rz(1.8289808) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2700513) q[0];
sx q[0];
rz(-1.7333374) q[0];
sx q[0];
rz(2.1664361) q[0];
rz(-pi) q[1];
rz(0.68676853) q[2];
sx q[2];
rz(-2.1109952) q[2];
sx q[2];
rz(1.6322002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.928682) q[1];
sx q[1];
rz(-0.31344673) q[1];
sx q[1];
rz(2.2261966) q[1];
rz(-pi) q[2];
rz(0.26793865) q[3];
sx q[3];
rz(-1.6959785) q[3];
sx q[3];
rz(0.94733688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7340362) q[2];
sx q[2];
rz(-2.4702256) q[2];
sx q[2];
rz(-0.35801312) q[2];
rz(0.19896209) q[3];
sx q[3];
rz(-2.3301061) q[3];
sx q[3];
rz(-3.0437886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1426549) q[0];
sx q[0];
rz(-2.1079347) q[0];
sx q[0];
rz(-0.79310966) q[0];
rz(-2.2456887) q[1];
sx q[1];
rz(-1.9205807) q[1];
sx q[1];
rz(-2.6633247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60837854) q[0];
sx q[0];
rz(-1.694931) q[0];
sx q[0];
rz(-0.45093765) q[0];
rz(-pi) q[1];
rz(1.9321279) q[2];
sx q[2];
rz(-1.1602957) q[2];
sx q[2];
rz(0.77958661) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7859808) q[1];
sx q[1];
rz(-2.0352239) q[1];
sx q[1];
rz(2.8465438) q[1];
x q[2];
rz(1.7630799) q[3];
sx q[3];
rz(-2.6182976) q[3];
sx q[3];
rz(1.3973325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6200977) q[2];
sx q[2];
rz(-1.564881) q[2];
sx q[2];
rz(-0.0077956789) q[2];
rz(1.1826285) q[3];
sx q[3];
rz(-1.7595485) q[3];
sx q[3];
rz(-1.8462605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1496534) q[0];
sx q[0];
rz(-2.6872334) q[0];
sx q[0];
rz(0.40865189) q[0];
rz(-1.0952134) q[1];
sx q[1];
rz(-1.4827671) q[1];
sx q[1];
rz(0.85831395) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3811144) q[0];
sx q[0];
rz(-1.71813) q[0];
sx q[0];
rz(3.0744746) q[0];
rz(-pi) q[1];
rz(3.0093391) q[2];
sx q[2];
rz(-2.0200854) q[2];
sx q[2];
rz(1.1960891) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62095642) q[1];
sx q[1];
rz(-0.22527105) q[1];
sx q[1];
rz(-2.0849289) q[1];
rz(-pi) q[2];
rz(2.5439203) q[3];
sx q[3];
rz(-1.7558756) q[3];
sx q[3];
rz(1.5782034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7909164) q[2];
sx q[2];
rz(-0.75296593) q[2];
sx q[2];
rz(2.47993) q[2];
rz(-0.800313) q[3];
sx q[3];
rz(-1.1626264) q[3];
sx q[3];
rz(-0.77320981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3043587) q[0];
sx q[0];
rz(-2.9999314) q[0];
sx q[0];
rz(0.71453553) q[0];
rz(-0.47572687) q[1];
sx q[1];
rz(-2.5560502) q[1];
sx q[1];
rz(-2.661396) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5479832) q[0];
sx q[0];
rz(-2.3385323) q[0];
sx q[0];
rz(0.70269967) q[0];
rz(-pi) q[1];
rz(-2.8101361) q[2];
sx q[2];
rz(-0.94539795) q[2];
sx q[2];
rz(-1.7632222) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3298397) q[1];
sx q[1];
rz(-1.2851213) q[1];
sx q[1];
rz(-2.4655254) q[1];
rz(-pi) q[2];
rz(-1.1350182) q[3];
sx q[3];
rz(-1.3131071) q[3];
sx q[3];
rz(0.99761744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6314038) q[2];
sx q[2];
rz(-1.9437342) q[2];
sx q[2];
rz(-2.8655444) q[2];
rz(2.7035642) q[3];
sx q[3];
rz(-1.5166401) q[3];
sx q[3];
rz(-0.62694222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.940687) q[0];
sx q[0];
rz(-1.154366) q[0];
sx q[0];
rz(0.05703297) q[0];
rz(0.036046473) q[1];
sx q[1];
rz(-1.5279952) q[1];
sx q[1];
rz(-1.5855018) q[1];
rz(-1.0895928) q[2];
sx q[2];
rz(-0.54554987) q[2];
sx q[2];
rz(0.050532374) q[2];
rz(1.7071758) q[3];
sx q[3];
rz(-1.3188667) q[3];
sx q[3];
rz(1.8018166) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
