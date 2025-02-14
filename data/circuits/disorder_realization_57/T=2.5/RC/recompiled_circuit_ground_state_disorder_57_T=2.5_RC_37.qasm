OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8785414) q[0];
sx q[0];
rz(-0.6113373) q[0];
sx q[0];
rz(-1.6482469) q[0];
rz(-2.2892294) q[1];
sx q[1];
rz(-1.747793) q[1];
sx q[1];
rz(1.2300904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85948259) q[0];
sx q[0];
rz(-1.9766269) q[0];
sx q[0];
rz(2.6187483) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1793433) q[2];
sx q[2];
rz(-1.9916996) q[2];
sx q[2];
rz(-0.98513033) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7609357) q[1];
sx q[1];
rz(-2.0208997) q[1];
sx q[1];
rz(0.37827079) q[1];
x q[2];
rz(1.3268747) q[3];
sx q[3];
rz(-0.29196843) q[3];
sx q[3];
rz(1.0002182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.40296951) q[2];
sx q[2];
rz(-1.7454742) q[2];
sx q[2];
rz(1.1633066) q[2];
rz(-2.2033384) q[3];
sx q[3];
rz(-0.80621755) q[3];
sx q[3];
rz(1.4279385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52629483) q[0];
sx q[0];
rz(-1.9087003) q[0];
sx q[0];
rz(-1.9799318) q[0];
rz(1.42234) q[1];
sx q[1];
rz(-1.4931581) q[1];
sx q[1];
rz(0.086440451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7122262) q[0];
sx q[0];
rz(-1.7284358) q[0];
sx q[0];
rz(-1.0382354) q[0];
rz(-0.49098726) q[2];
sx q[2];
rz(-1.3730959) q[2];
sx q[2];
rz(-0.035015496) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.045043) q[1];
sx q[1];
rz(-0.031686671) q[1];
sx q[1];
rz(2.6331054) q[1];
x q[2];
rz(-1.3996313) q[3];
sx q[3];
rz(-1.6361208) q[3];
sx q[3];
rz(2.8166376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.57261673) q[2];
sx q[2];
rz(-2.6628351) q[2];
sx q[2];
rz(-3.0968481) q[2];
rz(-2.4217126) q[3];
sx q[3];
rz(-1.5174815) q[3];
sx q[3];
rz(1.6911136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2429263) q[0];
sx q[0];
rz(-1.1761605) q[0];
sx q[0];
rz(-1.9157008) q[0];
rz(0.9179081) q[1];
sx q[1];
rz(-1.8501015) q[1];
sx q[1];
rz(-1.8570522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63375134) q[0];
sx q[0];
rz(-1.0695262) q[0];
sx q[0];
rz(-3.0847163) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.011083885) q[2];
sx q[2];
rz(-0.62168089) q[2];
sx q[2];
rz(1.7116837) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1333993) q[1];
sx q[1];
rz(-2.2643406) q[1];
sx q[1];
rz(0.72024017) q[1];
rz(-1.8281047) q[3];
sx q[3];
rz(-0.9633102) q[3];
sx q[3];
rz(-1.2953386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43117493) q[2];
sx q[2];
rz(-1.147889) q[2];
sx q[2];
rz(-0.61719027) q[2];
rz(0.91444531) q[3];
sx q[3];
rz(-2.4121598) q[3];
sx q[3];
rz(0.35703737) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0964088) q[0];
sx q[0];
rz(-3.0953188) q[0];
sx q[0];
rz(0.98840493) q[0];
rz(1.660396) q[1];
sx q[1];
rz(-2.587187) q[1];
sx q[1];
rz(-0.60715094) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2961126) q[0];
sx q[0];
rz(-0.37808266) q[0];
sx q[0];
rz(-2.3251359) q[0];
x q[1];
rz(0.49272196) q[2];
sx q[2];
rz(-0.70198503) q[2];
sx q[2];
rz(1.5240325) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7242291) q[1];
sx q[1];
rz(-2.0290142) q[1];
sx q[1];
rz(2.7643137) q[1];
rz(-pi) q[2];
rz(-0.097154599) q[3];
sx q[3];
rz(-1.7504255) q[3];
sx q[3];
rz(3.0727974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0561169) q[2];
sx q[2];
rz(-1.3885219) q[2];
sx q[2];
rz(-2.910854) q[2];
rz(-2.3750316) q[3];
sx q[3];
rz(-0.14403382) q[3];
sx q[3];
rz(-1.2994331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6799927) q[0];
sx q[0];
rz(-1.5138641) q[0];
sx q[0];
rz(-0.067721279) q[0];
rz(-1.8374775) q[1];
sx q[1];
rz(-1.8505406) q[1];
sx q[1];
rz(-1.4363323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50698603) q[0];
sx q[0];
rz(-1.9546967) q[0];
sx q[0];
rz(0.83457729) q[0];
rz(-1.0772938) q[2];
sx q[2];
rz(-1.4736325) q[2];
sx q[2];
rz(0.0070126931) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.41714121) q[1];
sx q[1];
rz(-2.0421237) q[1];
sx q[1];
rz(0.36496867) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0338289) q[3];
sx q[3];
rz(-2.7269533) q[3];
sx q[3];
rz(1.6920167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0453673) q[2];
sx q[2];
rz(-2.2787091) q[2];
sx q[2];
rz(-2.6048342) q[2];
rz(1.9219575) q[3];
sx q[3];
rz(-0.9044956) q[3];
sx q[3];
rz(1.1902683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9959975) q[0];
sx q[0];
rz(-0.97766101) q[0];
sx q[0];
rz(-1.4033432) q[0];
rz(2.6373236) q[1];
sx q[1];
rz(-1.206617) q[1];
sx q[1];
rz(0.23955841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5546416) q[0];
sx q[0];
rz(-1.3901781) q[0];
sx q[0];
rz(0.15968542) q[0];
rz(-pi) q[1];
rz(-0.94953612) q[2];
sx q[2];
rz(-2.4944948) q[2];
sx q[2];
rz(2.4280649) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9592683) q[1];
sx q[1];
rz(-2.2367918) q[1];
sx q[1];
rz(1.6823425) q[1];
rz(0.91383818) q[3];
sx q[3];
rz(-2.5322134) q[3];
sx q[3];
rz(2.1003124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2704894) q[2];
sx q[2];
rz(-1.6051382) q[2];
sx q[2];
rz(-0.25320539) q[2];
rz(-2.1829055) q[3];
sx q[3];
rz(-0.91512338) q[3];
sx q[3];
rz(-1.8668713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7764928) q[0];
sx q[0];
rz(-0.72117844) q[0];
sx q[0];
rz(-1.391063) q[0];
rz(-2.5513388) q[1];
sx q[1];
rz(-1.2140467) q[1];
sx q[1];
rz(-2.6388157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4462379) q[0];
sx q[0];
rz(-0.10646755) q[0];
sx q[0];
rz(-1.4346532) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4847772) q[2];
sx q[2];
rz(-1.2716827) q[2];
sx q[2];
rz(2.6292036) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5107837) q[1];
sx q[1];
rz(-1.3250588) q[1];
sx q[1];
rz(2.3151957) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2572978) q[3];
sx q[3];
rz(-1.0491706) q[3];
sx q[3];
rz(2.9056463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7949924) q[2];
sx q[2];
rz(-1.9575926) q[2];
sx q[2];
rz(-2.6896175) q[2];
rz(0.31451264) q[3];
sx q[3];
rz(-1.9032685) q[3];
sx q[3];
rz(0.048129169) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6314342) q[0];
sx q[0];
rz(-0.87961125) q[0];
sx q[0];
rz(-1.6360224) q[0];
rz(3.0360041) q[1];
sx q[1];
rz(-1.0786723) q[1];
sx q[1];
rz(0.67965913) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1435232) q[0];
sx q[0];
rz(-0.98833246) q[0];
sx q[0];
rz(-2.4650455) q[0];
x q[1];
rz(-0.40186581) q[2];
sx q[2];
rz(-2.1500714) q[2];
sx q[2];
rz(-1.4039672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1464771) q[1];
sx q[1];
rz(-1.5676284) q[1];
sx q[1];
rz(2.0407487) q[1];
rz(-pi) q[2];
rz(0.18631492) q[3];
sx q[3];
rz(-1.3491231) q[3];
sx q[3];
rz(-1.1849537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52123657) q[2];
sx q[2];
rz(-0.53400365) q[2];
sx q[2];
rz(2.3109069) q[2];
rz(2.1031117) q[3];
sx q[3];
rz(-0.65282789) q[3];
sx q[3];
rz(-1.7415107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.1164383) q[0];
sx q[0];
rz(-1.9898299) q[0];
sx q[0];
rz(-2.108216) q[0];
rz(2.0694464) q[1];
sx q[1];
rz(-1.9629581) q[1];
sx q[1];
rz(2.5796366) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852405) q[0];
sx q[0];
rz(-0.7228719) q[0];
sx q[0];
rz(1.7815352) q[0];
x q[1];
rz(0.35530151) q[2];
sx q[2];
rz(-0.94569478) q[2];
sx q[2];
rz(1.6297766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.452637) q[1];
sx q[1];
rz(-1.9105043) q[1];
sx q[1];
rz(0.79381408) q[1];
rz(-pi) q[2];
rz(-1.6283243) q[3];
sx q[3];
rz(-1.5794601) q[3];
sx q[3];
rz(1.5271562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0434887) q[2];
sx q[2];
rz(-1.4926814) q[2];
sx q[2];
rz(-1.4081504) q[2];
rz(0.060338542) q[3];
sx q[3];
rz(-1.8294168) q[3];
sx q[3];
rz(1.8797125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9982346) q[0];
sx q[0];
rz(-0.77333212) q[0];
sx q[0];
rz(-1.3704569) q[0];
rz(-1.0073608) q[1];
sx q[1];
rz(-1.3887364) q[1];
sx q[1];
rz(-0.14152424) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62011224) q[0];
sx q[0];
rz(-0.65516383) q[0];
sx q[0];
rz(-2.377009) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.510364) q[2];
sx q[2];
rz(-1.6360893) q[2];
sx q[2];
rz(-3.1217965) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2685769) q[1];
sx q[1];
rz(-2.0076224) q[1];
sx q[1];
rz(-0.42307968) q[1];
rz(-pi) q[2];
rz(0.022079682) q[3];
sx q[3];
rz(-1.0754657) q[3];
sx q[3];
rz(-1.2091523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0274028) q[2];
sx q[2];
rz(-1.1833311) q[2];
sx q[2];
rz(-1.8751422) q[2];
rz(0.13628422) q[3];
sx q[3];
rz(-2.2253021) q[3];
sx q[3];
rz(2.9436679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4895353) q[0];
sx q[0];
rz(-1.0132402) q[0];
sx q[0];
rz(-1.5860438) q[0];
rz(-2.1736705) q[1];
sx q[1];
rz(-2.617901) q[1];
sx q[1];
rz(-1.0760457) q[1];
rz(-0.10836149) q[2];
sx q[2];
rz(-0.81365267) q[2];
sx q[2];
rz(-1.1968422) q[2];
rz(1.7607197) q[3];
sx q[3];
rz(-1.6217762) q[3];
sx q[3];
rz(-3.1358596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
