OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.03798401355743) q[0];
sx q[0];
rz(4.82595661480958) q[0];
sx q[0];
rz(11.2136002540509) q[0];
rz(-5.09851312637329) q[1];
sx q[1];
rz(2.46896222432191) q[1];
sx q[1];
rz(11.0506347179334) q[1];
cx q[1],q[0];
rz(0.112134508788586) q[0];
sx q[0];
rz(1.37282756169374) q[0];
sx q[0];
rz(10.8056934833448) q[0];
rz(-2.04987287521362) q[2];
sx q[2];
rz(7.14833894570405) q[2];
sx q[2];
rz(2.39126632212802) q[2];
cx q[2],q[1];
rz(1.37990725040436) q[1];
sx q[1];
rz(5.68220058281953) q[1];
sx q[1];
rz(11.6488220453183) q[1];
rz(9.71189022064209) q[3];
sx q[3];
rz(1.9762044270807) q[3];
sx q[3];
rz(7.73041985034152) q[3];
cx q[3],q[2];
rz(4.1451940536499) q[2];
sx q[2];
rz(3.03339520295198) q[2];
sx q[2];
rz(12.3803744077603) q[2];
rz(-3.16023182868958) q[3];
sx q[3];
rz(7.57747093041474) q[3];
sx q[3];
rz(11.4960241079251) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.58285343647003) q[0];
sx q[0];
rz(5.71696558793122) q[0];
sx q[0];
rz(9.10641685723468) q[0];
rz(3.77862763404846) q[1];
sx q[1];
rz(5.50326267083222) q[1];
sx q[1];
rz(6.75242612361118) q[1];
cx q[1],q[0];
rz(5.31175804138184) q[0];
sx q[0];
rz(4.88323965867097) q[0];
sx q[0];
rz(8.05592606066867) q[0];
rz(-1.65474569797516) q[2];
sx q[2];
rz(4.33779660065705) q[2];
sx q[2];
rz(8.52590993642017) q[2];
cx q[2],q[1];
rz(-1.79828953742981) q[1];
sx q[1];
rz(5.03158429463441) q[1];
sx q[1];
rz(12.3273429632108) q[1];
rz(-0.00692650070413947) q[3];
sx q[3];
rz(0.68994656403596) q[3];
sx q[3];
rz(10.9258488178174) q[3];
cx q[3],q[2];
rz(1.31767535209656) q[2];
sx q[2];
rz(3.35183413525159) q[2];
sx q[2];
rz(13.2596447229306) q[2];
rz(-1.82158255577087) q[3];
sx q[3];
rz(4.46746054490144) q[3];
sx q[3];
rz(8.32938060759708) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.00051176548004) q[0];
sx q[0];
rz(3.77917537291581) q[0];
sx q[0];
rz(9.90432832240268) q[0];
rz(2.63748216629028) q[1];
sx q[1];
rz(5.13503006299073) q[1];
sx q[1];
rz(12.6246716737668) q[1];
cx q[1],q[0];
rz(3.82063889503479) q[0];
sx q[0];
rz(3.71473857958848) q[0];
sx q[0];
rz(8.34389374255344) q[0];
rz(2.17987895011902) q[2];
sx q[2];
rz(7.07138148148591) q[2];
sx q[2];
rz(11.9217154741208) q[2];
cx q[2],q[1];
rz(5.32467031478882) q[1];
sx q[1];
rz(0.564991863565989) q[1];
sx q[1];
rz(6.51135585307285) q[1];
rz(-1.06483864784241) q[3];
sx q[3];
rz(-1.33835968176787) q[3];
sx q[3];
rz(8.47876570223972) q[3];
cx q[3],q[2];
rz(-3.52799940109253) q[2];
sx q[2];
rz(5.10158804257447) q[2];
sx q[2];
rz(12.6036371946256) q[2];
rz(1.62183284759521) q[3];
sx q[3];
rz(3.59950629075105) q[3];
sx q[3];
rz(12.1849164724271) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.19504904747009) q[0];
sx q[0];
rz(6.75168076355989) q[0];
sx q[0];
rz(12.6329195260923) q[0];
rz(-1.52441239356995) q[1];
sx q[1];
rz(5.03972986538941) q[1];
sx q[1];
rz(10.1597195625226) q[1];
cx q[1],q[0];
rz(1.06828689575195) q[0];
sx q[0];
rz(3.77358445723588) q[0];
sx q[0];
rz(10.5962080717008) q[0];
rz(-3.37479782104492) q[2];
sx q[2];
rz(0.8542903979593) q[2];
sx q[2];
rz(8.66715673207446) q[2];
cx q[2],q[1];
rz(1.36147034168243) q[1];
sx q[1];
rz(4.35828688939149) q[1];
sx q[1];
rz(6.58169457911655) q[1];
rz(1.71646130084991) q[3];
sx q[3];
rz(3.99707177479798) q[3];
sx q[3];
rz(9.5551672488372) q[3];
cx q[3],q[2];
rz(-0.288960039615631) q[2];
sx q[2];
rz(4.70549753506715) q[2];
sx q[2];
rz(15.8166680097501) q[2];
rz(2.21478009223938) q[3];
sx q[3];
rz(5.31253233750398) q[3];
sx q[3];
rz(7.84699699877902) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.0842592716217) q[0];
sx q[0];
rz(3.77864608367021) q[0];
sx q[0];
rz(11.4755005597989) q[0];
rz(0.572516560554504) q[1];
sx q[1];
rz(5.09360173543031) q[1];
sx q[1];
rz(7.10749409197971) q[1];
cx q[1],q[0];
rz(-0.594127237796783) q[0];
sx q[0];
rz(3.59004205663735) q[0];
sx q[0];
rz(9.9889918923299) q[0];
rz(0.666070997714996) q[2];
sx q[2];
rz(2.38026294310624) q[2];
sx q[2];
rz(12.2728297471921) q[2];
cx q[2],q[1];
rz(2.90597724914551) q[1];
sx q[1];
rz(0.698254974680491) q[1];
sx q[1];
rz(8.78361270426914) q[1];
rz(-1.97639620304108) q[3];
sx q[3];
rz(4.31009438832337) q[3];
sx q[3];
rz(10.5989267587583) q[3];
cx q[3],q[2];
rz(-4.80222129821777) q[2];
sx q[2];
rz(2.16925600369508) q[2];
sx q[2];
rz(12.2751676797788) q[2];
rz(0.636998414993286) q[3];
sx q[3];
rz(4.8189344723993) q[3];
sx q[3];
rz(7.53566727637454) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.37851917743683) q[0];
sx q[0];
rz(5.99324169953401) q[0];
sx q[0];
rz(9.29372563063308) q[0];
rz(4.30849266052246) q[1];
sx q[1];
rz(8.44102254708345) q[1];
sx q[1];
rz(12.543750023834) q[1];
cx q[1],q[0];
rz(5.16377544403076) q[0];
sx q[0];
rz(7.50204578240449) q[0];
sx q[0];
rz(11.5763356447141) q[0];
rz(4.19289302825928) q[2];
sx q[2];
rz(6.03860154946382) q[2];
sx q[2];
rz(11.556116795532) q[2];
cx q[2],q[1];
rz(1.2056667804718) q[1];
sx q[1];
rz(2.49564108450944) q[1];
sx q[1];
rz(8.89606741666003) q[1];
rz(-0.698542416095734) q[3];
sx q[3];
rz(6.0333777983957) q[3];
sx q[3];
rz(8.16424701213046) q[3];
cx q[3],q[2];
rz(-1.87339460849762) q[2];
sx q[2];
rz(2.43842730124528) q[2];
sx q[2];
rz(13.6510996580045) q[2];
rz(1.14494132995605) q[3];
sx q[3];
rz(4.99651792843873) q[3];
sx q[3];
rz(8.40284607409641) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.55220830440521) q[0];
sx q[0];
rz(7.80199113686616) q[0];
sx q[0];
rz(12.0928485155027) q[0];
rz(1.22089684009552) q[1];
sx q[1];
rz(6.69561734993989) q[1];
sx q[1];
rz(10.6932829379956) q[1];
cx q[1],q[0];
rz(3.94777297973633) q[0];
sx q[0];
rz(2.07831934292848) q[0];
sx q[0];
rz(10.0270275235097) q[0];
rz(-2.24816846847534) q[2];
sx q[2];
rz(5.42240610917146) q[2];
sx q[2];
rz(10.8206786870877) q[2];
cx q[2],q[1];
rz(1.44226765632629) q[1];
sx q[1];
rz(5.21670833428437) q[1];
sx q[1];
rz(12.0858318567197) q[1];
rz(-0.176277324557304) q[3];
sx q[3];
rz(5.84120670159394) q[3];
sx q[3];
rz(9.79739636778041) q[3];
cx q[3],q[2];
rz(0.977416932582855) q[2];
sx q[2];
rz(4.7109319289499) q[2];
sx q[2];
rz(9.67119780778095) q[2];
rz(-0.942708015441895) q[3];
sx q[3];
rz(4.22753408749635) q[3];
sx q[3];
rz(10.1882087945859) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.81342685222626) q[0];
sx q[0];
rz(5.01281431515748) q[0];
sx q[0];
rz(9.48536993785902) q[0];
rz(4.64403676986694) q[1];
sx q[1];
rz(4.92012730439241) q[1];
sx q[1];
rz(7.87699720858737) q[1];
cx q[1],q[0];
rz(-0.12822026014328) q[0];
sx q[0];
rz(3.13760790008167) q[0];
sx q[0];
rz(11.2559746265332) q[0];
rz(-2.95340704917908) q[2];
sx q[2];
rz(7.28686872323091) q[2];
sx q[2];
rz(11.9724225759427) q[2];
cx q[2],q[1];
rz(-2.10515761375427) q[1];
sx q[1];
rz(2.16868248780305) q[1];
sx q[1];
rz(11.5330428838651) q[1];
rz(2.2864031791687) q[3];
sx q[3];
rz(7.50002876122529) q[3];
sx q[3];
rz(9.01096672414943) q[3];
cx q[3],q[2];
rz(-0.265552014112473) q[2];
sx q[2];
rz(4.06250110467012) q[2];
sx q[2];
rz(10.8820019721906) q[2];
rz(-0.11988104134798) q[3];
sx q[3];
rz(3.34617008467252) q[3];
sx q[3];
rz(11.5376887082975) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.847517132759094) q[0];
sx q[0];
rz(4.13599309523637) q[0];
sx q[0];
rz(10.5601308107297) q[0];
rz(-3.03822326660156) q[1];
sx q[1];
rz(4.87828746636445) q[1];
sx q[1];
rz(11.6187162160794) q[1];
cx q[1],q[0];
rz(-2.45548677444458) q[0];
sx q[0];
rz(2.49998328288133) q[0];
sx q[0];
rz(10.8836402654569) q[0];
rz(3.83002710342407) q[2];
sx q[2];
rz(-1.89700063864654) q[2];
sx q[2];
rz(12.1854793786924) q[2];
cx q[2],q[1];
rz(-2.74587416648865) q[1];
sx q[1];
rz(-1.3048790375418) q[1];
sx q[1];
rz(10.1201608538549) q[1];
rz(3.99870109558105) q[3];
sx q[3];
rz(7.16975560982759) q[3];
sx q[3];
rz(9.76585299371883) q[3];
cx q[3],q[2];
rz(-0.424429804086685) q[2];
sx q[2];
rz(3.76976767380769) q[2];
sx q[2];
rz(8.52044985293552) q[2];
rz(1.26339113712311) q[3];
sx q[3];
rz(5.04619398911531) q[3];
sx q[3];
rz(6.79967091082736) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.27102696895599) q[0];
sx q[0];
rz(3.89974698622758) q[0];
sx q[0];
rz(10.4101385235707) q[0];
rz(2.78997802734375) q[1];
sx q[1];
rz(5.32305136521394) q[1];
sx q[1];
rz(9.77164167761012) q[1];
cx q[1],q[0];
rz(3.20652198791504) q[0];
sx q[0];
rz(3.81492278178269) q[0];
sx q[0];
rz(9.0702002108018) q[0];
rz(-4.59984445571899) q[2];
sx q[2];
rz(4.26892438729341) q[2];
sx q[2];
rz(6.64960119723483) q[2];
cx q[2],q[1];
rz(-1.73694396018982) q[1];
sx q[1];
rz(0.882786663370677) q[1];
sx q[1];
rz(9.30956195890113) q[1];
rz(-0.366380304098129) q[3];
sx q[3];
rz(5.12702575524385) q[3];
sx q[3];
rz(12.2006807088773) q[3];
cx q[3],q[2];
rz(-0.88006192445755) q[2];
sx q[2];
rz(7.0021130164438) q[2];
sx q[2];
rz(9.28985617160007) q[2];
rz(3.01235818862915) q[3];
sx q[3];
rz(2.74154693086679) q[3];
sx q[3];
rz(8.38607392310306) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(4.73279476165771) q[0];
sx q[0];
rz(1.69065836270387) q[0];
sx q[0];
rz(9.10695699452564) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-1.14218950271606) q[1];
sx q[1];
rz(1.38154295285279) q[1];
sx q[1];
rz(12.3926071882169) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.560724198818207) q[2];
sx q[2];
rz(-0.816535321874074) q[2];
sx q[2];
rz(14.1406373739164) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.08561384677887) q[3];
sx q[3];
rz(4.8660079558664) q[3];
sx q[3];
rz(8.11786005496188) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
