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
rz(1.8062502) q[0];
sx q[0];
rz(-0.36646068) q[0];
sx q[0];
rz(-2.6824644) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(4.876457) q[1];
sx q[1];
rz(9.193037) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739959) q[0];
sx q[0];
rz(-1.2277506) q[0];
sx q[0];
rz(-1.6438707) q[0];
rz(1.628078) q[2];
sx q[2];
rz(-1.0250895) q[2];
sx q[2];
rz(-1.6451665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3331795) q[1];
sx q[1];
rz(-1.1929379) q[1];
sx q[1];
rz(-0.66587944) q[1];
rz(-pi) q[2];
rz(-0.11756331) q[3];
sx q[3];
rz(-2.8010578) q[3];
sx q[3];
rz(-1.3887608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0509402) q[2];
sx q[2];
rz(-0.53477627) q[2];
sx q[2];
rz(1.2789307) q[2];
rz(-0.88979641) q[3];
sx q[3];
rz(-1.5225478) q[3];
sx q[3];
rz(-0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5559674) q[0];
sx q[0];
rz(-0.90921679) q[0];
sx q[0];
rz(-2.2157748) q[0];
rz(-1.0649118) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(-0.085478641) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2538214) q[0];
sx q[0];
rz(-1.5675052) q[0];
sx q[0];
rz(-2.7517767) q[0];
rz(-pi) q[1];
rz(2.0485092) q[2];
sx q[2];
rz(-0.6685377) q[2];
sx q[2];
rz(-2.3195409) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1053727) q[1];
sx q[1];
rz(-0.58507996) q[1];
sx q[1];
rz(1.9963422) q[1];
x q[2];
rz(-1.1033589) q[3];
sx q[3];
rz(-0.93884838) q[3];
sx q[3];
rz(-1.8935209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.815879) q[2];
sx q[2];
rz(-1.6625983) q[2];
sx q[2];
rz(2.7549287) q[2];
rz(-1.9246842) q[3];
sx q[3];
rz(-2.7972126) q[3];
sx q[3];
rz(-0.2200505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.3306408) q[0];
sx q[0];
rz(-0.37136677) q[0];
sx q[0];
rz(-0.49705848) q[0];
rz(1.0423543) q[1];
sx q[1];
rz(-2.1940239) q[1];
sx q[1];
rz(0.29464468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.153107) q[0];
sx q[0];
rz(-1.9026347) q[0];
sx q[0];
rz(1.6436623) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9109086) q[2];
sx q[2];
rz(-2.4896224) q[2];
sx q[2];
rz(-0.99562746) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4528447) q[1];
sx q[1];
rz(-2.0032681) q[1];
sx q[1];
rz(2.2562863) q[1];
rz(-pi) q[2];
rz(-1.3623275) q[3];
sx q[3];
rz(-2.5304171) q[3];
sx q[3];
rz(-1.1217505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3543388) q[2];
sx q[2];
rz(-0.86091176) q[2];
sx q[2];
rz(1.5617237) q[2];
rz(2.0653557) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(-0.92897433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8266066) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(-2.9122747) q[0];
rz(-1.5062821) q[1];
sx q[1];
rz(-1.458026) q[1];
sx q[1];
rz(3.07952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45907606) q[0];
sx q[0];
rz(-0.78732765) q[0];
sx q[0];
rz(-0.83700257) q[0];
rz(0.27917316) q[2];
sx q[2];
rz(-1.9491157) q[2];
sx q[2];
rz(-0.73374464) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2885372) q[1];
sx q[1];
rz(-1.2934577) q[1];
sx q[1];
rz(-2.4278575) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7113879) q[3];
sx q[3];
rz(-2.9211126) q[3];
sx q[3];
rz(1.1663495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.092992358) q[2];
sx q[2];
rz(-0.91602641) q[2];
sx q[2];
rz(1.830706) q[2];
rz(0.038289573) q[3];
sx q[3];
rz(-1.3524651) q[3];
sx q[3];
rz(0.37461764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49486092) q[0];
sx q[0];
rz(-0.72135389) q[0];
sx q[0];
rz(0.89299655) q[0];
rz(-1.0294186) q[1];
sx q[1];
rz(-0.72975492) q[1];
sx q[1];
rz(-3.0457048) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5487089) q[0];
sx q[0];
rz(-1.0619231) q[0];
sx q[0];
rz(-2.6611317) q[0];
rz(-1.65972) q[2];
sx q[2];
rz(-1.6409573) q[2];
sx q[2];
rz(-1.905575) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3600625) q[1];
sx q[1];
rz(-1.9340098) q[1];
sx q[1];
rz(-0.59477721) q[1];
rz(-pi) q[2];
rz(-0.60378243) q[3];
sx q[3];
rz(-2.5386435) q[3];
sx q[3];
rz(-0.83151885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9868077) q[2];
sx q[2];
rz(-1.3900577) q[2];
sx q[2];
rz(2.7161157) q[2];
rz(0.42243877) q[3];
sx q[3];
rz(-1.1047624) q[3];
sx q[3];
rz(-0.5031684) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4130037) q[0];
sx q[0];
rz(-2.509403) q[0];
sx q[0];
rz(3.0244306) q[0];
rz(1.6475742) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(-2.07043) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2334918) q[0];
sx q[0];
rz(-1.0885493) q[0];
sx q[0];
rz(1.2511176) q[0];
rz(-pi) q[1];
rz(-2.6980459) q[2];
sx q[2];
rz(-2.5534592) q[2];
sx q[2];
rz(3.0556222) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2342959) q[1];
sx q[1];
rz(-1.7959372) q[1];
sx q[1];
rz(2.2683558) q[1];
x q[2];
rz(-0.46573205) q[3];
sx q[3];
rz(-2.7890786) q[3];
sx q[3];
rz(-0.62437144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5693207) q[2];
sx q[2];
rz(-0.57491493) q[2];
sx q[2];
rz(0.039483698) q[2];
rz(0.076400541) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(2.6881645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70306784) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(2.1912498) q[0];
rz(-1.9901265) q[1];
sx q[1];
rz(-1.6116424) q[1];
sx q[1];
rz(1.3075525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3111585) q[0];
sx q[0];
rz(-1.7896928) q[0];
sx q[0];
rz(1.9040742) q[0];
x q[1];
rz(-2.9588863) q[2];
sx q[2];
rz(-2.7596843) q[2];
sx q[2];
rz(0.71268247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4129909) q[1];
sx q[1];
rz(-2.6003693) q[1];
sx q[1];
rz(-1.9899998) q[1];
x q[2];
rz(-0.96486196) q[3];
sx q[3];
rz(-1.4601225) q[3];
sx q[3];
rz(-1.0718653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3211956) q[2];
sx q[2];
rz(-0.68009192) q[2];
sx q[2];
rz(2.5879228) q[2];
rz(-2.4456444) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11947908) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(-2.8346862) q[0];
rz(-1.2762997) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(-1.2409522) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81344714) q[0];
sx q[0];
rz(-1.4147583) q[0];
sx q[0];
rz(-2.0202083) q[0];
rz(-2.5503301) q[2];
sx q[2];
rz(-2.4545672) q[2];
sx q[2];
rz(2.9369773) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9556632) q[1];
sx q[1];
rz(-2.7765882) q[1];
sx q[1];
rz(-0.55162736) q[1];
x q[2];
rz(-1.9718902) q[3];
sx q[3];
rz(-2.1425284) q[3];
sx q[3];
rz(-0.16756646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79157311) q[2];
sx q[2];
rz(-0.80524033) q[2];
sx q[2];
rz(1.8512858) q[2];
rz(-0.6811412) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(-1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8660368) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(1.9679029) q[0];
rz(-2.5545919) q[1];
sx q[1];
rz(-0.78158164) q[1];
sx q[1];
rz(0.079708286) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3143334) q[0];
sx q[0];
rz(-1.200935) q[0];
sx q[0];
rz(0.24263675) q[0];
rz(-0.22144145) q[2];
sx q[2];
rz(-1.4613073) q[2];
sx q[2];
rz(1.2564645) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6708128) q[1];
sx q[1];
rz(-0.97320405) q[1];
sx q[1];
rz(2.2714991) q[1];
rz(-pi) q[2];
rz(1.2740785) q[3];
sx q[3];
rz(-2.35459) q[3];
sx q[3];
rz(-2.3001461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6012663) q[2];
sx q[2];
rz(-1.9129632) q[2];
sx q[2];
rz(1.8076753) q[2];
rz(-1.5506844) q[3];
sx q[3];
rz(-1.0100789) q[3];
sx q[3];
rz(-0.18868119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48225668) q[0];
sx q[0];
rz(-1.5784669) q[0];
sx q[0];
rz(1.851086) q[0];
rz(-3.105063) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(2.0972924) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8528906) q[0];
sx q[0];
rz(-1.1452617) q[0];
sx q[0];
rz(-2.8136926) q[0];
rz(0.26728435) q[2];
sx q[2];
rz(-1.2517831) q[2];
sx q[2];
rz(0.50179447) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86856213) q[1];
sx q[1];
rz(-2.2924463) q[1];
sx q[1];
rz(-0.20593404) q[1];
rz(-pi) q[2];
rz(1.6855816) q[3];
sx q[3];
rz(-1.8978531) q[3];
sx q[3];
rz(-0.17010526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9145987) q[2];
sx q[2];
rz(-1.031216) q[2];
sx q[2];
rz(1.3247789) q[2];
rz(-0.75657183) q[3];
sx q[3];
rz(-2.2533267) q[3];
sx q[3];
rz(0.85048401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4751547) q[0];
sx q[0];
rz(-1.7378687) q[0];
sx q[0];
rz(0.73874656) q[0];
rz(-1.4229763) q[1];
sx q[1];
rz(-1.9444793) q[1];
sx q[1];
rz(0.3107298) q[1];
rz(0.46427609) q[2];
sx q[2];
rz(-2.3515679) q[2];
sx q[2];
rz(-0.49754561) q[2];
rz(2.5632756) q[3];
sx q[3];
rz(-2.3925647) q[3];
sx q[3];
rz(-1.1306764) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
