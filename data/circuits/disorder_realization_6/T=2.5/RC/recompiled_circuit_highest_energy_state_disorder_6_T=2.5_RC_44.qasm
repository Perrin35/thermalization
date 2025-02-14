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
rz(-1.8974798) q[0];
sx q[0];
rz(-0.31960684) q[0];
sx q[0];
rz(2.6407916) q[0];
rz(0.29419857) q[1];
sx q[1];
rz(-2.3994212) q[1];
sx q[1];
rz(0.54744005) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0078024) q[0];
sx q[0];
rz(-0.34388516) q[0];
sx q[0];
rz(-2.6965428) q[0];
rz(-pi) q[1];
rz(-2.6260761) q[2];
sx q[2];
rz(-0.52243865) q[2];
sx q[2];
rz(-1.9439657) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3708122) q[1];
sx q[1];
rz(-0.51518422) q[1];
sx q[1];
rz(-2.9574802) q[1];
x q[2];
rz(-2.7846401) q[3];
sx q[3];
rz(-1.0610749) q[3];
sx q[3];
rz(-1.3689643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9886292) q[2];
sx q[2];
rz(-2.2965778) q[2];
sx q[2];
rz(-0.60626924) q[2];
rz(1.2074977) q[3];
sx q[3];
rz(-1.2946125) q[3];
sx q[3];
rz(-1.5906364) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9015053) q[0];
sx q[0];
rz(-1.5272239) q[0];
sx q[0];
rz(2.4497581) q[0];
rz(-1.2884619) q[1];
sx q[1];
rz(-1.144578) q[1];
sx q[1];
rz(0.28786927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26507227) q[0];
sx q[0];
rz(-1.2848043) q[0];
sx q[0];
rz(3.0080481) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1580196) q[2];
sx q[2];
rz(-0.97146875) q[2];
sx q[2];
rz(2.9119208) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6791742) q[1];
sx q[1];
rz(-0.44826213) q[1];
sx q[1];
rz(-3.03469) q[1];
rz(1.6784662) q[3];
sx q[3];
rz(-1.7070908) q[3];
sx q[3];
rz(-1.8267269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79491478) q[2];
sx q[2];
rz(-1.6592462) q[2];
sx q[2];
rz(-0.97290009) q[2];
rz(-1.3341303) q[3];
sx q[3];
rz(-0.79136807) q[3];
sx q[3];
rz(-2.5202675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7767104) q[0];
sx q[0];
rz(-2.8948687) q[0];
sx q[0];
rz(0.2488939) q[0];
rz(-1.6913255) q[1];
sx q[1];
rz(-1.9906809) q[1];
sx q[1];
rz(2.2379564) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.929707) q[0];
sx q[0];
rz(-0.60962661) q[0];
sx q[0];
rz(0.044187688) q[0];
rz(-1.2799706) q[2];
sx q[2];
rz(-2.7827459) q[2];
sx q[2];
rz(-2.1339231) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.94374471) q[1];
sx q[1];
rz(-1.1863803) q[1];
sx q[1];
rz(2.3626087) q[1];
rz(-pi) q[2];
rz(0.55428688) q[3];
sx q[3];
rz(-0.6817556) q[3];
sx q[3];
rz(-1.1415014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.16477188) q[2];
sx q[2];
rz(-2.4103319) q[2];
sx q[2];
rz(-2.8110647) q[2];
rz(0.22322379) q[3];
sx q[3];
rz(-1.1130755) q[3];
sx q[3];
rz(0.37842512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44499236) q[0];
sx q[0];
rz(-1.3560504) q[0];
sx q[0];
rz(-0.14183216) q[0];
rz(-3.0524571) q[1];
sx q[1];
rz(-0.24239692) q[1];
sx q[1];
rz(2.1884122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2640799) q[0];
sx q[0];
rz(-0.82897128) q[0];
sx q[0];
rz(0.52998752) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3802657) q[2];
sx q[2];
rz(-1.6820246) q[2];
sx q[2];
rz(-2.5750772) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64884463) q[1];
sx q[1];
rz(-0.78300965) q[1];
sx q[1];
rz(-2.3451508) q[1];
rz(-pi) q[2];
rz(1.7715376) q[3];
sx q[3];
rz(-2.0588518) q[3];
sx q[3];
rz(2.1303915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7715093) q[2];
sx q[2];
rz(-0.13088317) q[2];
sx q[2];
rz(-2.6123602) q[2];
rz(-2.0970763) q[3];
sx q[3];
rz(-1.7525201) q[3];
sx q[3];
rz(0.15531003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.9411377) q[0];
sx q[0];
rz(-1.8998242) q[0];
sx q[0];
rz(0.45503765) q[0];
rz(3.1269238) q[1];
sx q[1];
rz(-0.55611062) q[1];
sx q[1];
rz(-0.98247772) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9519314) q[0];
sx q[0];
rz(-1.8704544) q[0];
sx q[0];
rz(-1.4656187) q[0];
rz(-pi) q[1];
rz(0.87329726) q[2];
sx q[2];
rz(-1.1050537) q[2];
sx q[2];
rz(-2.1312775) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8807672) q[1];
sx q[1];
rz(-1.8463355) q[1];
sx q[1];
rz(-1.4620305) q[1];
rz(-pi) q[2];
rz(-2.4508844) q[3];
sx q[3];
rz(-2.7793607) q[3];
sx q[3];
rz(-1.9286523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5386397) q[2];
sx q[2];
rz(-2.2621138) q[2];
sx q[2];
rz(-2.5243916) q[2];
rz(-1.7547102) q[3];
sx q[3];
rz(-0.91104561) q[3];
sx q[3];
rz(0.174218) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0571601) q[0];
sx q[0];
rz(-2.3800157) q[0];
sx q[0];
rz(0.97849751) q[0];
rz(1.017336) q[1];
sx q[1];
rz(-2.8637736) q[1];
sx q[1];
rz(0.4496347) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7788197) q[0];
sx q[0];
rz(-1.3032493) q[0];
sx q[0];
rz(0.31407251) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80763517) q[2];
sx q[2];
rz(-1.908632) q[2];
sx q[2];
rz(2.5011245) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79053959) q[1];
sx q[1];
rz(-1.0996523) q[1];
sx q[1];
rz(-3.0912249) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4586304) q[3];
sx q[3];
rz(-1.0691081) q[3];
sx q[3];
rz(1.7428118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6189239) q[2];
sx q[2];
rz(-1.7974682) q[2];
sx q[2];
rz(1.7632923) q[2];
rz(-1.1528692) q[3];
sx q[3];
rz(-1.1533573) q[3];
sx q[3];
rz(2.8809179) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5744837) q[0];
sx q[0];
rz(-0.672988) q[0];
sx q[0];
rz(-2.8048977) q[0];
rz(-2.4163272) q[1];
sx q[1];
rz(-1.5870353) q[1];
sx q[1];
rz(-0.38549647) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8367675) q[0];
sx q[0];
rz(-2.4186717) q[0];
sx q[0];
rz(-0.54534955) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0573559) q[2];
sx q[2];
rz(-0.2845736) q[2];
sx q[2];
rz(-1.3651266) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.25076097) q[1];
sx q[1];
rz(-2.4422993) q[1];
sx q[1];
rz(-2.0307402) q[1];
rz(-pi) q[2];
x q[2];
rz(2.518744) q[3];
sx q[3];
rz(-1.2174574) q[3];
sx q[3];
rz(1.6652157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54726279) q[2];
sx q[2];
rz(-1.9219834) q[2];
sx q[2];
rz(0.29898137) q[2];
rz(-0.79214054) q[3];
sx q[3];
rz(-0.39528379) q[3];
sx q[3];
rz(-1.7637174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2250273) q[0];
sx q[0];
rz(-1.7074317) q[0];
sx q[0];
rz(2.3634971) q[0];
rz(1.7447225) q[1];
sx q[1];
rz(-0.94051802) q[1];
sx q[1];
rz(-3.0527589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4230115) q[0];
sx q[0];
rz(-1.4664343) q[0];
sx q[0];
rz(2.7947438) q[0];
x q[1];
rz(-0.85529144) q[2];
sx q[2];
rz(-1.119985) q[2];
sx q[2];
rz(-0.26329783) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.962304) q[1];
sx q[1];
rz(-1.3109164) q[1];
sx q[1];
rz(1.13729) q[1];
rz(1.9384814) q[3];
sx q[3];
rz(-1.5993777) q[3];
sx q[3];
rz(2.1052484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54625154) q[2];
sx q[2];
rz(-2.310414) q[2];
sx q[2];
rz(2.0383932) q[2];
rz(-2.5698419) q[3];
sx q[3];
rz(-2.5751556) q[3];
sx q[3];
rz(-1.6167238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.53200805) q[0];
sx q[0];
rz(-3.0201865) q[0];
sx q[0];
rz(-0.65015656) q[0];
rz(-1.0981759) q[1];
sx q[1];
rz(-1.2547341) q[1];
sx q[1];
rz(3.0756899) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0528735) q[0];
sx q[0];
rz(-2.2744226) q[0];
sx q[0];
rz(2.5077453) q[0];
rz(-pi) q[1];
rz(1.5123821) q[2];
sx q[2];
rz(-1.3810535) q[2];
sx q[2];
rz(-2.3924215) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53532584) q[1];
sx q[1];
rz(-0.78731489) q[1];
sx q[1];
rz(-1.7892862) q[1];
x q[2];
rz(2.8206943) q[3];
sx q[3];
rz(-0.68701744) q[3];
sx q[3];
rz(1.4167757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99964511) q[2];
sx q[2];
rz(-1.1885175) q[2];
sx q[2];
rz(-0.56335062) q[2];
rz(2.3422286) q[3];
sx q[3];
rz(-2.1269709) q[3];
sx q[3];
rz(-2.8857901) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.946741) q[0];
sx q[0];
rz(-0.020337157) q[0];
sx q[0];
rz(2.643423) q[0];
rz(0.26214552) q[1];
sx q[1];
rz(-1.405895) q[1];
sx q[1];
rz(-1.3462876) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58459613) q[0];
sx q[0];
rz(-0.89476953) q[0];
sx q[0];
rz(-2.1732794) q[0];
rz(-3.1392241) q[2];
sx q[2];
rz(-1.5813229) q[2];
sx q[2];
rz(1.0536556) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.4075824) q[1];
sx q[1];
rz(-0.90208131) q[1];
sx q[1];
rz(-2.0120502) q[1];
x q[2];
rz(0.22110053) q[3];
sx q[3];
rz(-2.3490864) q[3];
sx q[3];
rz(1.4421876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4761618) q[2];
sx q[2];
rz(-0.68887812) q[2];
sx q[2];
rz(-0.8821744) q[2];
rz(2.4836508) q[3];
sx q[3];
rz(-2.0177757) q[3];
sx q[3];
rz(-0.9995802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5049725) q[0];
sx q[0];
rz(-1.5047147) q[0];
sx q[0];
rz(2.104105) q[0];
rz(2.6344015) q[1];
sx q[1];
rz(-2.3585944) q[1];
sx q[1];
rz(2.0814887) q[1];
rz(1.7569594) q[2];
sx q[2];
rz(-1.8590448) q[2];
sx q[2];
rz(0.94696905) q[2];
rz(-2.2118123) q[3];
sx q[3];
rz(-2.2140183) q[3];
sx q[3];
rz(2.1341013) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
