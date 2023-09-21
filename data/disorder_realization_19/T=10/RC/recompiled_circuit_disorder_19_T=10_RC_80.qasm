OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0269545) q[0];
sx q[0];
rz(-1.6898328) q[0];
sx q[0];
rz(0.64557689) q[0];
rz(-2.7627856) q[1];
sx q[1];
rz(-1.768755) q[1];
sx q[1];
rz(-1.6436613) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0092702) q[0];
sx q[0];
rz(-1.9267123) q[0];
sx q[0];
rz(0.15412553) q[0];
x q[1];
rz(3.0718578) q[2];
sx q[2];
rz(-2.2001534) q[2];
sx q[2];
rz(0.33915181) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1995391) q[1];
sx q[1];
rz(-2.374211) q[1];
sx q[1];
rz(-2.8295423) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1708158) q[3];
sx q[3];
rz(-2.1809289) q[3];
sx q[3];
rz(-2.4657472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(0.34040889) q[2];
rz(0.83299625) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.67265636) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-2.5966068) q[0];
rz(-0.90822059) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(-2.3166336) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9075515) q[0];
sx q[0];
rz(-1.7579161) q[0];
sx q[0];
rz(-1.9812816) q[0];
x q[1];
rz(1.5603793) q[2];
sx q[2];
rz(-1.9736819) q[2];
sx q[2];
rz(1.8889697) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0032012) q[1];
sx q[1];
rz(-1.0111286) q[1];
sx q[1];
rz(1.707934) q[1];
rz(1.2689781) q[3];
sx q[3];
rz(-0.58773731) q[3];
sx q[3];
rz(2.1835818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5543582) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(2.9651802) q[2];
rz(0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(-0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762887) q[0];
sx q[0];
rz(-0.58350199) q[0];
sx q[0];
rz(-2.6718455) q[0];
rz(1.6167971) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(-0.038539561) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32917133) q[0];
sx q[0];
rz(-3.1338131) q[0];
sx q[0];
rz(-1.6532941) q[0];
rz(0.58263393) q[2];
sx q[2];
rz(-2.8672672) q[2];
sx q[2];
rz(-2.3290616) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5650428) q[1];
sx q[1];
rz(-1.0817263) q[1];
sx q[1];
rz(-3.0719041) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79629691) q[3];
sx q[3];
rz(-2.3781804) q[3];
sx q[3];
rz(0.89087668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(0.91252404) q[2];
rz(1.3085261) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(-1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7581166) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(2.7918949) q[0];
rz(1.8967459) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53699025) q[0];
sx q[0];
rz(-1.6911117) q[0];
sx q[0];
rz(-0.093219482) q[0];
x q[1];
rz(-2.1148557) q[2];
sx q[2];
rz(-1.9907021) q[2];
sx q[2];
rz(1.3865711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8198204) q[1];
sx q[1];
rz(-1.1204549) q[1];
sx q[1];
rz(2.8798073) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5870729) q[3];
sx q[3];
rz(-2.2098594) q[3];
sx q[3];
rz(1.0234327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23094709) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(1.267743) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0347663) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(0.1396133) q[0];
rz(-1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(2.7635014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17079167) q[0];
sx q[0];
rz(-2.1149201) q[0];
sx q[0];
rz(-2.8118954) q[0];
x q[1];
rz(0.80412229) q[2];
sx q[2];
rz(-2.1466549) q[2];
sx q[2];
rz(-0.59876195) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2851583) q[1];
sx q[1];
rz(-0.57796961) q[1];
sx q[1];
rz(-2.2846089) q[1];
rz(0.46697163) q[3];
sx q[3];
rz(-1.3936371) q[3];
sx q[3];
rz(0.60764473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87171626) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(-0.88796973) q[2];
rz(2.1652083) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.75509214) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(-0.37762541) q[0];
rz(2.8185484) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(0.91517085) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5765502) q[0];
sx q[0];
rz(-0.33320198) q[0];
sx q[0];
rz(0.25175005) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2890131) q[2];
sx q[2];
rz(-2.4396509) q[2];
sx q[2];
rz(-2.728087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.23657654) q[1];
sx q[1];
rz(-1.6394776) q[1];
sx q[1];
rz(-2.828039) q[1];
x q[2];
rz(3.0056474) q[3];
sx q[3];
rz(-0.69909401) q[3];
sx q[3];
rz(-1.5232777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53720981) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(2.3510695) q[2];
rz(2.8516155) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73615605) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(0.36703584) q[0];
rz(-1.908318) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.4253915) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050497748) q[0];
sx q[0];
rz(-2.0031345) q[0];
sx q[0];
rz(-1.3652486) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0076214) q[2];
sx q[2];
rz(-2.6920762) q[2];
sx q[2];
rz(-1.1531032) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5126309) q[1];
sx q[1];
rz(-2.1753864) q[1];
sx q[1];
rz(2.0357893) q[1];
x q[2];
rz(3.1281934) q[3];
sx q[3];
rz(-2.074476) q[3];
sx q[3];
rz(-2.3434533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1620862) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(1.210775) q[2];
rz(-0.0018421729) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84412557) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-0.076518245) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(-0.75235596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6589688) q[0];
sx q[0];
rz(-2.0730744) q[0];
sx q[0];
rz(1.9320022) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1444507) q[2];
sx q[2];
rz(-2.2144197) q[2];
sx q[2];
rz(-2.2709536) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68554316) q[1];
sx q[1];
rz(-1.4577663) q[1];
sx q[1];
rz(-3.0526524) q[1];
rz(-pi) q[2];
rz(2.5520677) q[3];
sx q[3];
rz(-1.3988964) q[3];
sx q[3];
rz(0.37026438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(-0.00014649815) q[2];
rz(-2.032062) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(-2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14934854) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(-2.1355656) q[0];
rz(-0.26793119) q[1];
sx q[1];
rz(-1.8012828) q[1];
sx q[1];
rz(-1.3148274) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7816759) q[0];
sx q[0];
rz(-1.8209551) q[0];
sx q[0];
rz(0.71235384) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7039653) q[2];
sx q[2];
rz(-1.2097934) q[2];
sx q[2];
rz(-1.5683057) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4115574) q[1];
sx q[1];
rz(-2.3282603) q[1];
sx q[1];
rz(0.98585339) q[1];
rz(2.3141765) q[3];
sx q[3];
rz(-1.5762868) q[3];
sx q[3];
rz(-2.1237802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.43977794) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(2.8961199) q[2];
rz(-2.7108575) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7219287) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(-0.39189664) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6673198) q[0];
sx q[0];
rz(-2.5840238) q[0];
sx q[0];
rz(1.5865109) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13550831) q[2];
sx q[2];
rz(-2.2221609) q[2];
sx q[2];
rz(-1.1609921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74589409) q[1];
sx q[1];
rz(-2.3682526) q[1];
sx q[1];
rz(1.8383154) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8896905) q[3];
sx q[3];
rz(-1.0341757) q[3];
sx q[3];
rz(1.5434138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(-1.1575451) q[2];
rz(0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3863603) q[0];
sx q[0];
rz(-1.3437143) q[0];
sx q[0];
rz(1.2585826) q[0];
rz(1.339284) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(-2.3618868) q[2];
sx q[2];
rz(-1.4479326) q[2];
sx q[2];
rz(-0.89276199) q[2];
rz(-1.0212392) q[3];
sx q[3];
rz(-1.869535) q[3];
sx q[3];
rz(0.46365999) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
