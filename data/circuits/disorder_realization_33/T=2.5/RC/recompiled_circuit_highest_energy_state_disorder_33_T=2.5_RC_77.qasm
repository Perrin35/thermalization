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
rz(0.1641195) q[0];
sx q[0];
rz(-2.1174705) q[0];
sx q[0];
rz(0.48409387) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(4.7197309) q[1];
sx q[1];
rz(9.3698256) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23942102) q[0];
sx q[0];
rz(-2.4706984) q[0];
sx q[0];
rz(1.6210728) q[0];
rz(-0.058637549) q[2];
sx q[2];
rz(-1.8203041) q[2];
sx q[2];
rz(-0.80896689) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.41356) q[1];
sx q[1];
rz(-1.4384585) q[1];
sx q[1];
rz(1.669426) q[1];
rz(-pi) q[2];
rz(-2.0255768) q[3];
sx q[3];
rz(-0.6729799) q[3];
sx q[3];
rz(-2.2448886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.310828) q[2];
sx q[2];
rz(-0.97042933) q[2];
sx q[2];
rz(0.27466276) q[2];
rz(0.87485391) q[3];
sx q[3];
rz(-1.091205) q[3];
sx q[3];
rz(0.27354512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65983588) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(-2.7428395) q[0];
rz(-0.2440456) q[1];
sx q[1];
rz(-1.2363385) q[1];
sx q[1];
rz(2.0358548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2139521) q[0];
sx q[0];
rz(-1.1301148) q[0];
sx q[0];
rz(1.8233612) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5697487) q[2];
sx q[2];
rz(-0.99858054) q[2];
sx q[2];
rz(-0.68293152) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2011021) q[1];
sx q[1];
rz(-0.0041714287) q[1];
sx q[1];
rz(-0.68912403) q[1];
rz(-3.0209474) q[3];
sx q[3];
rz(-1.835726) q[3];
sx q[3];
rz(2.210833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70046052) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(2.0665118) q[2];
rz(2.4454146) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27580801) q[0];
sx q[0];
rz(-1.2514665) q[0];
sx q[0];
rz(0.21670093) q[0];
rz(-1.3801344) q[1];
sx q[1];
rz(-0.41788995) q[1];
sx q[1];
rz(2.5221141) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9443837) q[0];
sx q[0];
rz(-1.608662) q[0];
sx q[0];
rz(-0.029649563) q[0];
rz(0.046929788) q[2];
sx q[2];
rz(-1.6937733) q[2];
sx q[2];
rz(-0.3058946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8113831) q[1];
sx q[1];
rz(-0.37322497) q[1];
sx q[1];
rz(-2.1887145) q[1];
x q[2];
rz(0.61010078) q[3];
sx q[3];
rz(-1.0786295) q[3];
sx q[3];
rz(-0.51612332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8850024) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(0.090506434) q[2];
rz(-2.6835486) q[3];
sx q[3];
rz(-2.6543255) q[3];
sx q[3];
rz(2.0335782) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8103545) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(-0.93165398) q[0];
rz(-0.16904198) q[1];
sx q[1];
rz(-2.5773498) q[1];
sx q[1];
rz(-2.1021252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4783206) q[0];
sx q[0];
rz(-0.86000681) q[0];
sx q[0];
rz(-1.2511613) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3170856) q[2];
sx q[2];
rz(-1.7895964) q[2];
sx q[2];
rz(1.4267061) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89371496) q[1];
sx q[1];
rz(-1.7778686) q[1];
sx q[1];
rz(-0.061092579) q[1];
rz(-1.3576034) q[3];
sx q[3];
rz(-1.8490095) q[3];
sx q[3];
rz(1.1602064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6197551) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(2.6452765) q[2];
rz(-0.79658341) q[3];
sx q[3];
rz(-2.8556672) q[3];
sx q[3];
rz(-2.0485785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26688823) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(-2.1697178) q[0];
rz(3.019849) q[1];
sx q[1];
rz(-1.5309445) q[1];
sx q[1];
rz(1.3294539) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055453528) q[0];
sx q[0];
rz(-2.8290966) q[0];
sx q[0];
rz(-1.3700831) q[0];
x q[1];
rz(-0.41556032) q[2];
sx q[2];
rz(-2.0274799) q[2];
sx q[2];
rz(1.4765075) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7763627) q[1];
sx q[1];
rz(-0.91971469) q[1];
sx q[1];
rz(-1.0427703) q[1];
x q[2];
rz(-3.132344) q[3];
sx q[3];
rz(-1.6498749) q[3];
sx q[3];
rz(1.4656064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10151265) q[2];
sx q[2];
rz(-1.9209346) q[2];
sx q[2];
rz(0.15138781) q[2];
rz(-2.3769412) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5966655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.0720035) q[0];
sx q[0];
rz(-1.66865) q[0];
sx q[0];
rz(1.0714916) q[0];
rz(-0.53120652) q[1];
sx q[1];
rz(-1.4899645) q[1];
sx q[1];
rz(0.62613097) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3817978) q[0];
sx q[0];
rz(-0.017649895) q[0];
sx q[0];
rz(1.8438898) q[0];
rz(-pi) q[1];
rz(-1.8161723) q[2];
sx q[2];
rz(-2.1538556) q[2];
sx q[2];
rz(2.8193605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8257013) q[1];
sx q[1];
rz(-0.69676149) q[1];
sx q[1];
rz(-0.48721643) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49895309) q[3];
sx q[3];
rz(-1.9658943) q[3];
sx q[3];
rz(2.5753491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35817394) q[2];
sx q[2];
rz(-2.5219315) q[2];
sx q[2];
rz(1.0142856) q[2];
rz(1.8861534) q[3];
sx q[3];
rz(-1.3557326) q[3];
sx q[3];
rz(-2.0881418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805873) q[0];
sx q[0];
rz(-0.87562457) q[0];
sx q[0];
rz(2.2210806) q[0];
rz(3.0275184) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(-1.1309518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75926542) q[0];
sx q[0];
rz(-0.59774071) q[0];
sx q[0];
rz(3.1315342) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5860687) q[2];
sx q[2];
rz(-1.0706524) q[2];
sx q[2];
rz(1.9490567) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5798545) q[1];
sx q[1];
rz(-0.67719141) q[1];
sx q[1];
rz(1.3586292) q[1];
rz(-pi) q[2];
rz(1.6507963) q[3];
sx q[3];
rz(-2.0857852) q[3];
sx q[3];
rz(-1.966371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.34746927) q[2];
sx q[2];
rz(-0.64212126) q[2];
sx q[2];
rz(2.7327909) q[2];
rz(0.18320228) q[3];
sx q[3];
rz(-1.5902404) q[3];
sx q[3];
rz(-2.0028152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0278397) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(2.8507932) q[0];
rz(-2.193702) q[1];
sx q[1];
rz(-2.6746076) q[1];
sx q[1];
rz(-1.9416169) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5644218) q[0];
sx q[0];
rz(-2.0977712) q[0];
sx q[0];
rz(1.1061263) q[0];
rz(-pi) q[1];
rz(2.3340204) q[2];
sx q[2];
rz(-2.117273) q[2];
sx q[2];
rz(-1.228491) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6281451) q[1];
sx q[1];
rz(-1.2772868) q[1];
sx q[1];
rz(3.0931482) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.387595) q[3];
sx q[3];
rz(-1.4716513) q[3];
sx q[3];
rz(0.20618901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7802508) q[2];
sx q[2];
rz(-1.4486518) q[2];
sx q[2];
rz(-1.7928436) q[2];
rz(1.5032984) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(2.9564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0272442) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(2.0674904) q[0];
rz(-0.4298003) q[1];
sx q[1];
rz(-1.0440412) q[1];
sx q[1];
rz(0.46357402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060424711) q[0];
sx q[0];
rz(-1.7588108) q[0];
sx q[0];
rz(0.096676143) q[0];
rz(1.3691291) q[2];
sx q[2];
rz(-1.3140162) q[2];
sx q[2];
rz(-0.078291206) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1660171) q[1];
sx q[1];
rz(-2.3834627) q[1];
sx q[1];
rz(-1.0075955) q[1];
rz(-1.4809247) q[3];
sx q[3];
rz(-2.1276132) q[3];
sx q[3];
rz(-1.9720322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8672436) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(-2.3331433) q[2];
rz(0.48458734) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3633858) q[0];
sx q[0];
rz(-2.6860542) q[0];
sx q[0];
rz(2.0090012) q[0];
rz(3.1032108) q[1];
sx q[1];
rz(-1.3382341) q[1];
sx q[1];
rz(-2.8448232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1365741) q[0];
sx q[0];
rz(-1.8474192) q[0];
sx q[0];
rz(1.0199976) q[0];
x q[1];
rz(-0.44906868) q[2];
sx q[2];
rz(-1.1658323) q[2];
sx q[2];
rz(1.2580037) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3546875) q[1];
sx q[1];
rz(-0.2345492) q[1];
sx q[1];
rz(-2.812318) q[1];
rz(-pi) q[2];
rz(2.040801) q[3];
sx q[3];
rz(-0.74277011) q[3];
sx q[3];
rz(-2.7239885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9639637) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(-0.56619823) q[2];
rz(1.1244134) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(-1.1894777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-0.88846702) q[0];
sx q[0];
rz(-2.4430226) q[0];
sx q[0];
rz(-3.0342614) q[0];
rz(0.26168564) q[1];
sx q[1];
rz(-1.2005922) q[1];
sx q[1];
rz(-2.190879) q[1];
rz(-3.059584) q[2];
sx q[2];
rz(-1.5463943) q[2];
sx q[2];
rz(0.80873185) q[2];
rz(-0.64601267) q[3];
sx q[3];
rz(-0.96886841) q[3];
sx q[3];
rz(-2.3652707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
