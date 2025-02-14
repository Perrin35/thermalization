OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.13088432) q[0];
sx q[0];
rz(-0.67357981) q[0];
sx q[0];
rz(-0.82290617) q[0];
rz(-1.6910488) q[1];
sx q[1];
rz(-2.4658642) q[1];
sx q[1];
rz(2.9339209) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9638478) q[0];
sx q[0];
rz(-1.2262934) q[0];
sx q[0];
rz(2.3314407) q[0];
x q[1];
rz(1.0124341) q[2];
sx q[2];
rz(-1.2744546) q[2];
sx q[2];
rz(1.3989965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.339141) q[1];
sx q[1];
rz(-1.4467518) q[1];
sx q[1];
rz(-1.5088827) q[1];
rz(-pi) q[2];
rz(3.0700686) q[3];
sx q[3];
rz(-1.9404963) q[3];
sx q[3];
rz(-0.92310753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9660008) q[2];
sx q[2];
rz(-1.5660183) q[2];
sx q[2];
rz(-1.2292181) q[2];
rz(-1.843533) q[3];
sx q[3];
rz(-1.3763873) q[3];
sx q[3];
rz(-1.957533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9640279) q[0];
sx q[0];
rz(-0.25968817) q[0];
sx q[0];
rz(-0.92726707) q[0];
rz(2.941653) q[1];
sx q[1];
rz(-1.6688469) q[1];
sx q[1];
rz(-2.1520069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4662113) q[0];
sx q[0];
rz(-0.37801925) q[0];
sx q[0];
rz(-1.2751352) q[0];
x q[1];
rz(2.7071955) q[2];
sx q[2];
rz(-1.6100307) q[2];
sx q[2];
rz(-2.5920282) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.687508) q[1];
sx q[1];
rz(-1.6191779) q[1];
sx q[1];
rz(-1.494074) q[1];
x q[2];
rz(0.92156593) q[3];
sx q[3];
rz(-2.3147221) q[3];
sx q[3];
rz(-2.215883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2135311) q[2];
sx q[2];
rz(-1.4175043) q[2];
sx q[2];
rz(2.7878063) q[2];
rz(-0.95156041) q[3];
sx q[3];
rz(-0.74717251) q[3];
sx q[3];
rz(2.814754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78615171) q[0];
sx q[0];
rz(-2.1936301) q[0];
sx q[0];
rz(2.1308664) q[0];
rz(3.082869) q[1];
sx q[1];
rz(-1.6207691) q[1];
sx q[1];
rz(2.7322863) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5053913) q[0];
sx q[0];
rz(-2.276148) q[0];
sx q[0];
rz(-2.5558997) q[0];
rz(1.204973) q[2];
sx q[2];
rz(-2.7049148) q[2];
sx q[2];
rz(2.9819289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61534897) q[1];
sx q[1];
rz(-2.1110635) q[1];
sx q[1];
rz(-1.3413642) q[1];
rz(-pi) q[2];
rz(0.50893754) q[3];
sx q[3];
rz(-1.5526315) q[3];
sx q[3];
rz(0.048232676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0851486) q[2];
sx q[2];
rz(-1.8935545) q[2];
sx q[2];
rz(0.15963456) q[2];
rz(-1.2498445) q[3];
sx q[3];
rz(-1.4188473) q[3];
sx q[3];
rz(-2.0554525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98642629) q[0];
sx q[0];
rz(-2.7041589) q[0];
sx q[0];
rz(2.0470108) q[0];
rz(-1.1711586) q[1];
sx q[1];
rz(-1.1181701) q[1];
sx q[1];
rz(1.0135244) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60166925) q[0];
sx q[0];
rz(-0.058383103) q[0];
sx q[0];
rz(-2.4235382) q[0];
x q[1];
rz(-1.0763747) q[2];
sx q[2];
rz(-1.5740047) q[2];
sx q[2];
rz(-0.50960827) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5202433) q[1];
sx q[1];
rz(-1.3101556) q[1];
sx q[1];
rz(1.5288749) q[1];
rz(1.7988241) q[3];
sx q[3];
rz(-2.6963391) q[3];
sx q[3];
rz(-0.54603117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1363498) q[2];
sx q[2];
rz(-1.760027) q[2];
sx q[2];
rz(1.8278149) q[2];
rz(-0.93787307) q[3];
sx q[3];
rz(-1.2910941) q[3];
sx q[3];
rz(0.098793678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2336642) q[0];
sx q[0];
rz(-1.6133244) q[0];
sx q[0];
rz(2.0966356) q[0];
rz(-2.9772671) q[1];
sx q[1];
rz(-1.3427837) q[1];
sx q[1];
rz(2.9716861) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63498492) q[0];
sx q[0];
rz(-1.7568551) q[0];
sx q[0];
rz(-0.072441262) q[0];
rz(-pi) q[1];
rz(2.3218574) q[2];
sx q[2];
rz(-2.6728874) q[2];
sx q[2];
rz(0.021989487) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.53105866) q[1];
sx q[1];
rz(-1.8137168) q[1];
sx q[1];
rz(0.90621913) q[1];
x q[2];
rz(1.4395797) q[3];
sx q[3];
rz(-2.2926712) q[3];
sx q[3];
rz(0.038824507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5728411) q[2];
sx q[2];
rz(-2.4390287) q[2];
sx q[2];
rz(-1.0315726) q[2];
rz(2.0555563) q[3];
sx q[3];
rz(-2.3417754) q[3];
sx q[3];
rz(-1.731855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3373435) q[0];
sx q[0];
rz(-2.0954837) q[0];
sx q[0];
rz(-1.7403437) q[0];
rz(2.7119472) q[1];
sx q[1];
rz(-2.255217) q[1];
sx q[1];
rz(-1.8362129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84112924) q[0];
sx q[0];
rz(-2.3337737) q[0];
sx q[0];
rz(0.43231583) q[0];
rz(-1.13187) q[2];
sx q[2];
rz(-1.5666851) q[2];
sx q[2];
rz(-0.7059577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80524492) q[1];
sx q[1];
rz(-1.3224396) q[1];
sx q[1];
rz(-1.6026366) q[1];
rz(-pi) q[2];
rz(-2.9885653) q[3];
sx q[3];
rz(-1.8108441) q[3];
sx q[3];
rz(1.9757063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1244916) q[2];
sx q[2];
rz(-1.2064826) q[2];
sx q[2];
rz(0.98480946) q[2];
rz(1.5752327) q[3];
sx q[3];
rz(-1.5896348) q[3];
sx q[3];
rz(-2.8563833) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.259909) q[0];
sx q[0];
rz(-0.30170983) q[0];
sx q[0];
rz(1.786422) q[0];
rz(0.024070865) q[1];
sx q[1];
rz(-1.6203974) q[1];
sx q[1];
rz(-2.7640061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9474038) q[0];
sx q[0];
rz(-0.67127514) q[0];
sx q[0];
rz(2.8888974) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3578916) q[2];
sx q[2];
rz(-1.3981552) q[2];
sx q[2];
rz(1.6253234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.43289513) q[1];
sx q[1];
rz(-2.5835977) q[1];
sx q[1];
rz(-0.79496986) q[1];
rz(-pi) q[2];
rz(-1.1374279) q[3];
sx q[3];
rz(-1.876017) q[3];
sx q[3];
rz(-0.63852541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54921237) q[2];
sx q[2];
rz(-2.8040631) q[2];
sx q[2];
rz(-0.40360061) q[2];
rz(2.9523383) q[3];
sx q[3];
rz(-1.0523825) q[3];
sx q[3];
rz(2.6846867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52983487) q[0];
sx q[0];
rz(-0.54946041) q[0];
sx q[0];
rz(-2.8305565) q[0];
rz(-2.1850736) q[1];
sx q[1];
rz(-1.4776968) q[1];
sx q[1];
rz(2.8172353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2431902) q[0];
sx q[0];
rz(-1.0694188) q[0];
sx q[0];
rz(-1.7883975) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5115405) q[2];
sx q[2];
rz(-1.0612592) q[2];
sx q[2];
rz(1.6409724) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92128235) q[1];
sx q[1];
rz(-1.1571572) q[1];
sx q[1];
rz(-0.92343729) q[1];
x q[2];
rz(-2.2480184) q[3];
sx q[3];
rz(-1.8704318) q[3];
sx q[3];
rz(2.0183589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73584622) q[2];
sx q[2];
rz(-2.3509071) q[2];
sx q[2];
rz(0.22845593) q[2];
rz(-2.6089846) q[3];
sx q[3];
rz(-0.76627982) q[3];
sx q[3];
rz(0.84958357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5936977) q[0];
sx q[0];
rz(-2.6314647) q[0];
sx q[0];
rz(-1.8713895) q[0];
rz(2.5529329) q[1];
sx q[1];
rz(-1.207186) q[1];
sx q[1];
rz(-0.38280815) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7080806) q[0];
sx q[0];
rz(-1.545934) q[0];
sx q[0];
rz(1.3415501) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43395502) q[2];
sx q[2];
rz(-1.7979017) q[2];
sx q[2];
rz(-1.3004829) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36675699) q[1];
sx q[1];
rz(-2.2765571) q[1];
sx q[1];
rz(-1.4053759) q[1];
rz(-pi) q[2];
rz(1.0790014) q[3];
sx q[3];
rz(-1.331509) q[3];
sx q[3];
rz(-0.21904473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1560912) q[2];
sx q[2];
rz(-1.568855) q[2];
sx q[2];
rz(-0.4001948) q[2];
rz(-2.8943446) q[3];
sx q[3];
rz(-1.4618382) q[3];
sx q[3];
rz(2.7483773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8409214) q[0];
sx q[0];
rz(-1.3341757) q[0];
sx q[0];
rz(2.1260496) q[0];
rz(0.66889846) q[1];
sx q[1];
rz(-1.4543507) q[1];
sx q[1];
rz(1.6355754) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43011623) q[0];
sx q[0];
rz(-2.0680987) q[0];
sx q[0];
rz(-2.3946752) q[0];
rz(0.65176516) q[2];
sx q[2];
rz(-0.6547857) q[2];
sx q[2];
rz(-0.98818615) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4989711) q[1];
sx q[1];
rz(-1.6097415) q[1];
sx q[1];
rz(1.0907111) q[1];
x q[2];
rz(0.56699879) q[3];
sx q[3];
rz(-0.97494379) q[3];
sx q[3];
rz(-1.3064377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9754535) q[2];
sx q[2];
rz(-1.0756805) q[2];
sx q[2];
rz(-1.6157185) q[2];
rz(-0.29159355) q[3];
sx q[3];
rz(-0.44277954) q[3];
sx q[3];
rz(-0.35167882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6780554) q[0];
sx q[0];
rz(-2.1778477) q[0];
sx q[0];
rz(0.5974593) q[0];
rz(-2.8529104) q[1];
sx q[1];
rz(-0.79816993) q[1];
sx q[1];
rz(0.023963902) q[1];
rz(-0.409375) q[2];
sx q[2];
rz(-1.2843046) q[2];
sx q[2];
rz(1.6094111) q[2];
rz(0.47824974) q[3];
sx q[3];
rz(-2.6326055) q[3];
sx q[3];
rz(2.5647246) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
